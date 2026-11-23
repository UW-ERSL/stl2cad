"""Mesh segmentation into analytic primitives (efficient-RANSAC style, fixed seed),
followed by design-intent snapping and optional symmetric orbit refits.
"""

import os
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .primitives import (ALGEBRAIC, COMPLEXITY, N_SAMPLES, Prim, fit, fit_extrusion, fit_revolution, from_samples,
                         match_score, perp_basis, torus_init)


class Segmenter:
    def __init__(self, mesh, tol, kinds=("plane", "cylinder", "cone", "sphere", "torus"),
                 allow_extrusion=False, allow_facets=True, angle_tol=25.0, min_region_frac=2e-4,
                 samples_per_kind=80, seed=0, log=print, concave_only=False, r_max=None):
        self.m = mesh
        self.tol = tol
        self.kinds = [k for k in kinds if k in ALGEBRAIC]
        self.allow_extrusion = allow_extrusion
        self.allow_facets = allow_facets
        self.cos_ang = np.cos(np.radians(angle_tol))
        self.min_area = min_region_frac * mesh.area
        self.samples = samples_per_kind
        self.rng = np.random.default_rng(seed)
        self.log = log
        self.concave_only = concave_only
        self.r_max = r_max
        self.diag = float(np.linalg.norm(mesh.extents))
        F = len(mesh.faces)
        adj = mesh.face_adjacency
        self.adj = coo_matrix((np.ones(len(adj)), (adj[:, 0], adj[:, 1])), shape=(F, F)).tocsr()
        self.adj = self.adj + self.adj.T
        self.fc = mesh.triangles_center
        self.fn = mesh.face_normals
        self.fa = mesh.area_faces
        self.tree = cKDTree(self.fc)
        T = mesh.triangles
        self.edge_len = np.max(np.linalg.norm(T - np.roll(T, 1, axis=1), axis=2), axis=1)
        sm = mesh.face_adjacency_angles < np.radians(25)
        A = coo_matrix((np.ones(sm.sum()), (adj[sm, 0], adj[sm, 1])), shape=(F, F)).tocsr()
        self.adj_smooth = A + A.T

    @staticmethod
    def curv_radius(pr, X):
        k, p = pr.kind, pr.p
        if k == "plane":
            return np.full(len(X), np.inf)
        if k in ("cylinder", "sphere"):
            return np.full(len(X), p["r"])
        if k == "torus":
            return np.full(len(X), min(p["r"], abs(p["R"] - p["r"]) + 1e-12))
        if k == "cone":
            w = X - p["v"]
            h = w @ p["a"]
            rho = np.linalg.norm(w - np.outer(h, p["a"]), axis=1)
            return np.maximum(rho / np.cos(p["t"]), 1e-12)
        return np.full(len(X), np.inf)

    # ---- region evaluation ---------------------------------------------------
    def inliers(self, pr, cand_faces, scale=1.0):
        tol = self.tol * scale
        F = self.m.faces[cand_faces]
        # distances at the three corners of each candidate triangle (no np.unique: it dominated
        # the cost; shared vertices are evaluated more than once, with identical results)
        dv = np.abs(pr.distance(self.m.vertices[F.reshape(-1)])).reshape(-1, 3)
        ok = dv.max(1) < tol
        if not ok.any():
            return cand_faces[:0]
        sel = cand_faces[ok]
        # centroid test, allowing the chord sag of a correctly tessellated curved surface
        sag = self.edge_len[sel] ** 2 / (8 * self.curv_radius(pr, self.fc[sel]))
        ok = np.abs(pr.distance(self.fc[sel])) < tol + 1.05 * sag
        sel = sel[ok]
        if not len(sel):
            return sel
        nrm = pr.normal(self.fc[sel])
        ok2 = np.abs(np.einsum("ij,ij->i", nrm, self.fn[sel])) > self.cos_ang
        return sel[ok2]

    def largest_component(self, faces):
        if len(faces) == 0:
            return faces
        sub = self.adj[faces][:, faces]
        n, lab = connected_components(sub, directed=False)
        if n == 1:
            return faces
        areas = np.bincount(lab, weights=self.fa[faces])
        return faces[lab == np.argmax(areas)]

    def refine(self, pr, remaining_faces, iters=3):
        faces = self.largest_component(self.inliers(pr, remaining_faces))
        for _ in range(iters):
            if len(faces) < 2:
                return None, faces
            X = self.m.vertices[np.unique(self.m.faces[faces])]
            if pr.kind == "torus":
                q = fit("torus", X, init=pr)
            else:
                q = fit(pr.kind, X, init=pr)
            if q is None or not np.isfinite(q.rms):
                return None, faces
            new = self.largest_component(self.inliers(q, remaining_faces))
            if self.fa[new].sum() < 0.9 * self.fa[faces].sum():
                break
            pr, faces = q, new
        return pr, faces

    # ---- main loop -----------------------------------------------------------
    def run(self, max_prims=400, patience=0):
        mesh = self.m
        remaining = np.ones(len(mesh.faces), bool)
        prims = []
        ransac_kinds = [k for k in self.kinds if k in N_SAMPLES]
        local_r = 0.06 * self.diag
        empty_rounds = 0
        while remaining.any() and len(prims) < max_prims:
            rem = np.flatnonzero(remaining)
            if self.fa[rem].sum() < self.min_area:
                break
            prob = self.fa[rem] / self.fa[rem].sum()
            best = None
            for kind in ransac_kinds:
                k = N_SAMPLES[kind]
                for _ in range(self.samples if kind != "plane" else max(20, self.samples // 2)):
                    i = rem[self.rng.choice(len(rem), p=prob)]
                    idx = [i]
                    if k > 1:
                        nb = self.tree.query_ball_point(self.fc[i], local_r)
                        nb = [j for j in nb if remaining[j] and j != i]
                        if len(nb) < k - 1:
                            continue
                        idx += list(self.rng.choice(nb, k - 1, replace=False))
                    pr = from_samples(kind, self.fc[idx], self.fn[idx], self.diag)
                    if pr is None:
                        continue
                    if self.concave_only and kind == "cylinder":
                        # hole: the axis lies on the side the outward normal points to
                        if (pr.p["c"] - self.fc[i]) @ self.fn[i] <= 0 or (self.r_max and pr.p["r"] > self.r_max):
                            continue
                    sel = self.largest_component(self.inliers(pr, rem))
                    score = self.fa[sel].sum() * (1 - 0.03 * COMPLEXITY[kind])
                    if best is None or score > best[0]:
                        best = (score, pr, sel)
            if "torus" in self.kinds:
                nt = max(10, self.samples // 3)
                for _ in range(nt):
                    i = rem[self.rng.choice(len(rem), p=prob)]
                    nb = [j for j in self.tree.query_ball_point(self.fc[i], local_r) if remaining[j]]
                    if len(nb) < 10:
                        continue
                    nb = np.array(nb)[np.argsort(np.linalg.norm(self.fc[nb] - self.fc[i], axis=1))[:24]]
                    pr = torus_init(self.fc[nb], self.fn[nb])
                    if pr is None or not (0 < pr.p["r"] < 0.5 * self.diag) or pr.p["R"] > self.diag:
                        continue
                    sel = self.largest_component(self.inliers(pr, rem, scale=2.0))
                    score = self.fa[sel].sum() * (1 - 0.03 * COMPLEXITY["torus"])
                    if best is None or score > best[0]:
                        best = (score, pr, sel)
            if best is None or best[0] < self.min_area:
                empty_rounds += 1
                if empty_rounds > patience:
                    break
                continue
            empty_rounds = 0
            pr, faces = self.refine(best[1], rem)
            if pr is None or self.fa[faces].sum() < self.min_area:
                # accept the unrefined candidate region as unexplained, so the loop progresses
                bad = best[2]
                remaining[bad] = False
                self._unexplained = getattr(self, "_unexplained", []) + [bad]
                continue
            X = mesh.vertices[np.unique(mesh.faces[faces])]
            pr.faces, pr.pts, pr.area = faces, X, float(self.fa[faces].sum())
            pr.rms = float(np.sqrt(np.mean(pr.distance(X) ** 2)))
            prims.append(pr)
            remaining[faces] = False
        leftovers = np.flatnonzero(remaining)
        extra = getattr(self, "_unexplained", [])
        if extra:
            leftovers = np.unique(np.r_[leftovers, np.concatenate(extra)])
        prims, leftovers = self.absorb(prims, leftovers)
        prims = self.compete(prims)
        prims, leftovers = self.validate(prims, leftovers)
        if self.allow_extrusion:
            prims, leftovers = self.merge_extrusions(prims, leftovers)
        prims += self.explain_leftovers(leftovers)
        return prims

    def validate(self, prims, leftovers, scale=2.5, rounds=3):
        """Expel faces farther than scale*tol from their surface (drift after merging), refit;
        expelled faces become leftovers."""
        expelled = []
        for _ in range(rounds):
            changed = False
            for i, p in enumerate(prims):
                if p.kind == "facets" or not len(p.faces):
                    continue
                ok = np.isin(p.faces, self.inliers(p, p.faces, scale=scale))
                if ok.all():
                    continue
                changed = True
                expelled.append(p.faces[~ok])
                faces = p.faces[ok]
                if len(faces) < 2:
                    p.faces = faces
                    continue
                X = self.m.vertices[np.unique(self.m.faces[faces])]
                q = fit(p.kind, X, init=p) if p.kind != "plane" else fit("plane", X, init=p)
                if q is None or not np.isfinite(q.rms):
                    q = p
                q.faces, q.pts, q.area = faces, X, float(self.fa[faces].sum())
                prims[i] = q
            if not changed:
                break
        keep = [p for p in prims if len(p.faces) >= 2 or p.kind == "facets"]
        for p in prims:
            if not (len(p.faces) >= 2 or p.kind == "facets"):
                expelled.append(p.faces)
        if expelled:
            ex = np.concatenate(expelled)
            self.log(f"  validation expelled {len(ex)} faces from drifted surfaces")
            leftovers = np.unique(np.r_[leftovers, ex]).astype(int)
        return keep, leftovers

    def merge_extrusions(self, prims, leftovers):
        """Replace clusters of poorly fitting curved regions (and leftover faces) by spline extrusions."""
        LEFT = -2
        owner = np.full(len(self.m.faces), -1)
        for i, p in enumerate(prims):
            owner[p.faces] = i
        owner[leftovers] = LEFT
        poor = {i for i, p in enumerate(prims) if p.kind in ("cylinder", "cone", "sphere", "torus")
                and p.rms > 0.05 * self.tol}
        if not poor and not len(leftovers):
            return prims, leftovers
        axes = [prims[i].p["a"] for i in poor if prims[i].kind == "cylinder"]
        members = set(poor) | ({LEFT} if len(leftovers) else set())
        grown = set()
        while True:                      # grow through adjacent cylinders parallel to the poor ones
            faces = np.flatnonzero(np.isin(owner, list(members | grown)))
            nb = set(owner[self.adj[faces].indices]) - {-1, LEFT} - members - grown
            big = 4 * max([prims[i].area for i in poor] + [self.fa[leftovers].sum() if len(leftovers) else 0])
            add = {j for j in nb if (prims[j].kind == "cylinder" and
                                     any(abs(prims[j].p["a"] @ a) > 0.999 for a in axes)) or
                   (prims[j].kind == "plane" and prims[j].area < big and
                    any(abs(prims[j].p["n"] @ a) < 1e-3 for a in axes))}
            if not add:
                break
            grown |= add

        def clusters(groups):
            faces = np.flatnonzero(np.isin(owner, list(groups)))
            if not len(faces):
                return []
            n, lab = connected_components(self.adj[faces][:, faces], directed=False)
            return [faces[lab == c] for c in range(n)]

        new_prims, drop, used_left = [], set(), set()

        def try_fit(region):
            ids = set(owner[region]) - {LEFT}
            if len(ids) + int((owner[region] == LEFT).any()) < 2:
                return False
            X = self.m.vertices[np.unique(self.m.faces[region])]
            # a union that one cylinder explains stays a cylinder (fragmented round hole)
            cyl_ids = [i for i in ids if prims[i].kind == "cylinder"]
            if cyl_ids:
                init = max((prims[i] for i in cyl_ids), key=lambda p: p.area)
                q = fit("cylinder", X, init=init)
                if q is not None and np.abs(q.distance(X)).max() < 1.5 * self.tol:
                    q.faces, q.pts, q.area = region, X, float(self.fa[region].sum())
                    new_prims.append(q)
                    drop.update(ids)
                    used_left.update(region[owner[region] == LEFT].tolist())
                    return True
            best = None
            axes = [(p.p.get("c", p.p.get("v")), p.p["a"]) for p in prims
                    if p.kind in ("cylinder", "cone", "torus") and p.rms < 0.05 * self.tol]
            for fitter in (fit_extrusion, lambda m, f, t: fit_revolution(m, f, t, axes=axes)):
                try:
                    cand = fitter(self.m, region, self.tol)
                except Exception:
                    cand = None
                if cand is None:
                    continue
                err = np.abs(cand.distance(X)).max()
                if err < 1.5 * self.tol and (best is None or err < best[0]):
                    best = (err, cand)
            if best is None:
                return False
            q = best[1]
            q.faces, q.pts, q.area = region, X, float(self.fa[region].sum())
            new_prims.append(q)
            drop.update(ids)
            used_left.update(region[owner[region] == LEFT].tolist())
            return True

        def smooth_parts(region):
            n, lab = connected_components(self.adj_smooth[region][:, region], directed=False)
            return [region[lab == c] for c in range(n)]

        for region in clusters(members | grown):
            if try_fit(region):
                continue
            parts = [q for q in smooth_parts(region) if len(q) >= 4]
            if len(parts) > 1 and sum(try_fit(sub) for sub in parts) > 0:
                continue
            for sub in clusters(members):
                if np.isin(sub, region).all():
                    try_fit(sub)
        if new_prims:
            prims = [p for i, p in enumerate(prims) if i not in drop] + new_prims
            leftovers = np.array(sorted(set(leftovers.tolist()) - used_left), int)
            ne = sum(q.kind == "extrusion" for q in new_prims)
            nr = sum(q.kind == "revolution" for q in new_prims)
            self.log(f"  merged {len(drop)} curved fragments into {ne} spline extrusion(s), {nr} spline "
                     f"revolution(s) and {len(new_prims) - ne - nr} cylinder(s)")
        return prims, leftovers

    def absorb(self, prims, leftovers):
        """Merge small regions (and leftover faces) into larger surfaces that already explain them."""
        order = sorted(range(len(prims)), key=lambda i: -prims[i].area)
        alive = [True] * len(prims)
        n_abs = 0
        for i in reversed(order):                     # smallest first
            small = prims[i]
            for j in order:
                if j == i or not alive[j] or prims[j].area < 2 * small.area:
                    continue
                if len(self.inliers(prims[j], small.faces, scale=2.5)) == len(small.faces):
                    big = prims[j]
                    faces = np.r_[big.faces, small.faces]
                    X = self.m.vertices[np.unique(self.m.faces[faces])]
                    q = fit(big.kind, X, init=big) if big.kind != "plane" else fit("plane", X, init=big)
                    if q is None or np.abs(q.distance(X)).max() > 2.5 * self.tol:
                        continue
                    q.faces, q.pts, q.area = faces, X, float(self.fa[faces].sum())
                    prims[j] = q
                    alive[i] = False
                    n_abs += 1
                    break
        prims = [p for p, a in zip(prims, alive) if a]
        if len(leftovers):
            keep = np.ones(len(leftovers), bool)
            for j in sorted(range(len(prims)), key=lambda i: -prims[i].area):
                ok = self.inliers(prims[j], leftovers[keep], scale=2.5)
                if len(ok):
                    pr = prims[j]
                    faces = np.r_[pr.faces, ok]
                    X = self.m.vertices[np.unique(self.m.faces[faces])]
                    pr.faces, pr.pts, pr.area = faces, X, float(self.fa[faces].sum())
                    keep &= ~np.isin(leftovers, ok)
                    n_abs += len(ok)
            leftovers = leftovers[keep]
        self.log(f"  absorbed {n_abs} fragments into larger surfaces")
        return prims, leftovers

    def face_error(self, pr, faces):
        F = self.m.faces[faces]
        dv = np.abs(pr.distance(self.m.vertices[F.reshape(-1)])).reshape(-1, 3).max(1)
        dc = np.abs(pr.distance(self.fc[faces]))
        return np.maximum(dv, dc)

    def compete(self, prims, rounds=3):
        """Reassign boundary faces to the adjacent primitive that fits them best, then refit."""
        owner = np.full(len(self.m.faces), -1)
        for i, p in enumerate(prims):
            owner[p.faces] = i
        adj = self.m.face_adjacency
        moved_total = 0
        for _ in range(rounds):
            a, b = owner[adj[:, 0]], owner[adj[:, 1]]
            mask = (a != b) & (a >= 0) & (b >= 0)
            cand = {}
            for f, g in ((adj[mask, 0], b[mask]), (adj[mask, 1], a[mask])):
                for fi, gi in zip(f, g):
                    cand.setdefault(fi, set()).add(gi)
            moved = 0
            for f, others in cand.items():
                cur = owner[f]
                best, err = cur, self.face_error(prims[cur], np.array([f]))[0]
                for o in others:
                    e = self.face_error(prims[o], np.array([f]))[0]
                    if e < 0.5 * err and e < self.tol:
                        best, err = o, e
                if best != cur:
                    owner[f] = best
                    moved += 1
            moved_total += moved
            for i, p in enumerate(prims):
                faces = np.flatnonzero(owner == i)
                if len(faces) == 0:
                    continue
                X = self.m.vertices[np.unique(self.m.faces[faces])]
                q = fit(p.kind, X, init=p) if p.kind != "plane" else fit("plane", X, init=p)
                if q is not None and np.isfinite(q.rms):
                    q.faces, q.pts, q.area = faces, X, float(self.fa[faces].sum())
                    prims[i] = q
                else:
                    p.faces, p.pts, p.area = faces, X, float(self.fa[faces].sum())
            if moved == 0:
                break
        prims = [p for p in prims if len(p.faces)]
        self.log(f"  region competition moved {moved_total} boundary faces")
        return prims

    def explain_leftovers(self, faces):
        out = []
        if len(faces) == 0:
            return out
        sub = self.adj[faces][:, faces]
        n, lab = connected_components(sub, directed=False)
        for c in range(n):
            reg = faces[lab == c]
            X = self.m.vertices[np.unique(self.m.faces[reg])]
            done = False
            if "torus" in self.kinds and len(reg) >= 8:
                init = torus_init(self.fc[reg], self.fn[reg])
                if init is not None:
                    q = fit("torus", X, init=init)
                    if q is not None and np.isfinite(q.rms) and np.abs(q.distance(X)).max() < 2 * self.tol \
                            and q.p["r"] < self.diag and q.p["R"] < 2 * self.diag:
                        q.faces, q.pts, q.area = reg, X, float(self.fa[reg].sum())
                        out.append(q)
                        done = True
            if not done and self.allow_extrusion and len(reg) >= 8:
                q = fit_extrusion(self.m, reg, self.tol)
                if q is not None and np.abs(q.distance(X)).max() < 2 * self.tol:
                    q.faces, q.area = reg, float(self.fa[reg].sum())
                    out.append(q)
                    done = True
            if not done:
                q = Prim("facets")
                q.faces, q.pts, q.area, q.rms = reg, X, float(self.fa[reg].sum()), 0.0
                q.p["tris"] = self.m.triangles[reg]
                out.append(q)
        return out


# ---- design-intent snapping -------------------------------------------------
def reference_directions(axes, angle_step_deg=15.0):
    """World and part-frame axes, plus directions obtained by rotating one axis toward another
    in fixed angular steps."""
    dirs = [axes[:, i] for i in range(3)]
    for A in (np.eye(3),):
        dirs += [A[:, i] for i in range(3)]
        if angle_step_deg:
            for i in range(3):
                for j in range(3):
                    if i == j:
                        continue
                    for k in range(1, int(round(90 / angle_step_deg))):
                        t = np.radians(k * angle_step_deg)
                        dirs.append(np.cos(t) * A[:, i] + np.sin(t) * A[:, j])
    D = np.array(dirs)
    return D / np.linalg.norm(D, axis=1, keepdims=True)


def simplify_degenerate(prims, tol, log=print):
    """Replace degenerate fits by the simpler surface they approximate:
    tori with a vanishing major radius -> spheres, cones with a tiny half-angle -> cylinders,
    cones with a half-angle near 90 deg -> planes."""
    from .primitives import Prim
    n = 0
    for i, p in enumerate(prims):
        q = None
        if p.kind == "torus" and (p.p["R"] < 0.1 * p.p["r"] or p.p["R"] < 2 * tol):
            q = fit("sphere", p.pts, init=Prim("sphere", c=p.p["c"], r=p.p["r"]))
        elif p.kind == "cone" and p.p["t"] < np.radians(1.0):
            h = (p.pts - p.p["v"]) @ p.p["a"]
            r0 = float(np.median(h)) * np.tan(p.p["t"])
            q = fit("cylinder", p.pts, init=Prim("cylinder", c=p.p["v"] + np.median(h) * p.p["a"], a=p.p["a"], r=r0))
        elif p.kind == "cone" and p.p["t"] > np.radians(89.0):
            q = fit("plane", p.pts, init=Prim("plane", n=p.p["a"], d=float(p.pts.mean(0) @ p.p["a"])))
        if q is not None and np.isfinite(q.rms) and q.rms < 2 * tol:
            q.faces, q.pts, q.area = p.faces, p.pts, p.area
            prims[i] = q
            n += 1
    if n:
        log(f"  simplified {n} degenerate surface fits")
    return prims


def _fit_line2(Q):
    c = Q.mean(0)
    _, _, vt = np.linalg.svd(Q - c)
    return c, vt[0]


def _fit_circle2(Q):
    A = np.c_[2 * Q, np.ones(len(Q))]
    b = (Q ** 2).sum(1)
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    c = x[:2]
    return c, float(np.sqrt(max(x[2] + c @ c, 0.0)))


def _profile_pieces(Q, tol, max_pieces=5, seed=0):
    """Greedy RANSAC decomposition of 2D profile points (rho, z) into lines and circular arcs."""
    rng = np.random.default_rng(seed)
    rest = np.arange(len(Q))
    span = float(np.ptp(Q, axis=0).max()) + 1e-12
    pieces = []
    from scipy.spatial import cKDTree
    while len(rest) >= 3 and len(pieces) < max_pieces:
        R = Q[rest]
        best = None
        # minimal samples are drawn from a point and its nearest neighbours, so that every
        # candidate is a local fit (a circle through points of two different arcs is not proposed)
        _, nbr = cKDTree(R).query(R, k=min(len(R), 9))
        def local(n):
            i = rng.integers(len(R))
            return np.r_[i, rng.choice(nbr[i][1:], n - 1, replace=False)]
        for _ in range(300):
            i, j = local(2)
            d = R[j] - R[i]
            if np.linalg.norm(d) < tol:
                continue
            nrm = np.array([-d[1], d[0]]) / np.linalg.norm(d)
            inl = np.abs((R - R[i]) @ nrm) < tol
            if best is None or inl.sum() > best[1].sum():
                best = ("line", inl)
        for _ in range(300):
            ijk = local(3)
            try:
                c, r = _fit_circle2(R[ijk])
            except np.linalg.LinAlgError:
                continue
            if not np.isfinite(r) or r < tol or r > 2 * span:
                continue
            inl = np.abs(np.linalg.norm(R - c, axis=1) - r) < tol
            if best is None or inl.sum() > best[1].sum() + 1:      # prefer lines on ties
                best = ("arc", inl)
        if best is None or best[1].sum() < 3:
            break
        kind, inl = best
        for _ in range(3):                                         # refit and re-collect inliers
            if kind == "line":
                c, d = _fit_line2(R[inl])
                nrm = np.array([-d[1], d[0]])
                new = np.abs((R - c) @ nrm) < tol
            else:
                c, r = _fit_circle2(R[inl])
                new = np.abs(np.linalg.norm(R - c, axis=1) - r) < tol
            if new.sum() < 3:
                break
            inl = new
        if kind == "line":
            c, d = _fit_line2(R[inl])
            pieces.append(("line", c, d))
        else:
            c, r = _fit_circle2(R[inl])
            pieces.append(("arc", c, r))
        rest = rest[~inl]
    return pieces


def _piece_to_prim(piece, c3, a, tol):
    from .primitives import Prim
    kind = piece[0]
    if kind == "line":
        q, d = piece[1], piece[2]
        if abs(d[0]) < np.sin(np.radians(0.5)):                    # rho constant
            return Prim("cylinder", c=c3, a=a, r=float(q[0]))
        if abs(d[1]) < np.sin(np.radians(0.5)):                    # z constant
            return Prim("plane", n=a, d=float(a @ c3 + q[1]))
        if d[1] < 0:
            d = -d
        z_apex = q[1] - q[0] * d[1] / d[0]
        ax = a if d[0] > 0 else -a
        return Prim("cone", v=c3 + z_apex * a, a=ax, t=float(np.arctan(abs(d[0] / d[1]))))
    cc, r = piece[1], piece[2]
    if abs(cc[0]) < tol:
        return Prim("sphere", c=c3 + cc[1] * a, r=r)
    if r >= cc[0]:
        return None
    return Prim("torus", c=c3 + cc[1] * a, a=a, R=float(cc[0]), r=r)


def _piece_to_prim_extrusion(piece, o, u, w, a):
    from .primitives import Prim
    if piece[0] == "line":
        q, d = piece[1], piece[2]
        n2 = np.array([-d[1], d[0]])
        n = n2[0] * u + n2[1] * w
        n /= np.linalg.norm(n)
        return Prim("plane", n=n, d=float(n @ (o + q[0] * u + q[1] * w)))
    cc, r = piece[1], piece[2]
    return Prim("cylinder", c=o + cc[0] * u + cc[1] * w, a=a, r=r)


def facets_to_cylinders(prims, mesh, tol, max_step_deg=15.0, log=print):
    """A coarsely tessellated cylinder whose triangles span its full length is segmented as a ring of
    narrow planes. Adjacent plane regions whose normals turn by at most max_step_deg about one common
    axis, and whose vertices lie on one cylinder within tol, are replaced by that cylinder. The step
    limit keeps genuine polygonal prisms (hexagons, octagons) as planes."""
    import networkx as nx
    from .primitives import Prim
    own = np.full(len(mesh.faces), -1)
    for i, p in enumerate(prims):
        own[p.faces] = i
    adj = mesh.face_adjacency
    a, b = own[adj[:, 0]], own[adj[:, 1]]
    G = nx.Graph()
    step = np.cos(np.radians(max_step_deg))
    for x, y in set(zip(a.tolist(), b.tolist())):
        if x < 0 or y < 0 or x == y or prims[x].kind != "plane" or prims[y].kind != "plane":
            continue
        c = abs(prims[x].p["n"] @ prims[y].p["n"])
        if step <= c < np.cos(np.radians(0.5)):
            G.add_edge(x, y)
    drop, new = set(), []
    for comp in nx.connected_components(G):
        comp = sorted(comp)
        if len(comp) < 4:
            continue
        N = np.array([prims[i].p["n"] for i in comp])
        _, sv, vt = np.linalg.svd(N, full_matrices=False)
        ax = vt[-1]
        if np.abs(N @ ax).max() > np.sin(np.radians(1.0)):
            continue
        F = np.concatenate([prims[i].faces for i in comp])
        X = mesh.vertices[np.unique(mesh.faces[F])]
        Pn = X - np.outer(X @ ax, ax)
        c0 = Pn.mean(0)
        r0 = float(np.linalg.norm(Pn - c0, axis=1).mean())
        q = fit("cylinder", X, init=Prim("cylinder", c=c0, a=ax, r=r0))
        if q is None or not np.isfinite(q.rms) or np.abs(q.distance(X)).max() > tol:
            continue
        C = mesh.triangles_center[F]
        if np.abs(q.distance(C)).max() > 2 * tol:
            continue
        q.faces, q.pts, q.area = F, X, float(mesh.area_faces[F].sum())
        drop.update(comp)
        new.append(q)
    if not new:
        return prims
    log(f"  replaced {len(drop)} facet planes by {len(new)} cylinder(s)")
    return [p for i, p in enumerate(prims) if i not in drop] + new


def spheres_to_cylinders(prims, mesh, tol, log=print):
    """A short bore tessellated by triangles spanning its full length has all vertices on two circles,
    which a sphere fits as well as the true cylinder. If the region normals are all perpendicular to
    one axis and a cylinder fits vertices and centroids better, use the cylinder."""
    from .primitives import Prim
    n_rep = 0
    for i, p in enumerate(prims):
        if p.kind != "sphere" or len(p.faces) < 6:
            continue
        N = mesh.face_normals[p.faces]
        _, sv, vt = np.linalg.svd(N, full_matrices=False)
        if sv[-1] / np.sqrt(len(N)) > 0.05:          # normals not all perpendicular to one axis
            continue
        a = vt[-1]
        C = mesh.triangles_center[p.faces]
        X = np.vstack([p.pts, C])
        q = fit("cylinder", X, init=Prim("cylinder", c=p.p["c"], a=a, r=p.p["r"]))
        if q is None or not np.isfinite(q.rms):
            continue
        if np.abs(q.distance(X)).max() < tol and np.abs(q.distance(C)).max() < np.abs(p.distance(C)).max():
            q.faces, q.pts, q.area = p.faces, p.pts, p.area
            prims[i] = q
            n_rep += 1
    if n_rep:
        log(f"  replaced {n_rep} sphere fit(s) by cylinders")
    return prims


def _decompose_one(p, mesh, tol, tp):
    """Decompose one spline revolution / extrusion with 2D pieces fitted at tolerance tp.
    Returns (pieces as prims with faces, worst face error) or None."""
    V = mesh.vertices[np.unique(mesh.faces[p.faces])]
    if p.kind == "revolution":
        c3, a = p.p["c"], p.p["a"]
        W = V - c3
        z = W @ a
        Q2 = np.c_[np.linalg.norm(W - np.outer(z, a), axis=1), z]
    else:
        o, u, w, a = p.p["o"], p.p["u"], p.p["w"], p.p["a"]
        Q2 = np.c_[(V - o) @ u, (V - o) @ w]
    Q = np.unique(np.round(Q2 / (0.05 * tp)) * (0.05 * tp), axis=0)
    pieces = _profile_pieces(Q, tp, max_pieces=8)
    if p.kind == "revolution":
        new = [q for q in (_piece_to_prim(pc, c3, a, tol) for pc in pieces) if q is not None]
    else:
        new = [_piece_to_prim_extrusion(pc, o, u, w, a) for pc in pieces]
    if not (len(new) == len(pieces) and len(new) > 0):
        return None
    F = p.faces
    T = mesh.vertices[mesh.faces[F]]                                    # (n, 3, 3)

    def assign(new):
        D = np.stack([np.abs(q.distance(T.reshape(-1, 3))).reshape(-1, 3).max(1) for q in new], 1)
        own = D.argmin(1)
        err = D[np.arange(len(F)), own]
        return own, err
    own, err = assign(new)
    good = err < 1.5 * tol
    if not good.all():
        # faces whose vertices all lie on piece ends (a flat strip spanned by a few long
        # triangles) give too few profile points for RANSAC: fit one line through them
        Vb = mesh.vertices[np.unique(mesh.faces[F[~good]])]
        if p.kind == "revolution":
            Wb = Vb - c3
            zb = Wb @ a
            Qb = np.c_[np.linalg.norm(Wb - np.outer(zb, a), axis=1), zb]
        else:
            Qb = np.c_[(Vb - o) @ u, (Vb - o) @ w]
        cl, dl = _fit_line2(Qb)
        if np.abs((Qb - cl) @ np.array([-dl[1], dl[0]])).max() < tp:
            q = _piece_to_prim(("line", cl, dl), c3, a, tol) if p.kind == "revolution" \
                else _piece_to_prim_extrusion(("line", cl, dl), o, u, w, a)
            if q is not None:
                new = new + [q]
                own, err = assign(new)
                good = err < 1.5 * tol
    if not good.all():
        return None
    kept = []
    for k, q in enumerate(new):
        f = F[own == k]
        if len(f) == 0:
            continue
        q.faces = f
        q.pts = mesh.vertices[np.unique(mesh.faces[f])]
        q.area = float(mesh.area_faces[f].sum())
        q.rms = float(np.sqrt(np.mean(q.distance(q.pts) ** 2)))
        kept.append(q)
    return kept, float(err.max())


def _has_micro_arc(pieces, tol):
    """Arcs of radius below 5 tol are fragments of a noisy profile, not features."""
    for q in pieces:
        if q.kind in ("cylinder", "sphere") and q.p["r"] < 5 * tol:
            return True
        if q.kind == "torus" and q.p["r"] < 5 * tol:
            return True
    return False


def decompose_revolutions(prims, mesh, tol, log=print):
    """Replace a spline revolution or spline extrusion by analytic surfaces when its profile
    consists of lines and circular arcs (revolution: cylinder, plane, cone, torus, sphere;
    extrusion: planes and cylinders).

    The 2D pieces are fitted at tol, and again at tol/4. At tol one large arc can stand in for a
    line-fillet-line sequence (slot ends with small corner fillets); the tight result is used when
    it explains every face with less than half the error of the loose one."""
    out, n_rep, n_new = [], 0, 0
    for p in prims:
        if p.kind not in ("revolution", "extrusion"):
            out.append(p)
            continue
        loose = _decompose_one(p, mesh, tol, tol)
        tight = _decompose_one(p, mesh, tol, 0.25 * tol)
        best = loose
        if tight is not None and not _has_micro_arc(tight[0], tol) \
                and (loose is None or tight[1] < 0.5 * loose[1]):
            best = tight
        if best is None:
            out.append(p)
            continue
        out += best[0]
        n_rep += 1
        n_new += len(best[0])
    if n_rep:
        log(f"  decomposed {n_rep} spline surface(s) into {n_new} analytic surfaces")
    return out


def _world_aligned(a, deg=1.0):
    return np.max(np.abs(np.asarray(a, float))) > np.cos(np.radians(deg))


def _image_of(p, q, T, tol):
    """True if q is (nearly) the image of p under T."""
    if p.kind != q.kind:
        return False
    m = p.transformed(T)
    if p.kind == "plane":
        return abs(m.p["n"] @ q.p["n"]) > np.cos(np.radians(1)) and abs(m.p["d"] * np.sign(m.p["n"] @ q.p["n"]) - q.p["d"]) < 3 * tol
    if abs(abs(m.p["a"] @ q.p["a"]) - 1) > 1 - np.cos(np.radians(1)):
        return False
    key = "v" if p.kind == "cone" else "c"
    w = q.p[key] - m.p[key]
    a = q.p["a"]
    if np.linalg.norm(w - (w @ a) * a) > 3 * tol:
        return False
    if p.kind == "cylinder":
        return abs(p.p["r"] - q.p["r"]) < 3 * tol
    if p.kind == "cone":
        return abs(w @ a) < 3 * tol and abs(p.p["t"] - q.p["t"]) < 0.02
    if p.kind == "torus":
        return abs(w @ a) < 3 * tol and abs(p.p["R"] - q.p["R"]) < 3 * tol and abs(p.p["r"] - q.p["r"]) < 3 * tol
    return False


def symmetric_consensus(prims, tol, group, log=print):
    """Surfaces that are images of each other under the detected symmetry are fitted jointly:
    the regions are mapped into one frame, fitted once, and the result mapped back. Applied only
    when the joint fit stays within tol, so the symmetry is used, not imposed."""
    if not group:
        return prims
    kinds = ("plane", "cylinder", "cone", "torus")
    idx = [i for i, p in enumerate(prims) if p.kind in kinds]
    parent = {i: i for i in idx}
    link = {}

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for a_, i in enumerate(idx):
        for j in idx[a_ + 1:]:
            if find(i) == find(j):
                continue
            for T in group:
                if _image_of(prims[i], prims[j], T, tol):
                    parent[find(j)] = find(i)
                    link.setdefault(i, []).append((j, T))
                    break
    classes = {}
    for i in idx:
        classes.setdefault(find(i), []).append(i)
    n = 0
    for r, members in classes.items():
        if len(members) < 2:
            continue
        # transform from the representative to each member (breadth-first over the links)
        Tm = {r: np.eye(4)}
        frontier = [r]
        while frontier:
            i = frontier.pop()
            for j, T in link.get(i, []):
                if j not in Tm:
                    Tm[j] = T @ Tm[i]
                    frontier.append(j)
        if len(Tm) < len(members):
            continue
        # every copy already fitted cleanly: nothing to gain, and a joint fit would only perturb
        # surfaces that the later snapping steps make coaxial or coplanar exactly
        if max(prims[m].rms for m in members) < 0.25 * tol:
            continue
        P = []
        for m in members:
            Ti = np.linalg.inv(Tm[m])
            P.append(prims[m].pts @ Ti[:3, :3].T + Ti[:3, 3])
        X = np.vstack(P)
        # candidates: the joint fit and every member's own fit mapped into the representative's
        # frame; the one that explains all copies best wins (a clean copy beats a joint fit
        # polluted by seam faces of another copy)
        cands = []
        try:
            q = fit(prims[r].kind, X, init=prims[r])
            if q is not None and np.isfinite(q.rms):
                cands.append(q)
        except Exception:
            pass
        for m in members:
            cands.append(prims[m].transformed(np.linalg.inv(Tm[m])))
        scored = [(np.abs(c.distance(X)).max(), k) for k, c in enumerate(cands)]
        dbest, kbest = min(scored)
        if dbest > tol:
            continue
        q = cands[kbest]
        for m in members:
            mm = q.transformed(Tm[m])
            if "a" in mm.p and "a" in prims[m].p and mm.p["a"] @ prims[m].p["a"] < 0 and mm.kind != "cone":
                mm.p["a"] = -mm.p["a"]
            prims[m].p = mm.p
            prims[m].rms = float(np.sqrt(np.mean(prims[m].distance(prims[m].pts) ** 2)))
        n += 1
    if n:
        log(f"  symmetric consensus: {n} class(es) of mirror-image surfaces fitted jointly")
    return prims


def _arc_span_deg(p):
    from .primitives import perp_basis
    a = p.p["a"]
    u, w = perp_basis(a)
    W = p.pts - p.p["c"]
    ang = np.sort(np.mod(np.arctan2(W @ w, W @ u), 2 * np.pi))
    gaps = np.diff(np.r_[ang, ang[0] + 2 * np.pi])
    return float(np.degrees(2 * np.pi - gaps.max()))


def refit_odd_regions(prims, mesh, tol, log=print, seed=0, min_cover=0.9, group=None):
    """Re-segment odd regions with world-axis planes and cylinders.

    Coarse tessellations (long triangles across a short vertical band, a rounded corner of a few
    facets) let RANSAC explain a vertical cylinder or a side plane with spheres, tilted cones or
    facets. The faces of such regions are pooled and searched again for planes and cylinders whose
    normal or axis is a world axis. A region is replaced only when the new surfaces explain at
    least `min_cover` of its faces; otherwise it is kept, so genuinely tilted features stay."""
    from scipy.sparse.csgraph import connected_components
    from .primitives import Prim
    odd = []
    for i, p in enumerate(prims):
        if p.kind == "facets" or p.kind == "sphere":
            odd.append(i)
        elif p.kind in ("cone", "cylinder") and not _world_aligned(p.p["a"]):
            odd.append(i)
        elif p.kind == "plane" and len(p.faces) <= 3 and not _world_aligned(p.p["n"]):
            odd.append(i)
        elif p.kind == "cylinder" and len(p.pts) and _arc_span_deg(p) < 20 and p.p["r"] > 10 * tol:
            odd.append(i)                                    # a narrow arc band: plane-like, suspicious
        elif p.kind in ("cylinder", "cone", "torus") and len(p.pts):
            # mesh vertices of an exact export lie on the true surfaces: a fit with a sizeable
            # residual mixes faces of different surfaces
            if np.sqrt(np.mean(p.distance(p.pts) ** 2)) > 0.25 * tol:
                odd.append(i)
    # symmetric images of the good surfaces are hypotheses too (the group is a hint, not imposed)
    images = []
    if group:
        for i, p in enumerate(prims):
            if i in odd or p.kind not in ("plane", "cylinder", "cone", "torus"):
                continue
            for T in group:
                try:
                    images.append(p.transformed(T))
                except Exception:
                    pass
    if not odd:
        return prims
    pool = np.unique(np.concatenate([prims[i].faces for i in odd]))
    if len(pool) < 3:
        return prims
    rng = np.random.default_rng(seed)
    V, F = mesh.vertices, mesh.faces
    fn, fa = mesh.face_normals, mesh.area_faces
    diag = float(np.linalg.norm(mesh.extents))
    adj = mesh.face_adjacency
    in_pool = np.zeros(len(F), bool)
    in_pool[pool] = True
    cosang = np.cos(np.radians(25))
    E = np.eye(3)

    def inliers(pr, cand):
        d = np.abs(pr.distance(V[F[cand].reshape(-1)])).reshape(-1, 3).max(1)
        nr = pr.normal(mesh.triangles_center[cand])
        ok = (d < tol) & (np.abs(np.einsum("ij,ij->i", nr, fn[cand])) > cosang)
        sel = cand[ok]
        if len(sel) < 2:
            return sel
        loc = {f: k for k, f in enumerate(sel)}
        m = np.isin(adj[:, 0], sel) & np.isin(adj[:, 1], sel)
        from scipy.sparse import coo_matrix
        a_ = np.array([loc[x] for x in adj[m, 0]], int)
        b_ = np.array([loc[x] for x in adj[m, 1]], int)
        g = coo_matrix((np.ones(len(a_)), (a_, b_)), shape=(len(sel), len(sel)))
        n, lab = connected_components(g, directed=False)
        areas = np.bincount(lab, weights=fa[sel])
        return sel[lab == np.argmax(areas)]

    remaining = pool.copy()
    found = []
    good = [i for i, p in enumerate(prims) if i not in odd and p.kind in ("plane", "cylinder", "cone", "torus")]
    absorbed = {}
    while len(remaining) >= 2:
        best = None
        for pr in images:
            sel = inliers(pr, remaining)
            if len(sel) >= 2:
                area = fa[sel].sum()
                if best is None or area > best[0]:
                    best = (area, pr, sel)
        for gi in good:                                       # existing surfaces absorb pool faces
            sel = inliers(prims[gi], remaining)
            if len(sel) >= 1:
                area = fa[sel].sum()
                if best is None or area > best[0]:
                    best = (area, gi, sel)
        if best is not None and isinstance(best[1], (int, np.integer)):
            gi, sel = best[1], best[2]
            absorbed[gi] = np.r_[absorbed.get(gi, np.array([], int)), sel]
            remaining = np.setdiff1d(remaining, sel)
            continue
        for _ in range(200):
            i = remaining[rng.integers(len(remaining))]
            k = rng.integers(3)
            ax = E[k]
            if abs(fn[i] @ ax) > cosang:                      # plane with a world normal
                pr = Prim("plane", n=ax * np.sign(fn[i] @ ax), d=float(V[F[i]].mean(0) @ ax * np.sign(fn[i] @ ax)))
            elif abs(fn[i] @ ax) < np.sin(np.radians(25)):    # cylinder about a world axis
                near = remaining[np.linalg.norm(mesh.triangles_center[remaining] - mesh.triangles_center[i], axis=1) < 0.08 * diag]
                if len(near) < 2:
                    continue
                j = near[rng.integers(len(near))]
                P = np.vstack([V[F[i]], V[F[j]]])
                Q = P - np.outer(P @ ax, ax)
                u = E[(k + 1) % 3]
                w = E[(k + 2) % 3]
                q2 = np.c_[Q @ u, Q @ w]
                A = np.c_[2 * q2, np.ones(len(q2))]
                sol, *_ = np.linalg.lstsq(A, (q2 ** 2).sum(1), rcond=None)
                c2 = sol[:2]
                r = float(np.sqrt(max(sol[2] + c2 @ c2, 0)))
                if not (tol < r < diag):
                    continue
                pr = Prim("cylinder", c=c2[0] * u + c2[1] * w, a=ax, r=r)
            else:
                continue
            sel = inliers(pr, remaining)
            if len(sel) < 2:
                continue
            if pr.kind == "cylinder" and pr.p["r"] > 10 * tol:
                pr.pts = V[np.unique(F[sel])]
                if _arc_span_deg(pr) < 20:                    # a narrow band of a large cylinder is a plane
                    continue
            area = fa[sel].sum()
            if best is None or area > best[0]:
                best = (area, pr, sel)
        if best is None:
            break
        _, pr, sel = best
        import copy
        pr = copy.copy(pr)
        pr.p = dict(pr.p)
        X = V[np.unique(F[sel])]
        q = fit(pr.kind, X, init=pr, fixed_axis=True) if pr.kind in ("plane", "cylinder") else None
        if q is not None and np.isfinite(q.rms) and np.abs(q.distance(X)).max() < tol:
            q2_sel = inliers(q, remaining)
            if len(q2_sel) >= len(sel):
                pr, sel = q, q2_sel
        pr.faces = sel
        pr.pts = V[np.unique(F[sel])]
        pr.area = float(fa[sel].sum())
        pr.rms = float(np.sqrt(np.mean(pr.distance(pr.pts) ** 2)))
        found.append(pr)
        remaining = np.setdiff1d(remaining, sel)
        if fa[remaining].sum() < 1e-4 * mesh.area:
            break
    if not found and not absorbed:
        return prims
    def settle(cands):
        """Ownership with the regions not in cands kept; a found surface survives only with
        enough faces (two triangles fit any cylinder of large radius)."""
        own = np.full(len(F), -1)
        for k, q in enumerate(found):
            own[q.faces] = k
        for gi, sel in absorbed.items():
            own[sel] = len(found) + gi                        # marks: absorbed by an existing surface
        kept = [prims[i].faces for i in odd if i not in cands]
        if kept:
            own[np.concatenate(kept)] = -1                    # faces of kept regions stay with them
        for k, q in enumerate(found):
            if (own == k).sum() < (3 if q.kind == "cylinder" else 2):
                own[own == k] = -1
        return own
    # a region is replaced only if the surfaces that finally survive cover it; dropping one region
    # can starve a surface that covered another, so iterate to a fixed point
    replaced = list(odd)
    while True:
        owner = settle(replaced)
        keep = [i for i in replaced if (owner[prims[i].faces] >= 0).mean() < min_cover]
        if not keep:
            break
        replaced = [i for i in replaced if i not in keep]
    owner = settle(replaced)
    new = []
    for k, q in enumerate(found):
        f = np.flatnonzero(owner == k)
        if not len(f):
            continue
        q.faces = f
        q.pts = V[np.unique(F[f])]
        q.area = float(fa[f].sum())
        new.append(q)
    n_abs = 0
    for gi in absorbed:
        f = np.flatnonzero(owner == len(found) + gi)
        if len(f):
            g = prims[gi]
            g.faces = np.unique(np.r_[g.faces, f])
            g.pts = V[np.unique(F[g.faces])]
            g.area = float(fa[g.faces].sum())
            n_abs += len(f)
    if not replaced:
        return prims
    # faces of replaced regions that no new surface took become one facet region
    rest = np.concatenate([prims[i].faces for i in replaced])
    taken = np.concatenate([q.faces for q in new] + [prims[gi].faces for gi in absorbed])
    rest = np.setdiff1d(rest, taken)
    out = [p for i, p in enumerate(prims) if i not in replaced] + new
    if len(rest):
        fp = Prim("facets", tris=mesh.triangles[rest])
        fp.faces, fp.pts, fp.area, fp.rms = rest, V[np.unique(F[rest])], float(fa[rest].sum()), 0.0
        out.append(fp)
    log(f"  refit {len(replaced)} odd region(s) as {len(new)} plane(s)/cylinder(s)/tori"
        + (f", {n_abs} triangles absorbed by existing surfaces" if n_abs else "")
        + (f"; {len(rest)} triangles left as facets" if len(rest) else ""))
    return out


def cones_to_tori(prims, mesh, tol, log=print):
    """Replace coaxial cones (and spline revolutions) that together lie on one torus by that torus.

    A fillet of small tube radius is coarsely tessellated, and RANSAC often explains it with two
    cones of different half-angles; the least-squares torus through their union is exact."""
    from .primitives import Prim, torus_init
    cand = [i for i, p in enumerate(prims) if p.kind in ("cone", "revolution")]
    groups, seen = [], set()
    for i in cand:
        if i in seen or prims[i].kind != "cone":
            continue
        a, v = prims[i].p["a"], prims[i].p["v"]
        grp = [j for j in cand if j not in seen and prims[j].kind == "cone"
               and abs(abs(prims[j].p["a"] @ a) - 1) < 1e-3
               and np.linalg.norm(np.cross(prims[j].p["v"] - v, a)) < 2 * tol]
        seen.update(grp)
        groups.append(grp)
    groups += [[i] for i in cand if prims[i].kind == "revolution"]
    fc, fn = mesh.triangles_center, mesh.face_normals
    replace, drop = {}, set()
    for grp in groups:
        F = np.unique(np.concatenate([prims[j].faces for j in grp]))
        if len(F) < 8:
            continue
        init = torus_init(fc[F], fn[F])
        if init is None:
            continue
        X = mesh.vertices[np.unique(mesh.faces[F])]
        q = fit("torus", X, init=init)
        if q is None or not np.isfinite(q.rms) or q.p["R"] < 2 * tol or q.p["r"] < tol:
            continue
        old = np.sqrt(np.mean(np.concatenate([prims[j].distance(prims[j].pts) ** 2 for j in grp])))
        qmax = np.abs(q.distance(X)).max()
        # one exact torus for several cones is accepted when it fits all of them closely, even if
        # each cone fits its own narrow ring exactly
        if qmax < tol and (q.rms < 0.5 * old or (len(grp) >= 2 and qmax < 0.5 * tol)):
            q.faces, q.pts, q.area = F, X, float(sum(prims[j].area for j in grp))
            replace[grp[0]] = q
            drop.update(grp[1:])
    if not replace:
        return prims
    out = [replace.get(i, p) for i, p in enumerate(prims) if i not in drop]
    out, _ = merge_identical(out, tol)
    out = _grow_tori(out, mesh, tol)
    log(f"  replaced {len(replace) + len(drop)} cone/revolution fits by tori")
    return out


def _grow_tori(prims, mesh, tol, rounds=4):
    """Move faces next to a torus into it when they lie on it and their normal agrees with the
    torus better than with their own surface (the first row of a small fillet is within tolerance of
    the adjacent plane and is often absorbed by it); then refit the surfaces that lost faces."""
    own = np.full(len(mesh.faces), -1)
    for i, p in enumerate(prims):
        own[p.faces] = i
    adj = mesh.face_adjacency
    fc, fn = mesh.triangles_center, mesh.face_normals
    tori = [i for i, p in enumerate(prims) if p.kind == "torus"]
    lost = set()
    for _ in range(rounds):
        moved = 0
        for t in tori:
            T = prims[t]
            e = adj[(own[adj[:, 0]] == t) ^ (own[adj[:, 1]] == t)]
            cand = np.unique(np.where(own[e[:, 0]] == t, e[:, 1], e[:, 0]))
            cand = cand[own[cand] >= 0]
            for f in cand:
                o = prims[own[f]]
                if o.kind in ("facets", "extrusion", "revolution") or len(o.faces) <= 3:
                    continue
                if abs(T.distance(fc[f:f + 1])[0]) > tol:
                    continue
                ct = abs(T.normal(fc[f:f + 1])[0] @ fn[f])
                co = abs(o.normal(fc[f:f + 1])[0] @ fn[f])
                if ct > co + 1e-6:
                    lost.add(int(own[f]))
                    o.faces = o.faces[o.faces != f]
                    T.faces = np.r_[T.faces, f]
                    own[f] = t
                    moved += 1
        if not moved:
            break
    for i in lost | set(tori):
        p = prims[i]
        p.pts = mesh.vertices[np.unique(mesh.faces[p.faces])]
        p.area = float(mesh.area_faces[p.faces].sum())
        if i in lost and p.kind in ("plane", "cylinder", "cone", "sphere", "torus"):
            q = fit(p.kind, p.pts, init=p)
            if q is not None and np.isfinite(q.rms):
                p.p, p.rms = q.p, q.rms
    return prims


def snap_directions(prims, axes, snap_deg=1.0, angle_step_deg=15.0, log=print, tol=None):
    D = reference_directions(axes, angle_step_deg)
    n_snapped = 0
    for i, pr in enumerate(prims):
        key = "n" if pr.kind == "plane" else ("a" if pr.kind in ("cylinder", "cone", "torus") else None)
        if key is None:
            continue
        v = pr.p[key]
        dots = D @ v
        j = int(np.argmax(np.abs(dots)))
        if np.degrees(np.arccos(min(1.0, abs(dots[j])))) < snap_deg:
            pr.p[key] = D[j] * np.sign(dots[j])
            q = fit(pr.kind, pr.pts, init=pr, fixed_axis=True) if pr.kind != "plane" else \
                fit("plane", pr.pts, init=pr, fixed_axis=True)
            if q is not None and tol is not None and np.abs(q.distance(pr.pts)).max() > max(tol, 1.5 * np.abs(
                    pr.distance(pr.pts)).max()):
                pr.p[key] = v                       # the snapped direction does not fit: keep the fitted one
                continue
            if q is not None:
                q.faces, q.pts, q.area = pr.faces, pr.pts, pr.area
                prims[i] = q
                n_snapped += 1
    log(f"  snapped {n_snapped} surface directions to the part frame")
    return prims


def unify_parallel(prims, tol, log=print, angle_deg=1.0):
    """Merge near-duplicate surfaces into one shared surface: planes within angle_deg whose regions
    lie on each other's plane (within 2 tol) get one joint plane fit; cylinders/cones/tori with
    nearly parallel, nearly coincident axes share one axis (and cylinders of equal radius one fit)."""
    cos_a = np.cos(np.radians(angle_deg))
    n_merge = 0
    planes = [p for p in prims if p.kind == "plane"]
    groups, seen = [], set()
    for i, p in enumerate(planes):
        if i in seen:
            continue
        grp = [p]
        seen.add(i)
        for j in range(i + 1, len(planes)):
            q = planes[j]
            if j in seen or abs(p.p["n"] @ q.p["n"]) < cos_a:
                continue
            if np.abs(p.distance(q.pts)).max() < 2 * tol and np.abs(q.distance(p.pts)).max() < 2 * tol:
                grp.append(q)
                seen.add(j)
        if len(grp) > 1:
            groups.append(grp)
    for grp in groups:
        X = np.vstack([g.pts for g in grp])
        f = fit("plane", X, init=grp[0])
        if f is None or np.abs(f.distance(X)).max() > 2.5 * tol:
            continue
        for g in grp:
            sg = np.sign(g.p["n"] @ f.p["n"]) or 1.0
            g.p["n"], g.p["d"] = f.p["n"] * sg, f.p["d"] * sg
        n_merge += len(grp) - 1
    revs = [p for p in prims if p.kind in ("cylinder", "cone", "torus")]
    for i, p in enumerate(revs):
        cp = p.p.get("c", p.p.get("v"))
        for q in revs[i + 1:]:
            cq = q.p.get("c", q.p.get("v"))
            if abs(p.p["a"] @ q.p["a"]) < cos_a:
                continue
            w = cq - cp
            off = np.linalg.norm(w - (w @ p.p["a"]) * p.p["a"])
            if off >= 2 * tol:
                continue
            sg = np.sign(p.p["a"] @ q.p["a"]) or 1.0
            if p.kind == q.kind == "cylinder" and abs(p.p["r"] - q.p["r"]) < 2 * tol:
                X = np.vstack([p.pts, q.pts])
                f = fit("cylinder", X, init=p)
                if f is not None and np.abs(f.distance(X)).max() < 2.5 * tol:
                    for g in (p, q):
                        g.p["c"], g.p["a"], g.p["r"] = f.p["c"].copy(), f.p["a"] * (1 if g is p else sg), f.p["r"]
                    n_merge += 1
                    continue
            q.p["a"] = p.p["a"] * sg
            key = "c" if "c" in q.p else "v"
            q.p[key] = cp + ((cq - cp) @ p.p["a"]) * p.p["a"]
            n_merge += 1
    prims, n_dup = merge_identical(prims, tol)
    log(f"  unified {n_merge} coplanar/coaxial surfaces" + (f"; merged {n_dup} duplicate curved regions" if n_dup else ""))
    return prims


def _same_surface(p, q, tol):
    """Candidate test for two regions on one curved surface. Short regions are poorly conditioned,
    so the parameter test is loose; merge_identical accepts a pair only if a joint fit is tight."""
    if p.kind != q.kind or p.kind not in ("cylinder", "cone", "sphere", "torus"):
        return False
    r = float(p.p["r"]) if "r" in p.p else 1.0
    lt = max(5 * tol, 0.03 * r)
    if p.kind == "sphere":
        return np.linalg.norm(p.p["c"] - q.p["c"]) < lt and abs(p.p["r"] - q.p["r"]) < lt
    if abs(abs(p.p["a"] @ q.p["a"]) - 1) > 1 - np.cos(np.radians(1.0)):
        return False
    a = p.p["a"]
    key = "v" if p.kind == "cone" else "c"
    w = q.p[key] - p.p[key]
    if np.linalg.norm(w - (w @ a) * a) > lt:
        return False
    if p.kind == "cylinder":
        return abs(p.p["r"] - q.p["r"]) < lt
    if p.kind == "cone":
        return abs(w @ a) < lt and abs(p.p["t"] - q.p["t"]) < 0.02 and (p.p["a"] @ q.p["a"]) > 0
    return abs(w @ a) < lt and abs(p.p["R"] - q.p["R"]) < lt and abs(p.p["r"] - q.p["r"]) < lt


def merge_identical(prims, tol):
    """Regions lying on one curved surface (split apart by segmentation) become one region with
    one joint fit; separate copies would give coincident tool faces. Candidates are added one at a
    time and kept only while the joint fit stays within 2 tol; the (snapped) axis is kept fixed."""
    out, n = [], 0
    used = set()
    for i, p in enumerate(prims):
        if i in used:
            continue
        for j in range(i + 1, len(prims)):
            if j in used or not _same_surface(p, prims[j], tol):
                continue
            X = np.vstack([p.pts, prims[j].pts])
            init = p
            if "a" in p.p and prims[j].p["a"] @ p.p["a"] < 0 and p.kind == "cone":
                continue
            try:
                f = fit(p.kind, X, init=init, fixed_axis=p.kind in ("cylinder", "cone", "torus"))
            except TypeError:
                f = fit(p.kind, X, init=init)
            if f is None or not np.isfinite(f.rms) or np.abs(f.distance(X)).max() >= 2 * tol:
                continue
            if "a" in f.p and f.p["a"] @ p.p["a"] < 0:
                f.p["a"] = -f.p["a"]
            f.faces = np.unique(np.concatenate([p.faces, prims[j].faces]))
            f.pts, f.area = X, p.area + prims[j].area
            if hasattr(p, "is_hole"):
                f.is_hole = p.is_hole
            used.add(j)
            n += 1
            p = f
        out.append(p)
    return out, n


def snap_tangency(prims, mesh, tol, log=print):
    """Make cylinders exactly tangent to adjacent planes that are parallel to the cylinder axis."""
    face_owner = np.full(len(mesh.faces), -1)
    for i, p in enumerate(prims):
        face_owner[p.faces] = i
    adj = mesh.face_adjacency
    pairs = set()
    for f0, f1 in adj:
        a, b = face_owner[f0], face_owner[f1]
        if a >= 0 and b >= 0 and a != b:
            pairs.add((min(a, b), max(a, b)))
    # bridge through small regions: two regions that both touch a region of 6 faces or fewer
    # count as adjacent (coarse meshes leave thin junk regions along tangent seams)
    nbr = {}
    for (x, y) in list(pairs):
        nbr.setdefault(x, set()).add(y)
        nbr.setdefault(y, set()).add(x)
    for k, ns in nbr.items():
        if len(prims[k].faces) <= 6:
            ns = sorted(ns)
            for ii in range(len(ns)):
                for jj in range(ii + 1, len(ns)):
                    pairs.add((ns[ii], ns[jj]))
    count = 0
    tangent_pairs = []
    for i, p in enumerate(prims):
        if p.kind != "cylinder":
            continue
        a = p.p["a"]
        cons, cons_idx = [], []
        for (x, y) in pairs:
            if i not in (x, y):
                continue
            q = prims[y if x == i else x]
            if q.kind != "plane" or abs(q.p["n"] @ a) > np.sin(np.radians(0.5)):
                continue
            s = q.p["n"] @ p.p["c"] - q.p["d"]
            if abs(abs(s) - p.p["r"]) < 3 * tol:
                cons.append((q.p["n"], q.p["d"], np.sign(s)))
                cons_idx.append(y if x == i else x)
        if not cons:
            continue
        # make the axis exactly perpendicular to the tangent planes' normals
        a2 = a.copy()
        for n_, _, _ in cons:
            a2 = a2 - (a2 @ n_) * n_
        if np.linalg.norm(a2) < 0.5:
            continue
        a = a2 / np.linalg.norm(a2)
        p.p["a"] = a
        u, v = perp_basis(a)
        c0 = p.p["c"]
        # unknowns: (du, dv, r); minimise ||(du, dv, r - r0)||^2 s.t. sgn*(n.(c0 + du u + dv v) - d) = r
        A = np.array([[sg * (n @ u), sg * (n @ v), -1.0] for n, d, sg in cons])
        b = np.array([-sg * (n @ c0 - d) for n, d, sg in cons])
        x0 = np.array([0, 0, p.p["r"]])
        if len(cons) >= 3:
            x, *_ = np.linalg.lstsq(A, b, rcond=None)
        else:
            # minimum-norm correction from x0 satisfying A x = b (least squares when the two
            # constraints are dependent, e.g. two parallel planes: A A^T is then singular)
            dx, *_ = np.linalg.lstsq(A, b - A @ x0, rcond=None)
            x = x0 + dx
        if abs(x[2] - p.p["r"]) < 3 * tol and np.hypot(x[0], x[1]) < 3 * tol:
            p.p["c"] = c0 + x[0] * u + x[1] * v
            p.p["r"] = float(x[2])
            count += 1
            tangent_pairs += [(i, j) for j in cons_idx]
    log(f"  enforced tangency on {count} cylinders")
    return prims, tangent_pairs


def drop_seam_bridges(prims, mesh, pairs, log=print, max_faces=6):
    """Remove small regions (max_faces faces or fewer) that only bridge a tangent plane-cylinder
    pair. On coarse meshes the chordal seam triangles between a wall and its corner round deviate
    from both surfaces and are fitted as small tilted cylinders or tori; as tools they notch the
    solid. Returns (prims, pairs) with pairs reindexed."""
    owner = np.full(len(mesh.faces), -1)
    for i, p in enumerate(prims):
        owner[p.faces] = i
    nbr = {}
    for f0, f1 in mesh.face_adjacency:
        a, b = owner[f0], owner[f1]
        if a >= 0 and b >= 0 and a != b:
            nbr.setdefault(a, set()).add(b)
            nbr.setdefault(b, set()).add(a)
    tp = {tuple(sorted(x)) for x in pairs}
    in_pair = {i for x in tp for i in x}
    drop = set()
    for k, ns in nbr.items():
        if len(prims[k].faces) > max_faces or k in in_pair:
            continue
        if any(tuple(sorted((x, y))) in tp for x in ns for y in ns if x < y):
            drop.add(k)
    if not drop:
        return prims, pairs
    keep = [i for i in range(len(prims)) if i not in drop]
    new = {o: n for n, o in enumerate(keep)}
    log(f"  dropped {len(drop)} seam region(s) between tangent surfaces")
    return [prims[i] for i in keep], [(new[a], new[b]) for a, b in pairs]


def snap_grid(prims, origin, axes, grid, log=print):
    """Round radii and frame-space offsets of planes/axes to a grid (design dimensions)."""
    if not grid:
        return prims
    rnd = lambda x: float(np.round(x / grid) * grid)
    for p in prims:
        if p.kind == "plane":
            p.p["d"] = rnd(p.p["d"] - p.p["n"] @ origin) + p.p["n"] @ origin
        elif p.kind in ("cylinder", "torus", "sphere"):
            c = axes.T @ (p.p["c"] - origin)
            p.p["c"] = origin + axes @ np.array([rnd(x) for x in c])
            for k in ("r", "R"):
                if k in p.p:
                    p.p[k] = rnd(p.p[k])
    log(f"  snapped dimensions to a {grid} mm grid")
    return prims


# ---- symmetry ---------------------------------------------------------------
def symmetrize(prims, G, tol, diag, log=print, complete=True):
    """Refit each primitive on its whole orbit; with complete=True replace the set by complete orbits,
    otherwise keep only the images that matched an existing surface."""
    Ginv = [np.linalg.inv(T) for T in G]
    used = np.zeros(len(prims), bool)
    order = np.argsort([-p.area for p in prims])
    out = []
    matched_total = 0
    for i in order:
        if used[i]:
            continue
        P = prims[i]
        used[i] = True
        if P.kind in ("facets", "extrusion", "revolution"):
            if complete:
                for T in G:
                    out.append(P.transformed(T))
            else:
                out.append(P)
            continue
        pts = [P.pts]
        matched_T = [G[0]]
        for T, Ti in zip(G[1:], Ginv[1:]):
            img = P.transformed(T)
            best, bj = np.inf, -1
            for j in np.flatnonzero(~used):
                Q = prims[j]
                if Q.kind != P.kind:
                    continue
                s = match_score(img, Q, tol)
                if s < 2 * tol:
                    dd = cKDTree(Q.pts).query(img.pts)[0].mean()
                    if dd < 0.1 * diag and dd < best:
                        best, bj = dd, j
            if bj >= 0:
                used[bj] = True
                matched_total += 1
                Q = prims[bj]
                pts.append(Q.pts @ Ti[:3, :3].T + Ti[:3, 3])
                matched_T.append(T)
        X = np.vstack(pts)
        R = fit(P.kind, X, init=P) if P.kind != "plane" else fit("plane", X, init=P)
        if R is None:
            R = P
        R.pts, R.area, R.faces = X, P.area, P.faces
        orbit = []
        for T in (G if complete else matched_T):
            img = R.transformed(T)
            if not any(match_score(img, o, tol) < 0.5 * tol and
                       cKDTree(o.pts).query(img.pts)[0].mean() < 0.01 * diag for o in orbit):
                orbit.append(img)
        out += orbit
    log(f"  symmetric refit: {len(prims)} regions -> {len(out)} surfaces "
        f"({matched_total} matched images)")
    return out


def reassign_faces(prims, mesh, tol):
    """Give every mesh face to the closest surface (used after symmetrization, when images
    do not yet own mesh faces). Faces farther than 3 tol from all surfaces stay unowned."""
    fc = mesh.triangles_center
    best = np.full(len(fc), np.inf)
    owner = np.full(len(fc), -1)
    for i, p in enumerate(prims):
        if p.kind == "facets":
            continue
        d = np.abs(p.distance(fc))
        lo_, hi_ = p.pts.min(0) - 3 * tol, p.pts.max(0) + 3 * tol      # stay near the region
        near = np.all((fc >= lo_) & (fc <= hi_), axis=1)
        d[~near] = np.inf
        upd = d < best
        best[upd], owner[upd] = d[upd], i
    owner[best > 3 * tol] = -1
    for i, p in enumerate(prims):
        if p.kind == "facets":
            continue
        f = np.flatnonzero(owner == i)
        p.faces = f
        if len(f):
            p.area = float(mesh.area_faces[f].sum())
    return prims


def summarize(prims, mesh):
    from collections import Counter
    c = Counter(p.kind for p in prims)
    area = {k: sum(p.area for p in prims if p.kind == k) for k in c}
    tot = mesh.area
    s = ", ".join(f"{k} x{c[k]} ({100 * area[k] / tot:.1f}% area)" for k in sorted(c, key=lambda k: COMPLEXITY[k]))
    analytic = sum(area.get(k, 0) for k in ALGEBRAIC) / tot
    return s, analytic
