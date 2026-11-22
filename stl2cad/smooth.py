"""Smooth reconstruction for organic parts (topology optimization, lattices).

Pipeline:
  1. voxel occupancy -> Gaussian blur -> average over the symmetry group,
  2. skeleton of the symmetric field -> member graph (junctions, branches),
  3. bolt holes found by cylinder RANSAC; each becomes a round pad plus an exact hole,
  4. members are mapped into a fundamental domain of the symmetry group, copies are
     clustered and averaged, and each member is sized to the original's
     cross-sectional area (median over its symmetric copies),
  5. members become lofted tubes along smoothing splines, or (algebraic-only)
     chains of cones joined by spheres; junctions become spheres,
  6. the domain solid is patterned by the group, trimmed to the design domain
     (bounding slab), and holes are cut.
"""

import time

import networkx as nx
import numpy as np
import trimesh
from cadquery import Location, Plane, Shape, Solid, Vector
from cadquery import func as cf
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.gp import gp_Trsf
from scipy import ndimage
from scipy.interpolate import splev, splprep
from shapely.geometry import Point
from skimage.morphology import skeletonize

from .core import Result
from .primitives import perp_basis, unit
from .segment import Segmenter
from .symmetry import group_matrices


# ---- geometry helpers ---------------------------------------------------------
def occ_transform(shape, T):
    """Rigid or reflective transform; reflections are re-oriented so solids stay outward-facing."""
    tr = gp_Trsf()
    tr.SetValues(*[float(x) for x in T[:3, :4].reshape(-1)])
    out = Shape.cast(BRepBuilderAPI_Transform(shape.wrapped, tr, True).Shape())
    if np.linalg.det(T[:3, :3]) < 0 and out.Solids() and out.Volume() < 0:
        out = Shape.cast(out.wrapped.Reversed())
    return out


def tube(p, r):
    """Loft of circles normal to polyline p with radii r (parallel-transported frames)."""
    tang = np.gradient(p, axis=0)
    tang /= np.linalg.norm(tang, axis=1, keepdims=True)
    x = np.cross(tang[0], [0, 0, 1] if abs(tang[0][2]) < 0.9 else [1, 0, 0])
    secs = []
    for pi, ti, ri in zip(p, tang, r):
        x = x - np.dot(x, ti) * ti
        x /= np.linalg.norm(x)
        secs.append(cf.wire(cf.circle(float(ri))).moved(Location(Plane(Vector(*pi), Vector(*x), Vector(*ti)))))
    return cf.loft(secs, cap=True)


def cone_chain(p, r):
    """Algebraic-only member: truncated cones between stations, spheres at interior stations."""
    parts = []
    for i in range(len(p) - 1):
        d = p[i + 1] - p[i]
        L = float(np.linalg.norm(d))
        if L < 1e-6:
            continue
        parts.append(Solid.makeCone(float(r[i]), float(r[i + 1]), L, Vector(*p[i]), Vector(*(d / L))))
    for i in range(1, len(p) - 1):
        parts.append(cf.sphere(2 * float(r[i])).moved(Location(Vector(*p[i]))))
    out = parts[0]
    for q in parts[1:]:
        out = out + q
    return out


def fuse_all(pieces, log, label, fuzzy=None):
    """Sequential union; failed pieces are retried at the end with a fuzzy tolerance."""
    def ok(trial, base):
        return trial is not None and trial.isValid() and trial.Solids() and trial.Volume() >= base.Volume() - 1e-6
    out = pieces[0]
    failed = []
    for q in pieces[1:]:
        try:
            trial = out + q
        except Exception:
            trial = None
        if ok(trial, out):
            out = trial
        else:
            failed.append(q)
    skipped = 0
    for q in failed:
        trial = None
        for tol in ([fuzzy] if fuzzy else []) + [1e-5, 1e-4, 1e-3]:
            try:
                trial = cf.fuse(out, q, tol=tol)
            except Exception:
                trial = None
            if ok(trial, out):
                break
        if ok(trial, out):
            out = trial
        else:
            skipped += 1
    if skipped:
        log(f"  warning: {skipped} {label} skipped (boolean failure)")
    return out


# ---- stages ---------------------------------------------------------------------
def symmetric_field(mesh, G, pitch, blur, log):
    vg = mesh.voxelized(pitch).fill()
    F = ndimage.gaussian_filter(vg.matrix.astype(np.float32), sigma=blur / pitch)
    Tinv = np.linalg.inv(vg.transform)
    corners = np.array([[x, y, z] for x in mesh.bounds[:, 0] for y in mesh.bounds[:, 1] for z in mesh.bounds[:, 2]])
    allc = np.vstack([corners @ T[:3, :3].T + T[:3, 3] for T in G])
    lo, hi = allc.min(0) - 3 * blur, allc.max(0) + 3 * blur
    axes = [np.arange(lo[i], hi[i] + 1e-9, pitch) for i in range(3)]
    shape = tuple(len(a) for a in axes)
    S = np.zeros(shape, np.float32)
    X, Y = np.meshgrid(axes[0], axes[1], indexing="ij")
    for k, z in enumerate(axes[2]):            # one z-slab at a time keeps memory low
        W = np.stack([X.ravel(), Y.ravel(), np.full(X.size, z)], 1)
        acc = np.zeros(len(W), np.float32)
        for T in G:
            Ti = np.linalg.inv(T)
            M = Tinv @ Ti
            idx = W @ M[:3, :3].T + M[:3, 3]
            acc += ndimage.map_coordinates(F, idx.T, order=1, cval=0.0)
        S[:, :, k] = (acc / len(G)).reshape(X.shape)
    log(f"  symmetric field: grid {S.shape}, pitch {pitch:.3f} mm, blur {blur:.2f} mm, {len(G)} group elements")
    return S, np.array([a[0] for a in axes])


def skeleton_graph(S, origin, pitch, prune, log):
    B = S > 0.5
    dt = ndimage.distance_transform_edt(B) * pitch
    K = skeletonize(B)
    idx = np.argwhere(K)
    P = origin + idx * pitch
    key = {tuple(i): k for k, i in enumerate(idx)}
    G = nx.Graph()
    offs = np.argwhere(np.ones((3, 3, 3))) - 1
    for k, i in enumerate(idx):
        G.add_node(k, p=P[k], r=dt[tuple(i)])
        for o in offs:
            j = tuple(i + o)
            if j in key and key[j] > k:
                G.add_edge(k, key[j])
    if G.number_of_nodes() == 0:
        raise RuntimeError("empty skeleton: increase voxel resolution")
    comps = list(nx.connected_components(G))
    G = G.subgraph(max(comps, key=len)).copy()
    for _ in range(6):                     # prune spurs shorter than `prune` x local radius
        removed = 0
        for e in [n for n in G if G.degree(n) == 1]:
            if e not in G:
                continue
            path, cur = [e], e
            while True:
                nb = [x for x in G[cur] if x not in path]
                if len(nb) != 1:
                    break
                cur = nb[0]
                path.append(cur)
                if G.degree(cur) != 2:
                    break
            if G.degree(path[-1]) > 2 and len(path) * pitch < prune * G.nodes[path[-1]]["r"]:
                G.remove_nodes_from(path[:-1])
                removed += 1
        if not removed:
            break
    J = [n for n in G if G.degree(n) > 2]
    cl = list(nx.connected_components(G.subgraph(J)))
    cid = {n: i for i, c in enumerate(cl) for n in c}
    nodes = {f"J{i}": (np.mean([G.nodes[n]["p"] for n in c], 0), max(G.nodes[n]["r"] for n in c))
             for i, c in enumerate(cl)}
    branches, seen = [], set()
    for s in [n for n in G if G.degree(n) != 2]:
        for nb in G[s]:
            if s in cid and nb in cid and cid[s] == cid[nb]:
                continue
            if (s, nb) in seen:
                continue
            path, prev, cur = [s, nb], s, nb
            while G.degree(cur) == 2:
                nxt = [x for x in G[cur] if x != prev][0]
                path.append(nxt)
                prev, cur = cur, nxt
            seen |= {(s, nb), (path[-1], path[-2])}
            ends = []
            for n in (path[0], path[-1]):
                if n in cid:
                    ends.append(f"J{cid[n]}")
                else:
                    nodes[f"E{n}"] = (G.nodes[n]["p"], G.nodes[n]["r"])
                    ends.append(f"E{n}")
            if ends[0] == ends[1] and len(path) < 20:
                continue
            branches.append(dict(a=ends[0], z=ends[1], pts=np.array([G.nodes[n]["p"] for n in path]),
                                 r=np.array([G.nodes[n]["r"] for n in path])))
    log(f"  skeleton: {len(cl)} junctions, {sum(k[0] == 'E' for k in nodes)} free ends, {len(branches)} members")
    return nodes, branches


def detect_holes(mesh, tol, log, max_radius_frac=0.025):
    """Through/blind holes: concave cylinders with good angular coverage."""
    t = time.time()
    diag = float(np.linalg.norm(mesh.extents))
    seg = Segmenter(mesh, tol, kinds=("cylinder",), allow_extrusion=False, min_region_frac=5e-4,
                    samples_per_kind=300, log=lambda s: None, concave_only=True, r_max=max_radius_frac * diag)
    prims = [p for p in seg.run(max_prims=250, patience=6) if p.kind == "cylinder"]
    holes = []
    for p in prims:
        if p.p["r"] > max_radius_frac * diag or len(p.faces) < 6:
            continue
        a, c = p.p["a"], p.p["c"]
        fc = mesh.triangles_center[p.faces]
        w = fc - c
        rad = w - np.outer(w @ a, a)
        rad /= np.linalg.norm(rad, axis=1, keepdims=True)
        concave = np.mean(np.einsum("ij,ij->i", rad, mesh.face_normals[p.faces]) < 0) > 0.8
        u, v = perp_basis(a)
        ang = np.arctan2(rad @ v, rad @ u)
        coverage = len(np.unique(np.floor((ang + np.pi) / (2 * np.pi) * 36))) / 36
        if concave and coverage > 0.65:
            h = (p.pts - c) @ a
            holes.append(dict(c=c + 0.5 * (h.min() + h.max()) * a, a=a, r=float(p.p["r"]),
                              h=float(h.max() - h.min())))
    log(f"  holes: {len(holes)} found ({time.time() - t:.1f}s)")
    return holes


def pad_radius(mesh, hole, pitch):
    """Outer radius of the material ring around a hole (20th percentile over directions)."""
    a, c, r = hole["a"], hole["c"], hole["r"]
    u, v = perp_basis(a)
    rr = np.arange(r + 0.5 * pitch, r + 12 * r, 0.5 * pitch)
    out = []
    for th in np.linspace(0, 2 * np.pi, 48, endpoint=False):
        d = np.cos(th) * u + np.sin(th) * v
        q = c + np.outer(rr, d)
        ins = mesh.contains(q)
        if not ins.any():
            continue
        k = np.argmax(ins)
        j = k + np.argmax(~ins[k:]) if (~ins[k:]).any() else len(rr) - 1
        out.append(rr[j])
    return float(np.percentile(out, 20)) if out else 2 * r


def area_radius(mesh, p, t, r_guess):
    sec = mesh.section(plane_origin=p, plane_normal=t)
    if sec is None:
        return None
    try:
        planar, T = sec.to_2D()
    except Exception:
        return None
    q = trimesh.transform_points([p], np.linalg.inv(T))[0][:2]
    for poly in planar.polygons_full:
        if poly.contains(Point(q)) and poly.area < np.pi * (3 * r_guess) ** 2:
            return float(np.sqrt(poly.area / np.pi))
    return None


def resample(p, r, n=24):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, p[:, i]) for i in range(3)], 1), np.interp(t, s, r)


class Canon:
    """Map a point to its image in a fundamental domain (Dirichlet-style: maximize v . (g x))."""

    def __init__(self, G, spec):
        self.G = G
        self.c = spec.center if spec is not None else np.zeros(3)
        R = spec.axes if spec is not None else np.eye(3)
        self.v = R @ unit(np.array([0.83, 0.47, 0.29]))

    def element(self, x):
        vals = [self.v @ (T[:3, :3] @ x + T[:3, 3] - self.c) for T in self.G]
        return self.G[int(np.argmax(vals))]


# ---- main ------------------------------------------------------------------------
def reconstruct(mesh, opts, spec=None, log=print):
    t0 = time.time()
    diag = float(np.linalg.norm(mesh.extents))
    G = group_matrices(spec) if spec is not None else [np.eye(4)]
    pitch = opts.get("pitch") or max(diag / 450, (np.prod(mesh.extents + 4) / 6e6) ** (1 / 3))
    blur = opts.get("blur") or 3 * pitch
    log(f"smooth mode: {len(G)}-element symmetry group, pitch {pitch:.3f} mm")
    S, origin = symmetric_field(mesh, G, pitch, blur, log)
    nodes, branches = skeleton_graph(S, origin, pitch, opts.get("prune", 1.0), log)

    holes = detect_holes(mesh, opts.get("hole_tol") or 2e-3 * diag, log) if opts.get("holes", True) else []
    for h in holes:
        h["R"] = pad_radius(mesh, h, pitch)
    canon = Canon(G, spec)
    # symmetric hole set: canonicalize, cluster, average, pattern
    reps = []
    for h in holes:
        T = canon.element(h["c"])
        c, a = T[:3, :3] @ h["c"] + T[:3, 3], T[:3, :3] @ h["a"]
        for rp in reps:
            if np.linalg.norm(rp["c"][0] - c) < max(h["r"], 2 * pitch):
                rp["c"].append(c)
                rp["a"].append(a * np.sign(a @ rp["a"][0]))
                for k in ("r", "R", "h"):
                    rp[k].append(h[k])
                break
        else:
            reps.append(dict(c=[c], a=[a], r=[h["r"]], R=[h["R"]], h=[h["h"]]))
    pads = []
    for rp in reps:
        c, a = np.mean(rp["c"], 0), unit(np.mean(rp["a"], 0))
        r, R, h = (float(np.median(rp[k])) for k in ("r", "R", "h"))
        for T in G:
            cc = T[:3, :3] @ c + T[:3, 3]
            if not any(np.linalg.norm(cc - q["c"]) < r for q in pads):
                pads.append(dict(c=cc, a=T[:3, :3] @ a, r=r, R=max(R, 1.3 * r), h=h))
    if pads:
        dias = ", ".join(sorted({"%.2f" % (2 * q["r"]) for q in pads}))
        log(f"  pads: {len(pads)} (hole diameters {dias} mm)")

    def near_pad(x, margin):
        for q in pads:
            w = x - q["c"]
            ax = w @ q["a"]
            if abs(ax) < 0.5 * q["h"] + margin and np.linalg.norm(w - ax * q["a"]) < q["R"] + margin:
                return q["c"]
        return None

    # canonical member classes
    keep_free = opts.get("keep_free_ends", True)
    canon_b = []
    for b in branches:
        if not keep_free and (b["a"][0] == "E" or b["z"][0] == "E"):
            continue
        pts, rad = b["pts"], b["r"]
        if pads and all(near_pad(q, 2 * pitch) is not None for q in pts):
            continue                            # loop around a bolt hole: replaced by its pad
        ends, padflags = [], []
        for key in (b["a"], b["z"]):
            c = near_pad(nodes[key][0], pitch) if pads else None
            ends.append(c if c is not None else nodes[key][0])
            padflags.append(c is not None)
        P = np.vstack([ends[0], pts[2:-2], ends[1]]) if len(pts) > 5 else np.vstack([ends[0], ends[1]])
        Rr = np.r_[rad[0], rad[2:-2], rad[-1]] if len(pts) > 5 else np.r_[rad[0], rad[-1]]
        if np.linalg.norm(np.diff(P, axis=0), axis=1).sum() < 2 * pitch:
            continue
        P, Rr = resample(P, Rr)
        T = canon.element(P[len(P) // 2])
        canon_b.append(dict(p=P @ T[:3, :3].T + T[:3, 3], r=Rr, pad=tuple(padflags),
                            free=(b["a"][0] == "E", b["z"][0] == "E")))
    med_r = float(np.median(np.concatenate([b["r"] for b in branches]))) if branches else pitch
    match_tol = max(3 * pitch, 0.8 * med_r)
    classes = []
    for c in canon_b:
        best = None
        for cl in classes:
            ref = cl[0]["p"]
            for T in G:                      # align the whole curve, not just its midpoint
                q = c["p"] @ T[:3, :3].T + T[:3, 3]
                for rev in (False, True):
                    qq = q[::-1] if rev else q
                    d = np.linalg.norm(qq - ref, axis=1).mean()
                    if d < match_tol and (best is None or d < best[0]):
                        best = (d, cl, qq, rev)
        if best is None:
            classes.append([c])
        else:
            _, cl, qq, rev = best
            flip = (lambda t: t[::-1]) if rev else (lambda t: t)
            cl.append(dict(p=qq, r=flip(c["r"]), pad=flip(c["pad"]), free=flip(c["free"])))
    need = max(1, int(np.ceil(opts.get("min_copy_frac", 0.5) * len(G))))
    kept = [cl for cl in classes if len(cl) >= need]
    log(f"  member classes: {len(kept)} kept of {len(classes)} (need {need} copies of {len(G)}); "
        f"class sizes {sorted([len(c) for c in classes], reverse=True)}")

    # build the fundamental-domain solid
    algebraic = opts.get("algebraic_only", False)
    r_min = opts.get("r_min") or 1.5 * pitch
    smooth_s = opts.get("spline_smoothing", 1.0)
    pieces, joints = [], []
    for cl in kept:
        p = np.mean([c["p"] for c in cl], 0)
        r = np.mean([c["r"] for c in cl], 0)
        for e, i in ((0, 0), (1, -1)):
            if sum(c["pad"][e] for c in cl) * 2 >= len(cl) and pads:
                p[i] = min((q["c"] for q in pads), key=lambda x: np.linalg.norm(x - p[i]))
            else:
                joints.append((p[i].copy(), r[i]))
        s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
        L = s[-1]
        try:
            tck, _ = splprep(p.T, u=s / L, s=smooth_s * len(p) * pitch ** 2, k=3)
        except Exception:
            continue
        n = max(int(L / (6 * med_r) * (2 if not algebraic else 1.2)), 3)
        u = np.linspace(0, 1, n + 1)
        c = np.array(splev(u, tck)).T
        c[0], c[-1] = p[0], p[-1]
        r_fit = np.polyval(np.polyfit(np.linspace(0, 1, len(r)), r, 2), u)
        if opts.get("area_radius", True):
            us = np.linspace(0.15, 0.85, 7)
            cs, ts = np.array(splev(us, tck)).T, np.array(splev(us, tck, der=1)).T
            ok = []
            for ci, ti, ui in zip(cs, ts, us):
                vals = []
                for T in G:
                    v = area_radius(mesh, T[:3, :3] @ ci + T[:3, 3], T[:3, :3] @ unit(ti), np.interp(ui, u, r_fit))
                    if v is not None:
                        vals.append(v)
                if len(vals) >= max(1, len(G) // 2):
                    ok.append((ui, np.median(vals)))
            if len(ok) >= 3:
                uu, rr = np.array(ok).T
                r_fit = np.polyval(np.polyfit(uu, rr, 1 if len(ok) < 5 else 2), u)
        r_fit = np.clip(r_fit, r_min, 4 * med_r)
        try:
            pieces.append(cone_chain(c, r_fit) if algebraic else tube(c, r_fit))
        except Exception as e:
            log(f"  warning: member skipped ({e})")
    merged = []
    for q, rq in joints:
        for m_ in merged:
            if np.linalg.norm(m_[0] - q) < 2 * pitch:
                m_[1] = max(m_[1], rq)
                break
        else:
            merged.append([q, rq])
    pieces += [cf.sphere(2 * float(np.clip(rq, r_min, 4 * med_r))).moved(Location(Vector(*q))) for q, rq in merged]
    pad_solids, cutters = [], []
    seen_pad = []
    for q in pads:                               # pads: one per orbit in the domain, patterned below
        T = canon.element(q["c"])
        cc = T[:3, :3] @ q["c"] + T[:3, 3]
        if any(np.linalg.norm(cc - x) < q["r"] for x in seen_pad):
            continue
        seen_pad.append(cc)
        a = T[:3, :3] @ q["a"]
        base = cc - 0.5 * q["h"] * a
        disc = Solid.makeCylinder(q["R"], q["h"], Vector(*base), Vector(*a))
        try:
            disc = disc.fillet(min(opts.get("pad_fillet", 0.15) * q["h"], 0.3 * q["R"]), disc.Edges())
        except Exception:
            pass
        pad_solids.append(disc)
    for q in pads:
        cutters.append(Solid.makeCylinder(q["r"], q["h"] + 4 * pitch, Vector(*(q["c"] - (0.5 * q["h"] + 2 * pitch) * q["a"])),
                                          Vector(*q["a"])))
    if not pieces and not pad_solids:
        raise RuntimeError("no members reconstructed")
    log(f"  domain: {len(pieces)} members/joints, {len(pad_solids)} pads")

    # Dirichlet fundamental cell of the group: {x : v.(x - c) >= v.(g x - c) for all g}
    if len(G) > 1:
        big = 4 * diag
        cell = Solid.makeBox(2 * big, 2 * big, 2 * big, Vector(*(canon.c - big)))
        for T in G[1:]:
            nrm = canon.v - T[:3, :3].T @ canon.v
            if np.linalg.norm(nrm) < 1e-9:
                continue
            nrm = nrm / np.linalg.norm(nrm)
            u, w = perp_basis(nrm)
            half = Solid.makeBox(2 * big, 2 * big, big, Vector(0, 0, 0))
            M = np.eye(4)
            M[:3, :3] = np.c_[u, w, nrm]
            M[:3, 3] = canon.c - big * u - big * w
            cell = cell * occ_transform(half, M)
        parts = []
        for piece in pad_solids + pieces:
            for T in G:
                img = piece if np.allclose(T, np.eye(4)) else occ_transform(piece, T)
                try:
                    q = img * cell
                except Exception:
                    continue
                if q.Solids() and q.Volume() > 1e-9:
                    parts.append(q)
        wedge = None
        try:                                   # one general fuse of all fragments first
            cand = cf.fuse(*parts) if len(parts) > 1 else parts[0]
            if cand.isValid() and cand.Solids() and cand.Volume() >= max(q.Volume() for q in parts):
                wedge = cand
        except Exception:
            pass
        seq = fuse_all(parts, log, "domain pieces")
        if wedge is None or seq.Volume() > wedge.Volume():
            wedge = seq
        copies = [wedge if np.allclose(T, np.eye(4)) else occ_transform(wedge, T) for T in G]
        try:
            solid = cf.fuse(*copies, glue="full")
            if not (solid.isValid() and solid.Solids() and solid.Volume() > 0.9 * len(G) * wedge.Volume()):
                raise RuntimeError
        except Exception:
            solid = fuse_all(copies, log, "pattern copies")
        log(f"  patterned: fundamental cell from {len(parts)} piece fragments, {len(G)} copies")
    else:
        solid = fuse_all(pad_solids + pieces, log, "pieces")
    clip = opts.get("clip", "thin")
    if clip != "none":
        R = spec.axes if spec is not None else np.eye(3)
        lo, hi = mesh.bounds
        corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])
        fr = corners @ R
        flo, fhi = fr.min(0), fr.max(0)
        big = 10 * diag
        ext = np.where(np.arange(3) == np.argmin(mesh.extents @ np.abs(R)), fhi - flo, big) if clip == "thin" \
            else fhi - flo
        ctr = np.where(np.arange(3) == np.argmin(mesh.extents @ np.abs(R)), 0.5 * (flo + fhi), 0) if clip == "thin" \
            else 0.5 * (flo + fhi)
        box = Solid.makeBox(*ext, Vector(*(-ext / 2)))
        T = np.eye(4)
        T[:3, :3] = R
        T[:3, 3] = R @ ctr
        solid = solid * occ_transform(box, T)
    for cu in cutters:
        solid = solid - cu
    from .brep_helpers import safe_clean
    cleaned = safe_clean(solid)
    if cleaned.isValid() or not solid.isValid():
        solid = cleaned
    sols = solid.Solids()
    if len(sols) > 1:
        solid = max(sols, key=lambda x: x.Volume())
    log(f"  smooth model: valid={solid.isValid()}, volume {solid.Volume():.1f} mm^3 "
        f"(mesh {mesh.volume:.1f}) ({time.time() - t0:.1f}s)")
    res = Result(shape=solid, mode="smooth")
    res.extras = dict(n_members=len(kept), n_pads=len(pads))
    return res
