"""Round wires along helices: parts made of tubes of one circular section swept along helices about one
axis (stents, coil springs, helical lattices).

Recognition, all from the mesh:
  1. wire radius r: rays from surface points along the inward normal cross a round wire through its centre
     line, so the most common chord length is 2 r;
  2. centre-line points: surface points moved inward by r;
  3. axis: the centre-line points lie on one cylinder (radius R) about an axis of the part's frame;
  4. helices: unrolled about the axis, a helix is a straight line theta = phi + w z. For each sign of w the
     rate is found by maximising the sharpness of the histogram of phi = theta - w z; the peaks of that
     histogram are the wires (equally spaced phases are snapped to 360 / n);
  5. ends: a wire ends either on a flat cut normal to the axis (mesh faces facing along the axis at the
     part's end; the wire is trimmed there) or with a cap normal to the wire (the end of its centre-line
     points).

The tubes of crossing wires have equal radii and centre lines that meet, which OCC's Booleans do not fuse
(singular intersection curves); the model is kept as one solid per wire, and it is checked numerically: the
distance from a point to a helix is found by Newton's method, so the deviation both ways and the volume of
the union (Monte Carlo) need no Boolean. SolidWorks fuses the tubes itself.
"""

import math
import time

import numpy as np

MIN_WIRES = 1
CHORD_FRACTION = 0.45      # share of inward chords within 4 % of 2 r for the part to count as round wires
CYL_TOL = 0.15             # centre points on the cylinder: spread of their distance from the axis / r
EXPLAINED = 0.9            # share of centre points the helices must explain


def _frame(a):
    """(a, e1, e2), right-handed; e1 is the first of x, y, z not along a, made normal to a (angles are
    measured from it: from x for an axis along z)."""
    a = np.asarray(a, float) / np.linalg.norm(a)
    g = next(np.eye(3)[k] for k in range(3) if abs(a[k]) < 0.9)
    e1 = g - (g @ a) * a
    e1 /= np.linalg.norm(e1)
    return a, e1, np.cross(a, e1)


def _circle_fit(X):
    """Least-squares circle through 2D points (algebraic, then refined)."""
    A = np.c_[2 * X, np.ones(len(X))]
    b = (X ** 2).sum(1)
    s = np.linalg.lstsq(A, b, rcond=None)[0]
    c = s[:2]
    from scipy.optimize import least_squares
    c = least_squares(lambda c_: np.hypot(*(X - c_).T) - np.hypot(*(X - c_).T).mean(), c).x
    return c, float(np.hypot(*(X - c).T).mean())


def _snap(v, tol, rel=False):
    """v rounded to the coarsest of 1, 2, 3, ... significant digits that moves it by less than tol
    (0 when |v| is below tol)."""
    if v == 0 or (not rel and abs(v) < tol):
        return 0.0
    for d in range(1, 8):
        q = float(f"{v:.{d}g}")
        if abs(q - v) < (tol * abs(v) if rel else tol):
            return q
    return v


def _hist_sharpness(phi, nb=720):
    h = np.bincount(np.floor(phi / (2 * np.pi) * nb).astype(int) % nb, minlength=nb).astype(float)
    return float((h ** 2).sum()), h


def _best_rate(theta, z, sign, w_max):
    """Rate w (sign given) maximising the histogram sharpness of theta - w z; coarse grid, then refined."""
    best = (-1.0, 0.0)
    for w in sign * np.linspace(1e-3 * w_max, w_max, 6000):
        s, _ = _hist_sharpness((theta - w * z) % (2 * np.pi))
        if s > best[0]:
            best = (s, w)
    w0, step = best[1], w_max / 6000
    for _ in range(3):
        for w in np.linspace(w0 - step, w0 + step, 41):
            s, _ = _hist_sharpness((theta - w * z) % (2 * np.pi))
            if s > best[0]:
                best = (s, w)
        w0, step = best[1], step / 20
    return best[1], best[0]


def _peaks(h, frac=0.25):
    """Local maxima of a circular histogram above frac x the highest."""
    nb = len(h)
    hs = np.convolve(np.r_[h[-3:], h, h[:3]], np.ones(7) / 7, mode="same")[3:-3]
    thr = frac * hs.max()
    out = []
    for i in range(nb):
        if hs[i] >= thr and hs[i] >= hs[i - 1] and hs[i] > hs[(i + 1) % nb]:
            out.append(i)
    return out, hs


class Helices:
    """The fitted model: wires on helices about one axis."""

    def __init__(self, c, a, R, r, wires, cuts):
        self.c, (self.a, self.e1, self.e2) = np.asarray(c, float), _frame(a)
        self.R, self.r = float(R), float(r)
        self.wires = wires        # dict(w (rad per unit z), phi (rad at z = 0), z0, z1, cap0, cap1)
        self.cuts = cuts          # (lo, hi): flat cuts normal to the axis (None: no cut)

    def cyl(self, P):
        Q = np.asarray(P, float) - self.c
        z = Q @ self.a
        return np.hypot(Q @ self.e1, Q @ self.e2), np.arctan2(Q @ self.e2, Q @ self.e1), z

    def point(self, wire, z):
        t = wire["phi"] + wire["w"] * z
        return (self.c + np.outer(z, self.a) + self.R * (np.outer(np.cos(t), self.e1)
                                                         + np.outer(np.sin(t), self.e2)))

    def nearest(self, P, wire):
        """For points P: (squared distance to the helix, z of the nearest helix point); Newton on
        f(z) = R^2 + rho^2 - 2 R rho cos(phi + w z - theta) + (z - zp)^2 from z = zp."""
        rho, th, zp = self.cyl(P)
        w, phi, R = wire["w"], wire["phi"], self.R
        z = zp.copy()
        for _ in range(8):
            u = phi + w * z - th
            g = 2 * R * rho * w * np.sin(u) + 2 * (z - zp)
            H = 2 * R * rho * w * w * np.cos(u) + 2
            z = z - g / np.where(np.abs(H) > 1e-12, H, 1e-12)
        u = phi + w * z - th
        return R * R + rho * rho - 2 * R * rho * np.cos(u) + (z - zp) ** 2, z

    def inside(self, P):
        P = np.asarray(P, float)
        ins = np.zeros(len(P), bool)
        for wire in self.wires:
            d2, z = self.nearest(P, wire)
            ins |= (d2 <= self.r ** 2) & (z >= wire["z0"]) & (z <= wire["z1"])
        _, _, zp = self.cyl(P)
        lo, hi = self.cuts
        if lo is not None:
            ins &= zp >= lo
        if hi is not None:
            ins &= zp <= hi
        return ins

    def surface_distance(self, P):
        """Unsigned distance from points to the model surface (tube walls, caps normal to the wire, flat
        cuts), the closest wire taken for each point; accurate near the surface."""
        P = np.asarray(P, float)
        best = np.full(len(P), np.inf)
        _, _, zp = self.cyl(P)
        lo, hi = self.cuts
        for wire in self.wires:
            d2, z = self.nearest(P, wire)
            d = np.sqrt(np.maximum(d2, 0.0))
            wall = np.abs(d - self.r)
            # the end discs: beyond an end the distance to the disc; on the end face itself (within the
            # disc, near the end) the distance to its plane
            for zz, beyond in ((wire["z0"], z < wire["z0"]), (wire["z1"], z > wire["z1"])):
                near = beyond | (np.abs(z - zz) < 2 * self.r)
                if not near.any():
                    continue
                q = self.point(wire, np.full(near.sum(), zz))
                t = self.tangent(wire, np.full(near.sum(), zz))
                v = P[near] - q
                along = np.abs((v * t).sum(1))
                radial = np.linalg.norm(v - (v * t).sum(1)[:, None] * t, axis=1)
                disc = np.hypot(along, np.maximum(radial - self.r, 0.0))
                cur = wall[near]
                wall[near] = np.where(beyond[near], disc, np.minimum(cur, disc))
            # on a flat cut: distance to the cut plane where the point is within the tube
            for zc in (lo, hi):
                if zc is not None:
                    on = d <= self.r
                    wall = np.where(on, np.minimum(wall, np.abs(zp - zc)), wall)
            best = np.minimum(best, wall)
        return best

    def tangent(self, wire, z):
        t = wire["phi"] + wire["w"] * z
        v = (np.outer(np.ones_like(z), self.a) + self.R * wire["w"] * (np.outer(-np.sin(t), self.e1)
                                                                        + np.outer(np.cos(t), self.e2)))
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    def surface_samples(self, n=40000, seed=0):
        """Points on the union's surface: tube walls outside the other tubes and inside the cuts."""
        rng = np.random.default_rng(seed)
        pts = []
        per = max(n // len(self.wires), 100)
        for wire in self.wires:
            z = rng.uniform(wire["z0"], wire["z1"], per)
            psi = rng.uniform(0, 2 * np.pi, per)
            c = self.point(wire, z)
            t = self.tangent(wire, z)
            nrm = c - self.c - np.outer((c - self.c) @ self.a, self.a)
            nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
            b = np.cross(t, nrm)
            pts.append(c + self.r * (np.cos(psi)[:, None] * nrm + np.sin(psi)[:, None] * b))
        P = np.vstack(pts)
        # keep points not strictly inside another tube (shrink the radius slightly for the test)
        keep = np.ones(len(P), bool)
        for wire in self.wires:
            d2, z = self.nearest(P, wire)
            keep &= ~((d2 < (self.r * (1 - 1e-6)) ** 2 * 0.98) & (z > wire["z0"]) & (z < wire["z1"]))
        _, _, zp = self.cyl(P)
        lo, hi = self.cuts
        if lo is not None:
            keep &= zp >= lo
        if hi is not None:
            keep &= zp <= hi
        return P[keep]

    def describe(self):
        hand = {1: "right", -1: "left"}
        fam = {}
        for w in self.wires:
            fam.setdefault(int(np.sign(w["w"])), []).append(w)
        parts = [f"{len(v)} {hand[s]}-handed" for s, v in sorted(fam.items(), reverse=True)]
        p = 2 * math.pi / abs(self.wires[0]["w"])
        return (f"{len(self.wires)} wires ({', '.join(parts)}), wire diameter {2 * self.r:g}, helix radius "
                f"{self.R:g}, pitch {p:g}")


def detect(mesh, tol, log=print, n_samples=30000):
    """Fitted Helices model, or None when the part is not made of round wires on helices."""
    import trimesh
    from trimesh.ray.ray_triangle import RayMeshIntersector
    t0 = time.time()
    P, fi = trimesh.sample.sample_surface(mesh, n_samples, seed=0)
    N = mesh.face_normals[fi]
    eps = 1e-6 * float(np.linalg.norm(mesh.extents))
    O = P - eps * N
    d = np.full(len(P), np.inf)
    rm = RayMeshIntersector(mesh)
    for ch in np.array_split(np.arange(len(P)), max(1, len(P) // 1000)):   # chunks keep memory bounded
        loc, ri, _ = rm.intersects_location(O[ch], -N[ch], multiple_hits=False)
        d[ch[ri]] = np.linalg.norm(loc - O[ch][ri], axis=1)
    fin = d[np.isfinite(d)]
    if len(fin) < 0.5 * len(P):
        return None
    hist, edges = np.histogram(fin, bins=400, range=(0, np.percentile(fin, 90) * 1.5))
    k = int(np.argmax(np.convolve(hist, np.ones(5), "same")))
    chord = 0.5 * (edges[k] + edges[k + 1])
    near = np.abs(d - chord) < 0.04 * chord
    if near.mean() < CHORD_FRACTION:
        return None
    r = float(np.median(d[near])) / 2
    C = P[near] - r * N[near]                    # centre-line points

    # axis: a principal axis of the part about which the centre points lie on one cylinder
    from .symmetry import principal_frame
    try:
        frame = principal_frame(mesh)[1]
    except Exception:
        frame = np.eye(3)
    best = None
    frame = np.asarray(frame, float)
    for a in list(frame) + list(frame.T) + list(np.eye(3)):
        a, e1, e2 = _frame(a)
        X = np.c_[C @ e1, C @ e2]
        c2, R = _circle_fit(X)
        spread = float(np.std(np.hypot(*(X - c2).T)))
        if best is None or spread < best[0]:
            best = (spread, a, e1, e2, c2, R)
    spread, a, e1, e2, c2, R = best
    if spread > CYL_TOL * r or R < 2 * r:
        log(f"  wires: round section r {r:.4g}, but centre lines are not on one cylinder (spread "
            f"{spread:.3g})")
        return None
    c = c2[0] * e1 + c2[1] * e2                 # a point on the axis (z = 0 at the origin's level)
    Q = C - c
    z = Q @ a
    theta = np.arctan2(Q @ e2, Q @ e1)
    L = float(z.max() - z.min())
    w_max = math.tan(math.radians(88)) / R      # helix angle up to 88 degrees from the axis (coils)

    wires = []
    used = np.zeros(len(C), bool)
    zbins = max(int(L / (4 * r)), 4)
    for sign in (1, -1):
        free = ~used
        w, _ = _best_rate(theta[free], z[free], sign, w_max)
        width = 1.5 * r / R * math.sqrt(1 + (R * w) ** 2)    # r across the line, in phase units
        # greedy: strongest line first; its points are removed before the next peak is taken
        fam, taken, first = [], np.zeros(len(C), bool), None
        _, h0 = _hist_sharpness((theta[free] - w * z[free]) % (2 * np.pi))
        h0 = np.convolve(np.r_[h0[-3:], h0, h0[:3]], np.ones(7) / 7, mode="same")[3:-3]
        if np.median(h0) > 0.35 * h0.max():      # no distinct lines (a pipe wall, a filled surface)
            continue
        while True:
            m_ = free & ~taken
            if m_.sum() < 20:
                break
            phi = (theta[m_] - w * z[m_]) % (2 * np.pi)
            _, h = _hist_sharpness(phi)
            hs = np.convolve(np.r_[h[-3:], h, h[:3]], np.ones(7) / 7, mode="same")[3:-3]
            i_ = int(np.argmax(hs))
            if first is None:
                first = hs[i_]
            if hs[i_] < 0.3 * first:
                break
            ph = (i_ + 0.5) / len(h) * 2 * np.pi
            dphi = np.angle(np.exp(1j * (theta - w * z - ph)))
            sel = m_ & (np.abs(dphi) < width)
            ph = (ph + np.mean(dphi[sel])) % (2 * np.pi)
            dphi = np.angle(np.exp(1j * (theta - w * z - ph)))
            sel = m_ & (np.abs(dphi) < width)
            taken |= sel
            cover = len(np.unique(np.floor((z[sel] - z.min()) / L * zbins))) / zbins
            if cover < 0.5:                      # points only near crossings: not a wire
                continue
            fam.append((ph, sel))
        if not fam:
            continue
        # refit the rate with all the family's points (phases free), then snap equal spacing
        A_, b_ = [], []
        for k, (ph, sel) in enumerate(fam):
            dz = z[sel]
            row = np.zeros((sel.sum(), len(fam) + 1))
            row[:, 0] = dz
            row[:, k + 1] = 1.0
            A_.append(row)
            b_.append(ph + np.angle(np.exp(1j * (theta[sel] - w * dz - ph))) + w * dz)
        sol = np.linalg.lstsq(np.vstack(A_), np.concatenate(b_), rcond=None)[0]
        w = float(sol[0])
        phs = np.sort(np.mod(sol[1:], 2 * np.pi))
        n = len(phs)
        step = 2 * np.pi / n
        off = float(np.angle(np.mean(np.exp(1j * n * phs))) / n)
        dev = np.abs(np.angle(np.exp(1j * (phs - off - step * np.round((phs - off) / step)))))
        if n >= 2 and dev.max() * R < 0.5 * r:
            phs = np.mod(off + step * np.arange(n), 2 * np.pi)
        for ph in phs:
            dphi = np.angle(np.exp(1j * (theta - w * z - ph)))
            sel = np.abs(dphi) < width
            used |= sel
            wires.append(dict(w=w, phi=float(ph), zs=z[sel]))
        log(f"  wires: {n} {'right' if w > 0 else 'left'}-handed helices, pitch {2 * math.pi / abs(w):.5g}")
    if len(wires) < MIN_WIRES or used.mean() < EXPLAINED:
        if wires:
            log(f"  wires: helices explain {100 * used.mean():.0f} % of the centre points; not a wire part")
        return None

    # refinement on the mesh vertices (a faceted mesh has its vertices on the true surface): axis point,
    # R, r, and per family the rate and a phase offset, minimising |distance to the nearest helix| - r
    from scipy.optimize import least_squares
    V = mesh.vertices
    zv = (V - c) @ a
    keepv = (zv > zv.min() + 3 * r) & (zv < zv.max() - 3 * r)
    V = V[keepv]
    if len(V) > 15000:
        V = V[np.random.default_rng(0).choice(len(V), 15000, replace=False)]
    fams = sorted({int(np.sign(x["w"])) for x in wires})

    def unpack(q):
        cc = c + q[0] * e1 + q[1] * e2
        RR, rr = q[2], q[3]
        out = []
        for x in wires:
            k = fams.index(int(np.sign(x["w"])))
            out.append(dict(w=x["w"] + q[4 + 2 * k], phi=x["phi"] + q[5 + 2 * k], z0=-1e9, z1=1e9))
        return cc, RR, rr, out

    def resid(q):
        cc, RR, rr, ws = unpack(q)
        mdl = Helices(cc, a, RR, rr, ws, (None, None))
        d2 = np.min([mdl.nearest(V, x)[0] for x in ws], axis=0)
        return np.sqrt(np.maximum(d2, 0)) - rr

    q0 = np.r_[0.0, 0.0, R, r, np.zeros(2 * len(fams))]
    sol = least_squares(resid, q0, loss="soft_l1", f_scale=0.1 * r, x_scale="jac")
    c, R, r, ws = unpack(sol.x)
    for x, y in zip(wires, ws):
        x["w"], x["phi"] = y["w"], y["phi"] % (2 * np.pi)
    rs = resid(sol.x)
    log(f"  wires: refined on {len(V)} vertices: R {R:.5f}, r {r:.5f}; residual median "
        f"{np.median(np.abs(rs)):.2e}, p95 {np.percentile(np.abs(rs), 95):.2e}")

    # common rate per family magnitude: snap R, r, pitch
    R = _snap(R, 0.1 * tol)
    r = _snap(r, 0.1 * tol)
    for s in (1, -1):
        ws = [x for x in wires if np.sign(x["w"]) == s]
        if ws:
            p = _snap(2 * math.pi / abs(np.mean([x["w"] for x in ws])), 1e-4, rel=True)
            for x in ws:
                x["w"] = s * 2 * math.pi / p
    ps = sorted({round(2 * math.pi / abs(x["w"]), 9) for x in wires})
    if len(ps) == 2 and abs(ps[0] - ps[1]) < 1e-3 * ps[0]:
        for x in wires:
            x["w"] = np.sign(x["w"]) * 2 * math.pi / ps[0]
    # axis along a coordinate axis when within 1e-4 rad (angles stay measured from the same e1); the
    # axis point snapped like the sizes
    k = int(np.argmax(np.abs(a)))
    if np.linalg.norm(a - np.sign(a[k]) * np.eye(3)[k]) < 1e-4:
        a = np.sign(a[k]) * np.eye(3)[k]
        a, e1, e2 = _frame(a)
    c = c - (c @ a) * a
    c = np.array([_snap(float(v), 0.1 * tol) if abs(v) > 1e-9 else 0.0 for v in c])
    c = c - (c @ a) * a

    # ends: flat cuts normal to the axis at the part's ends, else caps normal to the wire
    Fz = (mesh.triangles_center - c) @ a
    Fn = mesh.face_normals @ a
    zv = (mesh.vertices - c) @ a
    zlo, zhi = _snap(float(zv.min()), 0.1 * tol), _snap(float(zv.max()), 0.1 * tol)   # the part's ends
    cut_lo = zlo if mesh.area_faces[(Fn < -0.999) & (np.abs(Fz - zlo) < tol)].sum() > 0.25 * math.pi * r * r \
        else None
    cut_hi = zhi if mesh.area_faces[(Fn > 0.999) & (np.abs(Fz - zhi) < tol)].sum() > 0.25 * math.pi * r * r \
        else None
    model = Helices(c, a, R, r, [], (cut_lo, cut_hi))
    tilt = math.atan(R * abs(wires[0]["w"]))
    for x in wires:
        zs = x.pop("zs")
        # cap normal to the wire: its lowest point is r sin(tilt) below the centre
        z0 = zlo + r * math.sin(tilt) if cut_lo is None else zlo - 2 * r
        z1 = zhi - r * math.sin(tilt) if cut_hi is None else zhi + 2 * r
        z0 = max(z0, float(zs.min()) - 2 * r) if cut_lo is None else z0
        z1 = min(z1, float(zs.max()) + 2 * r) if cut_hi is None else z1
        x.update(z0=_snap(z0, 0.1 * tol), z1=_snap(z1, 0.1 * tol))
        z0, z1 = x["z0"], x["z1"]
        # the angle at the wire's start in whole degrees when that is within the fit
        st = math.degrees(x["phi"] + x["w"] * z0)
        if abs(st - round(st)) * math.pi / 180 * R < 0.05 * r:
            x["phi"] = (math.radians(round(st)) - x["w"] * z0) % (2 * math.pi)
        model.wires.append(x)
    log(f"  wires: {model.describe()}; ends: {'cut' if cut_lo is not None else 'caps'} / "
        f"{'cut' if cut_hi is not None else 'caps'} ({time.time() - t0:.1f}s)")
    return model


def check(model, mesh, n_mc=2000000, seed=1):
    """Numerical comparison without Booleans: deviation both ways and the union's volume."""
    import trimesh
    rng = np.random.default_rng(seed)
    S, _ = trimesh.sample.sample_surface(mesh, 30000, seed=seed)
    d_ref = model.surface_distance(S)
    M = model.surface_samples(30000, seed=seed)
    d_mod = np.concatenate([trimesh.proximity.closest_point(mesh, M[ch])[1]
                            for ch in np.array_split(np.arange(len(M)), max(1, len(M) // 2000))])
    # volume: Monte Carlo in the annulus that holds every wire
    R, r = model.R, model.r
    zs = [w["z0"] for w in model.wires] + [w["z1"] for w in model.wires]
    lo, hi = model.cuts
    zlo = max(min(zs) - 2 * r, lo if lo is not None else -np.inf)
    zhi = min(max(zs) + 2 * r, hi if hi is not None else np.inf)
    rho = np.sqrt(rng.uniform((R - 1.2 * r) ** 2, (R + 1.2 * r) ** 2, n_mc))
    th = rng.uniform(0, 2 * np.pi, n_mc)
    zz = rng.uniform(zlo, zhi, n_mc)
    P = model.c + np.outer(zz, model.a) + rho[:, None] * (np.outer(np.cos(th), model.e1)
                                                          + np.outer(np.sin(th), model.e2))
    vol_shell = math.pi * ((R + 1.2 * r) ** 2 - (R - 1.2 * r) ** 2) * (zhi - zlo)
    frac = model.inside(P).mean()
    vol = vol_shell * frac
    se = vol_shell * math.sqrt(frac * (1 - frac) / n_mc)
    return dict(ref_to_model=d_ref, model_to_ref=d_mod, volume=vol, volume_se=se)


def format_check(chk, ref_volume):
    def st(d):
        return (f"mean {np.mean(d):.4f}  median {np.median(d):.4f}  p95 {np.percentile(d, 95):.4f}  "
                f"max {np.max(d):.4f} mm")
    v = chk["volume"]
    return (f"  ref_to_model  {st(chk['ref_to_model'])}\n  model_to_ref  {st(chk['model_to_ref'])}\n"
            f"  volume ref {ref_volume:.4f}  model {v:.4f} +- {chk['volume_se']:.4f} mm^3  "
            f"({100 * (v / ref_volume - 1):+.2f} %)")


def _r(v):
    return repr(float(f"{float(v):.10g}"))


def cadquery_code(model, name):
    """A CadQuery script that builds one solid per wire (tube swept along its helix, trimmed by the
    cuts) and exports them as one STEP compound."""
    w0 = model.wires[0]
    lines = [
        f'"""{name}: {model.describe()}.',
        "One solid per wire: OCC does not fuse tubes of equal radius whose centre lines meet; SolidWorks",
        'does (Combine or merged sweeps)."""',
        "import math",
        "import cadquery as cq",
        "",
        f"AXIS_POINT = ({_r(model.c[0])}, {_r(model.c[1])}, {_r(model.c[2])})",
        f"AXIS = ({_r(model.a[0])}, {_r(model.a[1])}, {_r(model.a[2])})",
        f"E1 = ({_r(model.e1[0])}, {_r(model.e1[1])}, {_r(model.e1[2])})      # angle 0 direction",
        f"R = {_r(model.R)}             # helix radius",
        f"WIRE_R = {_r(model.r)}        # wire radius",
        "CUTS = (" + ", ".join("None" if v is None else _r(v) for v in model.cuts) + ")   # flat cuts normal to the axis (axial coordinate)",
        "WIRES = [  # (pitch, right-handed, angle at the start z0 (deg, from E1), z0, z1)",
    ]
    for w in model.wires:
        st = math.degrees(w["phi"] + w["w"] * w["z0"]) % 360
        lines.append(f"    ({_r(2 * math.pi / abs(w['w']))}, {w['w'] > 0}, {_r(st)}, "
                     f"{_r(w['z0'])}, {_r(w['z1'])}),")
    lines += [
        "]",
        "",
        "",
        "def frame():",
        "    a = cq.Vector(*AXIS).normalized()",
        "    e1 = cq.Vector(*E1).normalized()",
        "    return a, e1, a.cross(e1)",
        "",
        "",
        "def wire_solid(pitch, right, start_deg, z0, z1):",
        "    a, e1, e2 = frame()",
            "    phi0 = math.radians(start_deg)                   # angle at the start of the wire",
        "    h = cq.Wire.makeHelix(pitch, z1 - z0, R, lefthand=not right)   # local: starts at (R, 0, 0), up z",
        "    # local frame -> part: x along the start angle, z along the axis, origin at the axis point + z0",
        "    xd = e1 * math.cos(phi0) + e2 * math.sin(phi0)",
        "    loc = cq.Location(cq.Plane(origin=cq.Vector(*AXIS_POINT) + a * z0, xDir=xd, normal=a))",
        "    h = h.moved(loc)",
        "    prof = cq.Wire.makeCircle(WIRE_R, h.startPoint(), h.tangentAt(0))",
        "    s = cq.Solid.sweep(prof, [], h, True, True)",
        "    lo, hi = CUTS",
        "    if lo is not None or hi is not None:",
        "        big = 10 * (R + pitch + abs(z1 - z0))",
        "        zl = lo if lo is not None else z0 - big",
        "        zh = hi if hi is not None else z1 + big",
        "        box = cq.Solid.makeBox(2 * big, 2 * big, zh - zl, pnt=cq.Vector(-big, -big, zl),",
        "                               dir=cq.Vector(0, 0, 1))",
        "        box = box.moved(cq.Location(cq.Plane(origin=cq.Vector(*AXIS_POINT), xDir=e1, normal=a)))",
        "        s = s.intersect(box)",
        "    return s",
        "",
        "",
        "def build():",
        "    return cq.Compound.makeCompound([wire_solid(*w) for w in WIRES])",
        "",
        "",
        "result = build()",
        "",
        'if __name__ == "__main__":',
        "    import sys",
        f"    out = sys.argv[1] if len(sys.argv) > 1 else {name + '.step'!r}",
        "    cq.exporters.export(result, out)",
        "    print('solids', len(result.Solids()), 'sum of volumes', round(sum(s.Volume() for s in result.Solids()), 6))",
    ]
    return "\n".join(lines) + "\n"


def reconstruct(mesh, tol, log=print):
    """Result with one solid per wire (a compound) and the numerical check, or None."""
    from .core import Result
    model = detect(mesh, tol, log=log)
    if model is None:
        return None
    ns = {}
    code = cadquery_code(model, "wires")
    t0 = time.time()
    exec(compile(code.replace('if __name__ == "__main__":', "if False:"), "wires_script", "exec"), ns)
    shape = ns["result"]
    chk = check(model, mesh)
    log(f"  wires: model built ({len(shape.Solids())} solids, {time.time() - t0:.1f}s); numerical check:\n"
        + format_check(chk, float(mesh.volume)))
    vol_err = abs(chk["volume"] / mesh.volume - 1)
    p95 = max(np.percentile(chk["ref_to_model"], 95), np.percentile(chk["model_to_ref"], 95))
    if vol_err > 0.02 or p95 > 2 * tol:
        log(f"  wires: rejected (volume {100 * vol_err:.2f} %, deviation p95 {p95:.4f} mm)")
        return None
    res = Result(shape=shape, mode="exact")
    res.extras = dict(wires=model, wires_check=chk, volume_error=vol_err, valid=True,
                      passed=bool(vol_err < 0.01 and p95 < tol))
    return res


def _volume_mc(model, cuts, n=3000000, seed=7):
    """Monte Carlo volume of the union of the wires with the given cuts (mm^3, standard error)."""
    m = Helices(model.c, model.a, model.R, model.r, model.wires, cuts)
    rng = np.random.default_rng(seed)
    R, r = model.R, model.r
    zl = min(w["z0"] for w in model.wires) - 2 * r
    zh = max(w["z1"] for w in model.wires) + 2 * r
    rho = np.sqrt(rng.uniform((R - 1.2 * r) ** 2, (R + 1.2 * r) ** 2, n))
    th = rng.uniform(0, 2 * np.pi, n)
    zz = rng.uniform(zl, zh, n)
    P = m.c + np.outer(zz, m.a) + rho[:, None] * (np.outer(np.cos(th), m.e1) + np.outer(np.sin(th), m.e2))
    shell = math.pi * ((R + 1.2 * r) ** 2 - (R - 1.2 * r) ** 2) * (zh - zl)
    f = m.inside(P).mean()
    return shell * f, shell * math.sqrt(f * (1 - f) / n)


def _families(model):
    """Wires grouped by hand, pitch and extent; a group with equally spaced start angles is one seed wire and
    a circular pattern."""
    groups = {}
    for w in model.wires:
        st = (math.degrees(w["phi"] + w["w"] * w["z0"])) % 360
        key = (w["w"] > 0, round(2 * math.pi / abs(w["w"]), 9), round(w["z0"], 9), round(w["z1"], 9))
        groups.setdefault(key, []).append(st)
    out = []
    for (right, pitch, z0, z1), sts in groups.items():
        sts = sorted(sts)
        n = len(sts)
        equal = n > 1 and all(abs(((sts[k] - sts[0]) - 360 * k / n + 180) % 360 - 180) < 1e-6 for k in range(n))
        out.append(dict(right=right, pitch=pitch, z0=z0, z1=z1, starts=sts, pattern=equal))
    return out


def solidworks_code(model, name):
    """SolidWorks script (pywin32): per family of wires a base circle on a plane at the wires' start, a helix
    (InsertHelix: height and pitch, start angle), a sweep with a circular profile (CreateDefinition(swFmSweep)),
    and a circular pattern about the axis when the family's start angles are equally spaced; then the flat
    cuts (FeatureCut4 from an offset). Calls checked with sw_probe_helix.py (SolidWorks 2025)."""
    a = model.a
    k = int(np.argmax(np.abs(a)))
    if np.linalg.norm(a - np.sign(a[k]) * np.eye(3)[k]) > 1e-9:
        raise ValueError("wires: the axis is not along a coordinate axis")
    # build frame: z along the axis, x along e1 (angles from x); the part is built in this frame
    T = np.array([model.e1, model.e2, model.a])
    c = T @ model.c
    v_cut, se_cut = _volume_mc(model, model.cuts)
    v_nocut, se_nocut = _volume_mc(model, (None, None))
    v_top, _ = _volume_mc(model, (None, model.cuts[1]))
    fams = _families(model)
    cuts = [None if v is None else float(v) for v in model.cuts]
    data = [
        f'"""{name}: {model.describe()}.',
        "SolidWorks script (pywin32), generated by stl2cad (wires). Feature tree per family of wires: base",
        "circle on a plane at the wires' start, Helix (height and pitch, start angle), Sweep with a circular",
        "profile, Circular pattern about the axis; then the flat cuts. Every wire's start and end are checked",
        "against the sweep's end faces; the part volume is checked against stl2cad's numerical model.",
        "The part is built with the wire axis along z (Front Plane normal), angles from x." +
        ("" if np.allclose(T, np.eye(3)) else " The STL frame is rotated: build = T @ stl, T = " +
         repr(np.round(T, 12).tolist())) + '"""',
        "import math",
        "import sys",
        "import time",
        "",
        "import pythoncom",
        "import win32com.client",
        "from win32com.client import VARIANT",
        "",
        "NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)",
        "MM = 0.001",
        f"AXIS = ({_r(c[0])}, {_r(c[1])})          # axis point (x, y), mm",
        f"R = {_r(model.R)}                    # helix radius, mm",
        f"WIRE_D = {_r(2 * model.r)}               # wire diameter, mm",
        f"CUTS = ({', '.join('None' if v is None else _r(v) for v in cuts)})    # flat cuts below / above (z, mm)",
        "FAMILIES = [  # right-handed, pitch, z0, z1, start angles (deg from x), circular pattern",
    ]
    for f in fams:
        data.append(f"    dict(right={f['right']}, pitch={_r(f['pitch'])}, z0={_r(f['z0'])}, z1={_r(f['z1'])}, "
                    f"starts={[_r(s) for s in f['starts']]!r}, pattern={f['pattern']}),".replace("'", ""))
    data += [
        "]",
        f"VOLUME = {v_cut:.6f}          # expected part volume (mm^3), stl2cad's model, +- {se_cut:.4f}",
        f"VOLUME_NO_CUTS = {v_nocut:.6f}  # before the cuts, +- {se_nocut:.4f}",
        f"VOLUME_TOP_CUT = {v_top:.6f}  # after the top cut only",
        f"VOLUME_TOL = {max(3 * max(se_cut, se_nocut), 0.002 * v_cut):.6f}",
    ]
    return "\n".join(data) + "\n" + _SW_BODY


_SW_BODY = r'''

def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


def settle(s=0.1):
    time.sleep(s)


def sw_constant(app, enum, member):
    import os
    exe = app.GetExecutablePath
    exe = str(exe() if callable(exe) else exe)
    for f in (exe, os.path.dirname(exe)):
        for sub in ((), ("api",)):
            path = os.path.join(f, *sub, "swconst.tlb")
            if not os.path.exists(path):
                continue
            tl = pythoncom.LoadTypeLib(path)
            for i in range(tl.GetTypeInfoCount()):
                if tl.GetDocumentation(i)[0] != enum:
                    continue
                ti = tl.GetTypeInfo(i)
                for j in range(ti.GetTypeAttr().cVars):
                    vd = ti.GetVarDesc(j)
                    if ti.GetNames(vd.memid)[0] == member:
                        return int(vd.value)
    raise LookupError(f"{enum}.{member} not found")


class SW:
    def __init__(self):
        self.app = win32com.client.Dispatch("SldWorks.Application")
        self.app.Visible = True
        try:
            self.app.SetUserPreferenceToggle(10, False)     # no Modify dialog when a dimension is created
        except Exception:
            pass
        self.model = call0(self.app, "NewPart")
        if self.model is None:
            self.model = self.app.NewDocument(self.app.GetUserPreferenceStringValue(8), 0, 0, 0)
        try:
            self.fm_sweep = sw_constant(self.app, "swFeatureNameID_e", "swFmSweep")
        except Exception:
            self.fm_sweep = 17                               # SolidWorks 2025 (probe)
        try:
            self.model.Extension.SetUserPreferenceInteger(
                sw_constant(self.app, "swUserPreferenceIntegerValue_e", "swUnitSystem"), 0,
                sw_constant(self.app, "swUnitSystem_e", "swUnitSystem_MMGS"))
        except Exception as ex:
            print(f"units: not set ({type(ex).__name__})")
        self.dims_ok, self.dims_failed, self.status = 0, [], []
        self.quiet(True)

    def quiet(self, on):
        """No graphics or feature-tree updates while building (faster); restored at the end."""
        for obj, attr in ((lambda: self.app, "CommandInProgress"),
                          (lambda: self.model.FeatureManager, "EnableFeatureTree"),
                          (lambda: call0(self.model, "ActiveView"), "EnableGraphicsUpdate")):
            try:
                setattr(obj(), attr, on if attr == "CommandInProgress" else not on)
            except Exception:
                pass

    def last(self):
        return self.model.FeatureByPositionReverse(0)

    def mass(self):
        self.model.ForceRebuild3(False)
        mp = call0(self.model.Extension, "CreateMassProperty")
        if mp is None:
            return 0.0, (0.0, 0.0, 0.0)
        return mp.Volume / MM ** 3, tuple(c / MM for c in mp.CenterOfMass)

    def delete(self, *feats):
        for f in feats:
            if f is None:
                continue
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(1)            # with absorbed features (the base sketch)

    def plane(self, name, base, dist, point):
        """base (by name) itself, or a reference plane at distance dist (mm) from it through point (mm),
        its side checked through a sketch's transform."""
        sel = lambda: self.model.Extension.SelectByID2(base, "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
        if abs(dist) < 1e-12:
            f = self.model.FeatureByName(base)
            return f
        for flip in (0, 256):                                  # swRefPlaneReferenceConstraint_OptionFlip
            self.model.ClearSelection2(True)
            ok = sel()
            p = self.model.FeatureManager.InsertRefPlane(8 | flip, abs(dist) * MM, 0, 0, 0, 0)
            if p is None:
                print(f"  {name}: InsertRefPlane (flag {8 | flip}, base {base} selected {ok}) returned None")
                continue
            f = self.last()
            T = self.sketch_transform(f)
            z = None if T is None else to_sk(T, point)[2]
            if z is not None and abs(z) < 1e-7:
                f.Name = name
                return f
            print(f"  {name}, flag {8 | flip}: " + ("no sketch transform" if T is None else
                                                    f"plane {z / MM:.6f} mm off the point; other side"))
            self.delete(f)
        raise RuntimeError(f"could not make {name}")

    def open_sketch(self, plane):
        """A sketch on the plane: (sketch manager, sketch, transform). The first sketch call after NewPart
        opened nothing in SolidWorks 2025 (swframe run): the plane is selected and the call repeated until a
        sketch is active."""
        sm = self.model.SketchManager
        for attempt in range(4):
            self.model.ClearSelection2(True)
            plane.Select2(False, 0)
            sm.InsertSketch(True)
            settle()
            sk = call0(sm, "ActiveSketch")
            if sk is None:
                sk = call0(self.model, "GetActiveSketch2")
            if sk is not None:
                if attempt:
                    print(f"  sketch opened on try {attempt + 1}")
                return sm, sk, list(sk.ModelToSketchTransform.ArrayData)
        raise RuntimeError(f"no sketch opened on {plane.Name} (4 tries)")

    def sketch_transform(self, plane):
        try:
            sm, sk, T = self.open_sketch(plane)
        except RuntimeError as ex:
            print(f"  {ex}")
            return None
        sm.InsertSketch(True)
        f = self.last()
        if call0(f, "GetTypeName2") == "ProfileFeature":
            self.delete(f)
        return T

    def circle_sketch(self, plane, centre, radius, z, name):
        """Circle (centre (x, y) at level z, mm) on the plane, with its diameter and its centre's distances
        from the origin (or a relation to the origin where a distance is zero)."""
        sm, sk, T = self.open_sketch(plane)
        u, v, _ = to_sk(T, (centre[0], centre[1], z))
        o = to_sk(T, (0.0, 0.0, z))
        sm.AddToDB = True
        g = sm.CreateCircleByRadius(u, v, 0.0, radius * MM)
        sm.AddToDB = False
        ext = self.model.Extension
        self.model.ClearSelection2(True)
        d = self.model.AddDimension2(u + 1.5 * radius * MM, v + 1.5 * radius * MM, 0.0) if g.Select4(False, NOTHING) else None
        self.count(d, f"{name} diameter")
        c = call0(g, "GetCenterPoint2")
        on_axis = [abs(u - o[0]) < 1e-12, abs(v - o[1]) < 1e-12]
        if all(on_axis):
            self.model.ClearSelection2(True)
            if c.Select4(False, NOTHING) and ext.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0, True, 0, NOTHING, 0):
                self.model.SketchAddConstraints("sgCOINCIDENT")
        else:
            for kind, zero in (("H", on_axis[0]), ("V", on_axis[1])):
                self.model.ClearSelection2(True)
                ok = c.Select4(False, NOTHING) and ext.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0,
                                                                   True, 0, NOTHING, 0)
                if not ok:
                    self.dims_failed.append(f"{name} centre {kind} (selection)")
                    continue
                if zero:
                    self.model.SketchAddConstraints("sgVERTICALPOINTS2D" if kind == "H" else "sgHORIZONTALPOINTS2D")
                    continue
                d = (self.model.AddHorizontalDimension2((u + o[0]) / 2, v - 2 * radius * MM, 0.0) if kind == "H"
                     else self.model.AddVerticalDimension2(u - 2 * radius * MM, (v + o[1]) / 2, 0.0))
                self.count(d, f"{name} centre {kind}")
        try:
            st = call0(sk, "GetConstrainedStatus")
        except Exception:
            st = None
        sm.InsertSketch(True)
        f = self.last()
        f.Name = name
        self.status.append({2: "under defined", 3: "fully defined", 4: "over defined"}.get(st, str(st)))
        return f

    def count(self, d, what):
        if d is None:
            self.dims_failed.append(what)
        else:
            self.dims_ok += 1

    def planar_centres(self, feature=None):
        """Centres (face box) of the planar faces of the feature (all bodies' faces if it gives none)."""
        out = []
        faces = []
        if feature is not None:
            try:
                faces = list(call0(feature, "GetFaces") or [])
            except Exception:
                faces = []
        if not faces:
            try:
                for b in (self.model.GetBodies2(0, True) or []):
                    faces += list(call0(b, "GetFaces") or [])
            except Exception:
                pass
        for f in faces:
            try:
                if call0(call0(f, "GetSurface"), "IsPlane"):
                    bx = call0(f, "GetBox")
                    out.append(tuple((bx[k] + bx[k + 3]) / 2 / MM for k in range(3)))
            except Exception:
                pass
        return out

    def wire(self, fam, start, k, plane):
        """Base circle, helix and sweep of one wire; its end faces must sit where the wire starts and ends.
        Variants (helix reversed, start angle mirrored) are tried in turn."""
        H = fam["z1"] - fam["z0"]
        th = math.radians(start)
        turn = 2 * math.pi * H / fam["pitch"] * (1 if fam["right"] else -1)
        S = (AXIS[0] + R * math.cos(th), AXIS[1] + R * math.sin(th), fam["z0"])
        E = (AXIS[0] + R * math.cos(th + turn), AXIS[1] + R * math.sin(th + turn), fam["z1"])
        fm = self.model.FeatureManager
        for rev, ang in ((False, th), (True, th), (False, -th), (True, -th), (False, math.pi - th),
                         (True, math.pi - th)):
            sk = self.circle_sketch(plane, AXIS, R, fam["z0"], f"Wire {k} circle")
            self.model.ClearSelection2(True)
            sk.Select2(False, 0)
            self.model.InsertHelix(rev, not fam["right"], False, False, 2, H * MM, fam["pitch"] * MM,
                                   H / fam["pitch"], 0.0, ang % (2 * math.pi))
            hx = self.last()
            if call0(hx, "GetTypeName2") != "Helix":
                raise RuntimeError(f"wire {k}: no helix")
            hx.Name = f"Wire {k} helix"
            self.model.ClearSelection2(True)
            hx.Select2(True, 4)
            d = fm.CreateDefinition(self.fm_sweep)
            d.CircularProfile = True
            d.CircularProfileDiameter = WIRE_D * MM
            sw_ = fm.CreateFeature(d)
            if sw_ is None:
                raise RuntimeError(f"wire {k}: no sweep")
            sw_.Name = f"Wire {k}"
            P = self.planar_centres(sw_)
            # an end cap overlapped by another wire's (wires starting at one node) is trimmed by the merge, so
            # its centre may move by up to about the wire radius; a wrong variant is off by far more
            near = lambda q: any(math.dist(p, q) < WIRE_D for p in P)
            if near(S) and near(E):
                return sw_
            print(f"  wire {k}: reversed={rev}, start angle {math.degrees(ang):.3f} deg: ends not at the wire's ends; "
                  "next variant")
            self.status.pop()
            self.delete(sw_, hx)
        raise RuntimeError(f"wire {k}: no helix variant puts the wire's ends in place")

    def axis(self):
        """Reference axis through the axis point, normal to the Front Plane: the Top and Right Planes, or
        planes offset from them."""
        px = self.plane("Axis plane x", "Right Plane", AXIS[0], (AXIS[0], 0.0, 0.0))
        py = self.plane("Axis plane y", "Top Plane", AXIS[1], (0.0, AXIS[1], 0.0))
        self.model.ClearSelection2(True)
        px.Select2(False, 0)
        py.Select2(True, 0)
        self.model.InsertAxis2(True)
        ax = self.last()
        if call0(ax, "GetTypeName2") != "RefAxis":
            raise RuntimeError("no axis")
        ax.Name = "Wire axis"
        return ax

    def pattern(self, ax, seeds, n, added, name):
        """Circular pattern of the seed wires about the axis; selection variants are tried in turn. The copies
        must add at least half of the volume expected from them (added, mm^3). None if SolidWorks makes no
        such pattern (a pattern whose copies would not touch the part makes none)."""
        fm = self.model.FeatureManager
        ext = self.model.Extension
        v0 = self.mass()[0]

        def pick_seeds(mark, append):
            ok = True
            for i_, sd in enumerate(seeds):
                ok = sd.Select2(append or i_ > 0, mark) and ok
            return ok
        variants = [
            ("axis Select2 mark 1, seeds mark 4", lambda: ax.Select2(False, 1) and pick_seeds(4, True)),
            ("axis by name (AXIS) mark 1, seeds mark 4",
             lambda: ext.SelectByID2(ax.Name, "AXIS", 0, 0, 0, False, 1, NOTHING, 0) and pick_seeds(4, True)),
            ("seeds mark 4 first, then axis mark 1", lambda: pick_seeds(4, False) and ax.Select2(True, 1)),
        ]
        calls = [("FeatureCircularPattern4", lambda: fm.FeatureCircularPattern4(n, 2 * math.pi, False, "NULL",
                                                                               False, True, False)),
                 ("FeatureCircularPattern5", lambda: fm.FeatureCircularPattern5(n, 2 * math.pi, False, "NULL",
                                                                               False, True, False, False, False,
                                                                               False, 1, 0.0, "NULL", False))]
        for label, sel in variants:
            for cname, call in calls:
                self.model.ClearSelection2(True)
                try:
                    ok = sel()
                    nsel = self.model.SelectionManager.GetSelectedObjectCount2(-1)
                    f = call()
                except Exception as ex:
                    print(f"  pattern ({label}, {cname}): {type(ex).__name__}: {ex}")
                    continue
                if f is None:
                    print(f"  pattern ({label}, {cname}): selected {ok} ({nsel} objects); no feature")
                    continue
                v = self.mass()[0]
                if v - v0 >= 0.5 * added:
                    f.Name = name
                    return f
                print(f"  pattern ({label}, {cname}): feature, but volume {v:.5f} mm^3 (before {v0:.5f}, expected "
                      f"about {added:.5f} more); removed")
                self.delete(f)
        return None

    def cut(self, z, up, name, expected):
        """Cut everything beyond level z (above if up): circle sketch on the Front Plane, Cut-Extrude from
        an offset of z; the direction giving the expected volume is kept."""
        big = R + 2 * WIRE_D
        sk = self.circle_sketch(self.model.FeatureByName("Front Plane"), AXIS, big, 0.0, f"{name} sketch")
        self.status.pop()
        depth = 2 * R + 10 * WIRE_D + max(abs(f["z1"] - f["z0"]) for f in FAMILIES)
        fm = self.model.FeatureManager
        for rev in (True, False):
            for flip in (False, True):
                self.model.ClearSelection2(True)
                sk.Select2(False, 0)
                f = fm.FeatureCut4(True, False, rev, 0, 0, depth * MM, 0.0, False, False, False, False,
                                   0.0, 0.0, False, False, False, False, False, True, True,
                                   False, False, False, 3, abs(z) * MM, flip, False)
                if f is None:
                    continue
                v, c = self.mass()
                print(f"  {name}: reverse={rev}, flip_start={flip}: volume {v:.5f} mm^3 (expected {expected:.5f})")
                if abs(v - expected) <= VOLUME_TOL:
                    f.Name = name
                    return v
                self.model.ClearSelection2(True)
                f.Select2(False, 0)
                self.model.Extension.DeleteSelection2(0)
        raise RuntimeError(f"{name}: no direction gives the expected volume")


def hide_construction(model):
    """Hide the reference planes, axes and sketches the script made (BlankRefGeom, BlankSketch), so that only
    the solid shows. They are hidden, not suppressed: suppressing a plane would suppress what is built on it."""
    n = {"planes": 0, "axes": 0, "sketches": 0}
    base = {"Front Plane", "Top Plane", "Right Plane"}
    f = call0(model, "FirstFeature")
    feats = []
    while f is not None and len(feats) < 100000:
        feats.append(f)
        f = call0(f, "GetNextFeature")
    for f in feats:
        kind = call0(f, "GetTypeName2")
        try:
            if kind in ("RefPlane", "RefAxis") and f.Name not in base:
                model.ClearSelection2(True)
                f.Select2(False, 0)
                model.BlankRefGeom()
                n["planes" if kind == "RefPlane" else "axes"] += 1
            elif kind in ("ProfileFeature", "3DProfileFeature"):
                model.ClearSelection2(True)
                f.Select2(False, 0)
                model.BlankSketch()
                n["sketches"] += 1
        except Exception:
            pass
        # sketches absorbed by a feature (a sweep's path, a loft's sections) are its sub-features
        sub = call0(f, "GetFirstSubFeature")
        guard = 0
        while sub is not None and guard < 1000:
            guard += 1
            try:
                if call0(sub, "GetTypeName2") in ("ProfileFeature", "3DProfileFeature"):
                    model.ClearSelection2(True)
                    sub.Select2(False, 0)
                    model.BlankSketch()
                    n["sketches"] += 1
            except Exception:
                pass
            sub = call0(sub, "GetNextSubFeature")
    model.ClearSelection2(True)
    print(f"hidden: {n['planes']} plane(s), {n['axes']} axis(es), {n['sketches']} sketch(es)")


def wire_volume(fam):
    """Volume of one wire of the family (mm^3): circle area x helix length."""
    H = fam["z1"] - fam["z0"]
    return math.pi * (WIRE_D / 2) ** 2 * H * math.sqrt(1 + (2 * math.pi * R / fam["pitch"]) ** 2)


def to_sk(T, p):
    """Model point (mm) -> sketch coordinates (m), by the sketch's ModelToSketchTransform."""
    x = [c * MM for c in p]
    return tuple((sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3))


def build():
    sw = SW()
    print(f"connected to SolidWorks; building {sum(len(f['starts']) for f in FAMILIES)} wire(s) in "
          f"{len(FAMILIES)} famil{'y' if len(FAMILIES) == 1 else 'ies'}")
    # one seed wire per family; the seeds of the families with equally spaced wires are patterned together
    # (a pattern of one family alone: its copies need not touch the part, and SolidWorks makes none)
    planes, seeds = {}, {}
    for k, fam in enumerate(FAMILIES, 1):
        same = next((planes[j] for j in planes if FAMILIES[j - 1]["z0"] == fam["z0"]), None)
        planes[k] = same or sw.plane(f"Wires start z {fam['z0']:g}", "Front Plane", fam["z0"],
                                     (AXIS[0], AXIS[1], fam["z0"]))
        st = fam["starts"][0]
        seeds[k] = sw.wire(fam, st, f"{k}", planes[k])
        print(f"  wire {k}: {'right' if fam['right'] else 'left'}-handed, pitch {fam['pitch']:g}, start {st:g} deg; "
              f"volume so far {sw.mass()[0]:.5f} mm^3")
    left = {k: list(fam["starts"][1:]) for k, fam in enumerate(FAMILIES, 1)}
    ax = None
    by_count = {}
    for k, fam in enumerate(FAMILIES, 1):
        if fam["pattern"]:
            by_count.setdefault(len(fam["starts"]), []).append(k)
    for n, ks in by_count.items():
        ax = ax or sw.axis()
        groups = [ks] + ([[k] for k in ks] if len(ks) > 1 else [])
        for g in groups:
            if any(not left[k] for k in g):
                continue
            added = (n - 1) * sum(wire_volume(FAMILIES[k - 1]) for k in g)
            name = "Wires pattern" if len(g) > 1 else f"Wires {g[0]} pattern"
            if sw.pattern(ax, [seeds[k] for k in g], n, added, name) is not None:
                print(f"  circular pattern x {n} of famil{'ies' if len(g) > 1 else 'y'} {', '.join(map(str, g))}: "
                      f"volume so far {sw.mass()[0]:.5f} mm^3")
                for k in g:
                    left[k] = []
    for k, fam in enumerate(FAMILIES, 1):
        if left[k] and fam["pattern"]:
            print(f"  family {k}: no circular pattern; the other {len(left[k])} wires one by one")
        for j, st in enumerate(left[k], 2):
            tag = f"{k}.{j}"
            sw.wire(fam, st, tag, planes[k])
            print(f"  wire {tag}: {'right' if fam['right'] else 'left'}-handed, pitch {fam['pitch']:g}, start "
                  f"{st:g} deg")
    v, c = sw.mass()
    print(f"all wires: volume {v:.5f} mm^3 (expected {VOLUME_NO_CUTS:.5f} +- {VOLUME_TOL:.5f})")
    if abs(v - VOLUME_NO_CUTS) > VOLUME_TOL:
        print("WARNING: the wires' volume differs from stl2cad's model")
    lo, hi = CUTS
    if hi is not None:
        v = sw.cut(hi, True, "Cut top", VOLUME_TOP_CUT)
    if lo is not None:
        v = sw.cut(lo, False, "Cut bottom", VOLUME)
    v, c = sw.mass()
    print(f"part: volume {v:.5f} mm^3 (expected {VOLUME:.5f} +- {VOLUME_TOL:.5f}), centre of mass "
          f"({c[0]:.4f}, {c[1]:.4f}, {c[2]:.4f}) mm")
    print(f"dimensions: {sw.dims_ok} added" + (f", failed: {sw.dims_failed}" if sw.dims_failed else ""))
    print(f"sketch status: {', '.join(sw.status)}")
    try:
        sw.app.SetUserPreferenceToggle(10, True)
    except Exception:
        pass
    try:
        hide_construction(sw.model)
    except Exception as ex:
        print(f"construction geometry not hidden: {type(ex).__name__}: {ex}")
    sw.quiet(False)
    try:
        sw.model.GraphicsRedraw2()
    except Exception:
        pass
    return sw


if __name__ == "__main__":
    sw = build()
    if len(sys.argv) > 1:
        sw.model.SaveAs3(sys.argv[1], 0, 0)
'''
