"""Analytic surface primitives: plane, cylinder, cone, sphere, torus, and the
non-algebraic extrusion (a spline curve swept along a direction).

Each primitive provides signed distance, unit normal, a minimal-sample
constructor (for RANSAC), a least-squares refit on points, and a rigid or
reflective transform.
"""

import numpy as np
from scipy.optimize import least_squares

ALGEBRAIC = ("plane", "cylinder", "cone", "sphere", "torus")
COMPLEXITY = {"plane": 0, "cylinder": 1, "cone": 2, "sphere": 2, "torus": 3, "extrusion": 4, "revolution": 4,
              "facets": 5}


def unit(v):
    return v / np.linalg.norm(v)


def perp_basis(a):
    t = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    u = unit(np.cross(a, t))
    return u, np.cross(a, u)


class Prim:
    def __init__(self, kind, **p):
        self.kind = kind
        self.p = {k: (np.asarray(v, float) if isinstance(v, (list, tuple, np.ndarray)) else v)
                  for k, v in p.items()}
        self.faces = np.zeros(0, int)   # mesh face indices of the region
        self.pts = np.zeros((0, 3))     # region vertices (used for extents and refits)
        self.rms = np.inf
        self.area = 0.0

    def __repr__(self):
        q = {k: (np.round(v, 4).tolist() if isinstance(v, np.ndarray) else round(float(v), 5))
             for k, v in self.p.items() if k not in ("curve",)}
        return f"{self.kind}({q}, area={self.area:.3f}, rms={self.rms:.2e})"

    # ---- geometry -------------------------------------------------------
    def distance(self, X):
        k, p = self.kind, self.p
        if k == "plane":
            return X @ p["n"] - p["d"]
        if k in ("cylinder", "cone", "torus"):
            w = X - (p["c"] if k != "cone" else p["v"])
            h = w @ p["a"]
            rho = np.linalg.norm(w - np.outer(h, p["a"]), axis=1)
            if k == "cylinder":
                return rho - p["r"]
            if k == "cone":
                return rho * np.cos(p["t"]) - h * np.sin(p["t"])
            return np.hypot(rho - p["R"], h) - p["r"]
        if k == "sphere":
            return np.linalg.norm(X - p["c"], axis=1) - p["r"]
        if k == "extrusion":
            return _extrusion_distance(self, X)
        if k == "revolution":
            return _revolution_eval(self, X)[0]
        raise ValueError(k)

    def normal(self, X):
        k, p = self.kind, self.p
        if k == "plane":
            return np.tile(p["n"], (len(X), 1))
        if k == "sphere":
            w = X - p["c"]
            return w / np.linalg.norm(w, axis=1, keepdims=True)
        if k in ("cylinder", "cone", "torus"):
            w = X - (p["c"] if k != "cone" else p["v"])
            h = w @ p["a"]
            rad = w - np.outer(h, p["a"])
            rho = np.linalg.norm(rad, axis=1, keepdims=True)
            e = rad / np.maximum(rho, 1e-12)
            if k == "cylinder":
                return e
            if k == "cone":
                t = p["t"]
                return e * np.cos(t) - np.outer(np.ones(len(X)), p["a"]) * np.sin(t)
            q = (rho - p["R"]) * e + h[:, None] * p["a"]
            return q / np.maximum(np.linalg.norm(q, axis=1, keepdims=True), 1e-12)
        if k == "extrusion":
            return _extrusion_normal(self, X)
        if k == "revolution":
            return _revolution_eval(self, X)[1]
        raise ValueError(k)

    def transformed(self, T):
        L, t = T[:3, :3], T[:3, 3]
        q = Prim(self.kind)
        for key, val in self.p.items():
            if key in ("c", "v", "o"):
                q.p[key] = L @ val + t
            elif key in ("n", "a", "u", "w"):
                q.p[key] = L @ val
            elif key == "curve":
                q.p[key] = val.copy()
            else:
                q.p[key] = val
        if self.kind == "plane":
            q.p["d"] = float(q.p["n"] @ (L @ (self.p["n"] * self.p["d"]) + t))
        if self.kind == "extrusion":                 # curve coordinates stay valid in the mapped (u, w) basis
            shift = float(q.p["a"] @ t)
            q.p["h0"], q.p["h1"] = self.p["h0"] + shift, self.p["h1"] + shift
        q.faces = self.faces
        q.pts = self.pts @ L.T + t if len(self.pts) else self.pts
        q.rms, q.area = self.rms, self.area
        return q

    def copy(self):
        return self.transformed(np.eye(4))


# ---- minimal-sample constructors (RANSAC) -----------------------------------
def from_samples(kind, P, N, diag):
    """Build a candidate from oriented points P (k x 3), N (k x 3). Returns None if degenerate."""
    if kind == "plane":
        n = unit(N[0])
        return Prim("plane", n=n, d=float(P[0] @ n))
    if kind == "cylinder":
        a = np.cross(N[0], N[1])
        if np.linalg.norm(a) < 0.05:
            return None
        a = unit(a)
        d = P[1] - P[0]
        d = d - (d @ a) * a
        A = np.stack([N[0] - (N[0] @ a) * a, -(N[1] - (N[1] @ a) * a)], 1)
        ts, *_ = np.linalg.lstsq(A, d, rcond=None)
        c = P[0] + ts[0] * (N[0] - (N[0] @ a) * a)
        r = 0.5 * (abs(ts[0]) + abs(ts[1]))
        if not (1e-4 * diag < r < 2 * diag):
            return None
        return Prim("cylinder", c=c, a=a, r=r)
    if kind == "sphere":
        w0 = P[1] - P[0]
        b = N[0] @ N[1]
        den = 1 - b * b
        if den < 1e-3:
            return None
        t = (w0 @ N[0] - b * (w0 @ N[1])) / den
        s = (b * (w0 @ N[0]) - w0 @ N[1]) / den
        c = 0.5 * (P[0] + t * N[0] + P[1] + s * N[1])
        r = 0.5 * (np.linalg.norm(P[0] - c) + np.linalg.norm(P[1] - c))
        if not (1e-4 * diag < r < diag):
            return None
        return Prim("sphere", c=c, r=r)
    if kind == "cone":
        A = N[:3]
        if abs(np.linalg.det(A)) < 1e-3:
            return None
        v = np.linalg.solve(A, np.einsum("ij,ij->i", N[:3], P[:3]))
        if np.linalg.norm(v - P[0]) > 5 * diag:
            return None
        D = np.array([unit(P[i] - v) for i in range(3)])
        a = np.cross(D[1] - D[0], D[2] - D[0])
        if np.linalg.norm(a) < 1e-6:
            return None
        a = unit(a)
        if (D @ a).mean() < 0:
            a = -a
        t = float(np.mean(np.arccos(np.clip(D @ a, -1, 1))))
        if not (np.radians(3) < t < np.radians(87)):
            return None
        return Prim("cone", v=v, a=a, t=t)
    raise ValueError(kind)


N_SAMPLES = {"plane": 1, "cylinder": 2, "sphere": 2, "cone": 3}


# ---- least-squares refits ---------------------------------------------------
def fit(kind, X, init=None, fixed_axis=False):
    """Least-squares refit of a primitive to points X, starting from `init`."""
    if kind == "plane":
        if init is not None and fixed_axis:
            n = init.p["n"]
        else:
            c = X.mean(0)
            _, _, vt = np.linalg.svd(X - c, full_matrices=False)
            n = vt[-1]
            if init is not None and n @ init.p["n"] < 0:
                n = -n
        pr = Prim("plane", n=n, d=float(np.mean(X @ n)))
        pr.rms = float(np.sqrt(np.mean(pr.distance(X) ** 2)))
        return pr
    if kind == "torus" and init is None:
        init = torus_init(X)
        if init is None:
            return None
    if init is None:
        raise ValueError("init required")
    a0 = init.p.get("a")
    u, v = perp_basis(a0) if a0 is not None else (None, None)

    def axis(s):
        return a0 if fixed_axis else unit(a0 + s[0] * u + s[1] * v)

    if kind == "cylinder":
        x0 = np.r_[init.p["c"], 0, 0, init.p["r"]]
        f = lambda x: Prim("cylinder", c=x[:3], a=axis(x[3:5]), r=x[5]).distance(X)
    elif kind == "cone":
        x0 = np.r_[init.p["v"], 0, 0, init.p["t"]]
        f = lambda x: Prim("cone", v=x[:3], a=axis(x[3:5]), t=x[5]).distance(X)
    elif kind == "sphere":
        x0 = np.r_[init.p["c"], init.p["r"]]
        f = lambda x: Prim("sphere", c=x[:3], r=x[3]).distance(X)
    elif kind == "torus":
        x0 = np.r_[init.p["c"], 0, 0, init.p["R"], init.p["r"]]
        f = lambda x: Prim("torus", c=x[:3], a=axis(x[3:5]), R=x[5], r=x[6]).distance(X)
    else:
        raise ValueError(kind)
    try:
        res = least_squares(f, x0, method="trf", max_nfev=200, x_scale="jac")
        x = res.x
    except Exception:
        return None
    if kind == "cylinder":
        a = axis(x[3:5])
        c = x[:3] - ((x[:3] - X.mean(0)) @ a) * a     # put c at the region's axial middle
        pr = Prim("cylinder", c=c, a=a, r=abs(x[5]))
    elif kind == "cone":
        pr = Prim("cone", v=x[:3], a=axis(x[3:5]), t=x[5])
    elif kind == "sphere":
        pr = Prim("sphere", c=x[:3], r=abs(x[3]))
    else:
        pr = Prim("torus", c=x[:3], a=axis(x[3:5]), R=x[5], r=abs(x[6]))
    pr.rms = float(np.sqrt(np.mean(pr.distance(X) ** 2)))
    return pr


def torus_init(X, N=None):
    """Torus initial guess: axis as the line met by all normal lines (Pluecker), profile circle fit."""
    if N is None:
        return None
    M = np.cross(X, N)
    A = np.hstack([M, N])                 # unknown (a, m): n_i . m + m_i . a = 0
    _, _, vt = np.linalg.svd(A, full_matrices=False)
    a, m = vt[-1][:3], vt[-1][3:]
    na = np.linalg.norm(a)
    if na < 1e-9:
        return None
    a, m = a / na, m / na
    p0 = np.cross(a, m)                   # closest point of the axis to the origin
    w = X - p0
    h = w @ a
    rho = np.linalg.norm(w - np.outer(h, a), axis=1)
    hc, rc, rr = circle_fit_2d(np.c_[h, rho])
    if rr <= 0 or rc <= 0:
        return None
    return Prim("torus", c=p0 + hc * a, a=a, R=rc, r=rr)


def circle_fit_2d(Q):
    """Algebraic (Kasa) circle fit; returns (cx, cy, r)."""
    A = np.c_[2 * Q, np.ones(len(Q))]
    b = (Q ** 2).sum(1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy = sol[:2]
    r2 = sol[2] + cx * cx + cy * cy
    return cx, cy, float(np.sqrt(max(r2, 0)))


# ---- extrusion (spline swept along a direction) -----------------------------
def _extrusion_frame(pr):
    return pr.p["o"], pr.p["u"], pr.p["w"], pr.p["a"]


def _curve_dense(pr):
    from scipy.interpolate import splev
    tck = pr.p["tck"]
    s = np.linspace(0, 1, 800)
    return np.array(splev(s, tck)).T, np.array(splev(s, tck, der=1)).T


def _extrusion_distance(pr, X):
    o, u, w, a = _extrusion_frame(pr)
    Q = np.c_[(X - o) @ u, (X - o) @ w]
    C, D = _curve_dense(pr)
    from scipy.spatial import cKDTree
    d, i = cKDTree(C).query(Q)
    tang = D[i] / np.linalg.norm(D[i], axis=1, keepdims=True)
    nrm = np.c_[tang[:, 1], -tang[:, 0]]
    return np.einsum("ij,ij->i", Q - C[i], nrm)


def _extrusion_normal(pr, X):
    o, u, w, a = _extrusion_frame(pr)
    Q = np.c_[(X - o) @ u, (X - o) @ w]
    C, D = _curve_dense(pr)
    from scipy.spatial import cKDTree
    _, i = cKDTree(C).query(Q)
    tang = D[i] / np.linalg.norm(D[i], axis=1, keepdims=True)
    return np.outer(tang[:, 1], u) - np.outer(tang[:, 0], w)


def fit_extrusion(mesh, faces, tol):
    """Fit a spline extrusion to a mesh region whose normals are all perpendicular to one direction."""
    from scipy.interpolate import splprep
    N = mesh.face_normals[faces]
    wts = mesh.area_faces[faces]
    C = (N * wts[:, None]).T @ N
    evals, evecs = np.linalg.eigh(C)
    a = evecs[:, 0]
    if np.max(np.abs(N @ a)) > np.sin(np.radians(8)):
        return None
    sub = mesh.submesh([faces], append=True)
    X = sub.vertices
    h = X @ a
    mid = 0.5 * (h.min() + h.max())
    P = section_polyline(sub, a * mid, a)
    if P is None:
        return None
    u, w = perp_basis(a)
    o = a * mid
    Q = np.c_[(P - o) @ u, (P - o) @ w]
    keep = np.r_[True, np.linalg.norm(np.diff(Q, axis=0), axis=1) > 1e-9]
    Q = Q[keep]
    if len(Q) < 6:
        return None
    closed = np.linalg.norm(Q[0] - Q[-1]) < 2 * tol and len(Q) > 8
    if closed:
        Q = np.vstack([Q[:-1], Q[:1]]) if np.linalg.norm(Q[0] - Q[-1]) > 1e-9 else Q
    try:
        tck, _ = splprep(Q.T, s=len(Q) * (0.25 * tol) ** 2, k=3, per=1 if closed else 0)
    except Exception:
        return None
    pr = Prim("extrusion", o=o, u=u, w=w, a=a, h0=float(h.min()), h1=float(h.max()), closed=closed)
    pr.p["tck"] = tck
    pr.p["curve"] = Q
    d = pr.distance(X)
    pr.rms = float(np.sqrt(np.mean(d ** 2)))
    pr.pts = X
    return pr


def section_polyline(mesh, origin, normal):
    """Longest chained polyline of a plane section of a (possibly open) mesh."""
    import networkx as nx
    import trimesh
    segs = trimesh.intersections.mesh_plane(mesh, normal, origin)
    if len(segs) == 0:
        return None
    pts = segs.reshape(-1, 3)
    key = np.round(pts / 1e-7).astype(np.int64)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    G = nx.Graph()
    for i in range(len(segs)):
        a, b = inv[2 * i], inv[2 * i + 1]
        if a != b:
            G.add_edge(a, b)
    coords = {inv[j]: pts[j] for j in range(len(pts))}
    best = None
    for comp in nx.connected_components(G):
        H = G.subgraph(comp)
        ends = [v for v in H if H.degree(v) == 1]
        start = ends[0] if ends else next(iter(H))
        order = list(nx.dfs_preorder_nodes(H, start))
        if best is None or len(order) > len(best):
            best = order
    return np.array([coords[v] for v in best])


def _revolution_eval(pr, X):
    """Signed distance and unit normal for a spline profile (h, rho) revolved about (c, a)."""
    from scipy.interpolate import splev
    from scipy.spatial import cKDTree
    c, a = pr.p["c"], pr.p["a"]
    w = X - c
    h = w @ a
    radial = w - np.outer(h, a)
    rho = np.linalg.norm(radial, axis=1)
    e = radial / np.maximum(rho, 1e-12)[:, None]
    s = np.linspace(0, 1, 1500)
    C = np.array(splev(s, pr.p["tck"])).T
    D = np.array(splev(s, pr.p["tck"], der=1)).T
    d, i = cKDTree(C).query(np.c_[h, rho])
    t = D[i] / np.linalg.norm(D[i], axis=1, keepdims=True)
    n2 = np.c_[t[:, 1], -t[:, 0]]                     # normal in the (h, rho) half-plane
    sd = np.einsum("ij,ij->i", np.c_[h, rho] - C[i], n2)
    N = np.outer(n2[:, 0], a) + n2[:, 1][:, None] * e
    return sd, N


def _revolution_about(sub, c0, a, tol):
    """Spline profile of a mesh region about the axis (c0, a); returns a Prim or None."""
    from scipy.interpolate import splprep
    X = sub.vertices
    cent = X.mean(0) - c0
    e = cent - (cent @ a) * a
    if np.linalg.norm(e) < 1e-9:
        e, _ = perp_basis(a)
    e = unit(e)
    P = None
    for ang in (0.0, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0, 3.0):
        ee = np.cos(ang) * e + np.sin(ang) * np.cross(a, e)
        poly = section_polyline(sub, c0, np.cross(a, ee))
        if poly is None:
            continue
        poly = poly[((poly - c0) @ ee) > 0]
        if len(poly) >= 6:
            P = poly
            break
    if P is None:
        return None
    h = (P - c0) @ a
    rho = np.linalg.norm((P - c0) - np.outer(h, a), axis=1)
    Q = np.c_[h, rho]
    keep = np.r_[True, np.linalg.norm(np.diff(Q, axis=0), axis=1) > 1e-9]
    Q = Q[keep]
    if len(Q) < 6:
        return None
    try:
        tck, _ = splprep(Q.T, s=len(Q) * (0.25 * tol) ** 2, k=3)
    except Exception:
        return None
    pr = Prim("revolution", c=c0, a=a)
    pr.p["tck"] = tck
    pr.p["curve"] = Q
    d = pr.distance(X)
    pr.rms = float(np.sqrt(np.mean(d ** 2)))
    pr.maxdev = float(np.abs(d).max())
    pr.pts = X
    return pr


def fit_revolution(mesh, faces, tol, axes=()):
    """Fit a surface of revolution with a spline profile to a mesh region.
    Candidate axes: the line met by the region's normal lines (Pluecker fit), that line snapped to
    the world axes, and any known exact axes passed in `axes` (list of (point, direction))."""
    sub = mesh.submesh([faces], append=True)
    cands = list(axes)
    init = torus_init(sub.triangles_center, sub.face_normals)
    if init is not None:
        a, c0 = init.p["a"], init.p["c"]
        cands.append((c0, a))
        k = int(np.argmax(np.abs(a)))
        if np.degrees(np.arccos(min(1.0, abs(a[k])))) < 2.0:
            a2 = np.eye(3)[k] * np.sign(a[k])
            cands.append((c0, a2))
    best = None
    size = float(np.linalg.norm(sub.extents))
    fc, fn = sub.triangles_center, sub.face_normals
    for c0, a in cands:
        a = unit(np.asarray(a, float))
        c0 = np.asarray(c0, float)
        cr = np.cross(fn, a)
        nrm = np.linalg.norm(cr, axis=1)
        ok = nrm > 1e-6
        if ok.sum() > 3:
            dist = np.abs(np.einsum("ij,ij->i", fc[ok] - c0, cr[ok]) / nrm[ok])
            if np.median(dist) > 0.05 * size + 3 * tol:
                continue                     # normal lines do not meet this axis
        pr = _revolution_about(sub, c0, a, tol)
        if pr is not None and (best is None or pr.maxdev < best.maxdev):
            best = pr
    return best


def match_score(P, Q, tol):
    """Median distance of P's region points to Q's surface, or inf if kinds differ."""
    if P.kind != Q.kind or not len(P.pts):
        return np.inf
    return float(np.median(np.abs(Q.distance(P.pts))))
