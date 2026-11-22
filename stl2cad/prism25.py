"""2.5D reconstruction: stacked prisms.

A part whose faces are all either perpendicular or parallel to one world axis is a stack of
prisms. The cap faces give the levels; between two levels the cross-section is constant. Each
section is fitted with lines and circular arcs, extruded over its height range, and the extrusions
are fused. Slabs with the same section are merged into one extrusion.

The result is accepted only if its volume and its distance to the mesh vertices are within
tolerance; otherwise the caller falls back to the surface arrangement method.
"""
import numpy as np

from .core import Result

AX = "xyz"


# ----------------------------------------------------------------------------- detection
def prism_axis(mesh, min_fraction=0.999, deg=1.0):
    """World axis k along which the mesh is 2.5D, or None. A face qualifies if its normal is
    parallel (cap) or perpendicular (wall) to the axis within `deg`."""
    n, a = mesh.face_normals, mesh.area_faces
    best, kbest = 0.0, None
    for k in range(3):
        c = np.abs(n[:, k])
        ok = (c > np.cos(np.radians(deg))) | (c < np.sin(np.radians(deg)))
        f = a[ok].sum() / a.sum()
        if f > best:
            best, kbest = f, k
    return (kbest, best) if best >= min_fraction else (None, best)


def levels(mesh, k, tol):
    """Heights of the cap faces along axis k, clustered within tol."""
    cap = np.abs(mesh.face_normals[:, k]) > np.cos(np.radians(1.0))
    h = np.sort(mesh.triangles_center[cap, k])
    if not len(h):
        return []
    out, cur = [], [h[0]]
    for x in h[1:]:
        if x - cur[-1] > tol:
            out.append(float(np.mean(cur)))
            cur = [x]
        else:
            cur.append(x)
    out.append(float(np.mean(cur)))
    return out


# ----------------------------------------------------------------------------- sections
def section_loops(mesh, k, z, tol):
    """Closed section loops at height z as ordered 2D points that are mesh vertices (exact
    surface points), not plane-edge crossings (chord points)."""
    from trimesh.intersections import mesh_plane
    i, j = [m for m in range(3) if m != k]
    nrm = np.zeros(3)
    nrm[k] = 1.0
    org = np.zeros(3)
    org[k] = z
    segs, fidx = mesh_plane(mesh, nrm, org, return_faces=True)
    if not len(segs):
        return []
    S = segs[:, :, [i, j]]
    scale = 1e-6 * float(np.linalg.norm(mesh.extents))
    key = lambda p: (int(round(p[0] / scale)), int(round(p[1] / scale)))
    ends = {}
    for s in range(len(S)):
        for e in range(2):
            ends.setdefault(key(S[s, e]), []).append((s, e))
    used = np.zeros(len(S), bool)
    loops = []
    for s0 in range(len(S)):
        if used[s0]:
            continue
        order, s, e = [], s0, 0
        while not used[s]:
            used[s] = True
            order.append((s, e))
            nxt = [(t, f) for t, f in ends[key(S[s, 1 - e])] if t != s and not used[t]]
            if not nxt:
                break
            s, e = nxt[0]
        if len(order) >= 3:
            loops.append(order)
    V = mesh.vertices[:, [i, j]]
    F = mesh.faces
    out = []
    for order in loops:
        poly = np.array([S[s, e] for s, e in order])
        pts, par = [], []
        for n_, (s, e) in enumerate(order):
            a, b = S[s, e], S[s, 1 - e]
            d = b - a
            L2 = float(d @ d) + 1e-30
            for v in V[F[fidx[s]]]:
                t = np.clip((v - a) @ d / L2, 0.0, 1.0)
                if np.linalg.norm(a + t * d - v) < 0.05 * np.sqrt(L2) + 4 * tol:
                    pts.append(v)
                    par.append(n_ + t)
        pts, par = np.array(pts), np.array(par)
        o = np.argsort(par, kind="stable")
        pts = pts[o]
        keep = [0]
        for m in range(1, len(pts)):
            if np.linalg.norm(pts[m] - pts[keep[-1]]) > scale * 10:
                keep.append(m)
        pts = pts[keep]
        if len(pts) > 2 and np.linalg.norm(pts[-1] - pts[0]) < scale * 10:
            pts = pts[:-1]
        # a sliver loop (mean width below tol) is a gap between two faceted shells that meet
        # along a curve, not a feature
        per = float(np.linalg.norm(np.diff(np.vstack([poly, poly[:1]]), axis=0), axis=1).sum())
        if len(pts) >= 3 and 2 * abs(_area(poly)) > tol * per:
            out.append((pts, poly))
    return out


def _area(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


# ----------------------------------------------------------------------------- 2D fitting
def _line_fit(Q):
    c = Q.mean(0)
    _, _, vt = np.linalg.svd(Q - c)
    d = vt[0]
    nrm = np.array([-d[1], d[0]])
    return c, d, float(np.abs((Q - c) @ nrm).max())


def _circle_fit(Q):
    A = np.c_[2 * Q, np.ones(len(Q))]
    x, *_ = np.linalg.lstsq(A, (Q ** 2).sum(1), rcond=None)
    c = x[:2]
    r = float(np.sqrt(max(x[2] + c @ c, 0.0)))
    return c, r, float(np.abs(np.linalg.norm(Q - c, axis=1) - r).max())


def _sweep(Q, c, max_step=40.0):
    """Signed angle swept along the ordered points about c (radians)."""
    ang = np.unwrap(np.arctan2(Q[:, 1] - c[1], Q[:, 0] - c[0]))
    d = np.diff(ang)
    if len(d) and not (np.all(d > -1e-9) or np.all(d < 1e-9)):
        return None                                   # not monotonic: not one arc
    if len(d) and np.abs(d).max() > np.radians(max_step):
        return None                                   # a tessellated arc has small facet steps
    if len(d) >= 2 and np.abs(d).max() > 3.0 * np.abs(d).min() + 1e-9:
        return None                                   # ... of similar size: a far point is not on it

    return float(ang[-1] - ang[0])


def _full_circle(P, tol):
    """Circle through a closed loop, or None. Points more than tol inside the circle are chord
    points (an STL merged from overlapping shells has vertices on the facets of the other shell);
    they are dropped and the circle refitted on the points that lie on it."""
    if len(P) < 6:
        return None
    keep = np.ones(len(P), bool)
    for _ in range(4):
        c, r, _ = _circle_fit(P[keep])
        d = np.linalg.norm(P - c, axis=1)
        new = d > r - tol
        if new.sum() < 6 or (new == keep).all():
            keep = new if new.sum() >= 6 else keep
            break
        keep = new
    c, r, _ = _circle_fit(P[keep])
    d = np.linalg.norm(P - c, axis=1)
    ang = np.sort(np.arctan2(P[keep, 1] - c[1], P[keep, 0] - c[0]))
    gap = np.diff(np.r_[ang, ang[0] + 2 * np.pi]).max()
    sag = r * (1 - np.cos(gap / 2))
    if (d - r).max() > tol or (r - d).max() > tol + sag or gap > np.radians(30):
        return None
    Q = np.vstack([P, P[:1]])
    a = np.unwrap(np.arctan2(Q[:, 1] - c[1], Q[:, 0] - c[0]))
    if abs(abs(a[-1] - a[0]) - 2 * np.pi) > 0.2:     # the loop goes once around the centre
        return None
    return dict(c=c, r=r)


def fit_chain(Q, tol):
    """Greedy lines and arcs along the ordered open point chain Q; pieces share end points."""
    M = len(Q) - 1
    span = float(np.ptp(Q, axis=0).max())
    pieces, i = [], 0
    while i < M:
        jl = i + 1                                     # longest line from i
        while jl + 1 <= M and _line_fit(Q[i:jl + 2])[2] < tol:
            jl += 1
        ja, arc = i + 1, None                          # longest arc from i (at least 3 points)
        j = i + 2
        while j <= M:
            cc, rr, ee = _circle_fit(Q[i:j + 1])
            if ee >= tol or rr > 20 * span or rr < 2 * tol:
                if arc is not None or j > i + 2:
                    break
                j += 1
                continue
            sw = _sweep(Q[i:j + 1], cc)
            if sw is None or abs(sw) > 1.9 * np.pi:
                break
            ja, arc = j, dict(c=cc, r=rr)
            j += 1
        if arc is not None and ja > jl:
            pieces.append(["arc", i, ja, arc])
            i = ja
        else:
            pieces.append(["line", i, jl, None])
            i = jl
    merged = [pieces[0]]                               # merge neighbours on one line or one circle
    for p in pieces[1:]:
        q = merged[-1]
        if p[0] == q[0]:
            Qm = Q[q[1]:p[2] + 1]
            if p[0] == "line" and _line_fit(Qm)[2] < tol:
                q[2] = p[2]
                continue
            if p[0] == "arc":
                cc, rr, ee = _circle_fit(Qm)
                sw = _sweep(Qm, cc)
                if ee < tol and sw is not None and abs(sw) < 1.9 * np.pi:
                    q[2], q[3] = p[2], dict(c=cc, r=rr)
                    continue
        merged.append(p)
    for p in merged:                                  # final parameters on the final ranges
        Qm = Q[p[1]:p[2] + 1]
        if p[0] == "line":
            c0, d0, _ = _line_fit(Qm)
            p[3] = dict(c=c0, d=d0)
        else:
            cc, rr, _ = _circle_fit(Qm)
            p[3] = dict(c=cc, r=rr)
    return [tuple(p) for p in merged]


def fit_loop(P, tol):
    """Lines and arcs through the ordered closed point list P. Returns a list of pieces
    (kind, i0, i1, params) over a rotated copy of P, and the rotated P."""
    M = len(P)
    full = _full_circle(P, tol)
    if full is not None:
        return [("circle", 0, M, full)], P
    # start at the sharpest corner, so that no piece wraps around the start
    a = P - np.roll(P, 1, axis=0)
    b = np.roll(P, -1, axis=0) - P
    cosang = np.einsum("ij,ij->i", a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-30)
    s = int(np.argmin(cosang))
    P = np.roll(P, -s, axis=0)
    return fit_chain(np.vstack([P, P[:1]]), tol), P


def _project(p, piece):
    kind, _, _, g = piece
    if kind == "line":
        return g["c"] + ((p - g["c"]) @ g["d"]) * g["d"]
    v = p - g["c"]
    return g["c"] + g["r"] * v / (np.linalg.norm(v) + 1e-30)


def _intersect(pa, pb, near):
    """Intersection of two fitted curves closest to `near`, or None."""
    ka, kb = pa[0], pb[0]
    ga, gb = pa[3], pb[3]
    if ka == "line" and kb == "line":
        A = np.c_[ga["d"], -gb["d"]]
        if abs(np.linalg.det(A)) < 1e-9:
            return None
        t = np.linalg.solve(A, gb["c"] - ga["c"])
        return ga["c"] + t[0] * ga["d"]
    if ka == "arc" and kb == "arc":
        c0, r0, c1, r1 = ga["c"], ga["r"], gb["c"], gb["r"]
        d = np.linalg.norm(c1 - c0)
        if d < 1e-12 or d > r0 + r1 or d < abs(r0 - r1):
            return None
        a = (r0 ** 2 - r1 ** 2 + d ** 2) / (2 * d)
        h = np.sqrt(max(r0 ** 2 - a ** 2, 0))
        u = (c1 - c0) / d
        m = c0 + a * u
        w = np.array([-u[1], u[0]])
        cand = [m + h * w, m - h * w]
    else:
        ln, ar = (ga, gb) if ka == "line" else (gb, ga)
        f = ln["c"] - ar["c"]
        b = f @ ln["d"]
        disc = b * b - (f @ f - ar["r"] ** 2)
        if disc < 0:
            return None
        cand = [ln["c"] + (-b + s * np.sqrt(disc)) * ln["d"] for s in (1, -1)]
    return min(cand, key=lambda q: np.linalg.norm(q - near))


def loop_vertices(pieces, P, tol):
    """Corner points between consecutive pieces: the curves' intersection where they cross at an
    angle, otherwise the mean of the data point's projections (tangent junctions)."""
    n = len(pieces)
    Q = np.vstack([P, P[:1]])
    verts = []
    for m in range(n):
        pa, pb = pieces[m - 1], pieces[m]
        near = Q[pb[1]]
        x = _intersect(pa, pb, near)
        if x is None or np.linalg.norm(x - near) > 3 * tol:
            x = 0.5 * (_project(near, pa) + _project(near, pb))
        verts.append(x)
    return verts                                       # verts[m] = start of piece m


# ----------------------------------------------------------------------------- profiles
def _rot(P, c, ang):
    R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    return (P - c) @ R.T + c


def rotational_period(P, tol, n_max=360):
    """(centre, N) when the closed loop P is unchanged by a rotation of 2 pi / N about its
    centre (N >= 3), else None. N is the strongest harmonic of the radius as a function of angle,
    confirmed by rotating every point onto the loop."""
    import shapely
    if len(P) < 30:
        return None
    x, y = P[:, 0], P[:, 1]                            # area centroid (point density is uneven)
    cr = x * np.roll(y, -1) - np.roll(x, -1) * y
    A = cr.sum() / 2
    c = np.array([((x + np.roll(x, -1)) * cr).sum(), ((y + np.roll(y, -1)) * cr).sum()]) / (6 * A)
    th = np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0])
    a = np.unwrap(np.r_[th, th[0]])
    if abs(abs(a[-1] - a[0]) - 2 * np.pi) > 0.2:
        return None                                    # the loop does not go once around c
    o = np.argsort(th)
    t, rho = th[o], np.linalg.norm(P[o] - c, axis=1)
    g = np.linspace(-np.pi, np.pi, 4096, endpoint=False)
    rg = np.interp(g, t, rho, period=2 * np.pi)
    spec = np.abs(np.fft.rfft(rg - rg.mean()))
    ring = shapely.LinearRing(P)
    for N in np.argsort(spec[3:n_max + 1])[::-1][:3] + 3:
        if spec[N] < 0.05 * spec[1:].max():
            break
        Q = _rot(P, c, 2 * np.pi / N)
        if np.max(shapely.distance(shapely.points(Q), ring)) < 2 * tol:
            return c, int(N)
    return None


def pattern_loop(P, tol):
    """One period of a rotationally repeated loop, fitted on the points of all periods pooled
    (a coarse mesh has few points per tooth; together they sample one tooth densely). Returns
    (centre, N, edges of one period) or None."""
    per = rotational_period(P, tol)
    if per is None:
        return None
    c, N = per
    T = 2 * np.pi / N
    th = np.mod(np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]), 2 * np.pi)
    k = np.floor(th / T)
    pooled = _rot(P, c, -k * 0)[:0]
    pooled = np.vstack([_rot(P[k == m], c, -m * T) for m in range(N) if np.any(k == m)])
    tp = np.mod(np.arctan2(pooled[:, 1] - c[1], pooled[:, 0] - c[0]), 2 * np.pi)
    o = np.argsort(tp)
    pooled, tp = pooled[o], tp[o]
    keep = np.r_[True, np.linalg.norm(np.diff(pooled, axis=0), axis=1) > 1e-3 * tol]
    pooled, tp = pooled[keep], tp[keep]
    # one and a half periods, so that the pieces at the period boundary are fitted with support
    # on both sides; the period is then cut at junctions one period apart
    ext = np.vstack([pooled, _rot(pooled[tp < 0.5 * T], c, T)])
    pieces = fit_chain(ext, tol)
    if len(pieces) < 3:
        return None
    verts = [None] + [_junction(pieces[m - 1], pieces[m], ext[pieces[m][1]], tol) for m in range(1, len(pieces))]
    angs = [None] + [np.mod(np.arctan2(v[1] - c[1], v[0] - c[0]), 2 * np.pi) for v in verts[1:]]
    best = None
    for a_ in range(1, len(pieces)):
        if not (0.1 * T < angs[a_] < 0.6 * T):
            continue
        for b_ in range(a_ + 1, len(pieces)):
            if np.linalg.norm(verts[b_] - _rot(verts[a_][None], c, T)[0]) < 2 * tol:
                if best is None or b_ - a_ < best[1] - best[0]:
                    best = (a_, b_)
                break
    if best is None:
        return None
    a_, b_ = best
    V = [verts[m] for m in range(a_, b_)] + [_rot(verts[a_][None], c, T)[0]]
    edges = []
    for m in range(a_, b_):
        p = pieces[m]
        A, B = V[m - a_], V[m - a_ + 1]
        if p[0] == "line":
            edges.append(("line", A, B))
        else:
            mid = _project(ext[(p[1] + p[2]) // 2], p)
            edges.append(("arc", A, mid, B))
    return c, N, edges


def _junction(pa, pb, near, tol):
    x = _intersect(pa, pb, near)
    if x is None or np.linalg.norm(x - near) > 3 * tol:
        x = 0.5 * (_project(near, pa) + _project(near, pb))
    return x


MIN_PATTERN = 6                       # rotational patterns of fewer repeats are fitted as plain loops


class Loop:
    """A fitted closed loop: list of edges ("line", a, b) / ("arc", a, mid, b) / ("circle", c, r).
    pattern: (centre, N, edges of one period) when the loop repeats N times about a centre."""

    def __init__(self, edges, area, pattern=None):
        self.edges, self.area, self.pattern = edges, area, pattern


def build_loop(P, tol):
    pieces, P = fit_loop(P, tol)
    if pieces[0][0] == "circle":
        g = pieces[0][3]
        return Loop([("circle", g["c"], g["r"])], _area(P))
    pat = pattern_loop(P, tol) if len(pieces) >= 12 else None
    if pat is not None and pat[1] < MIN_PATTERN:
        pat = None                                       # a few repeats are clearer written out
    if pat is not None:
        c, N, one = pat
        edges = []
        for m in range(N):
            a = 2 * np.pi * m / N
            for e in one:
                edges.append((e[0],) + tuple(_rot(np.asarray(q)[None], c, a)[0] for q in e[1:]))
        for n_ in range(1, len(edges)):                # identical shared vertices around the loop
            edges[n_] = (edges[n_][0], edges[n_ - 1][-1]) + edges[n_][2:]
        edges[-1] = edges[-1][:-1] + (edges[0][1],)
        return Loop(edges, _area(P), pattern=(c, N, one))
    verts = loop_vertices(pieces, P, tol)
    Q = np.vstack([P, P[:1]])
    edges = []
    for m, p in enumerate(pieces):
        a, b = verts[m], verts[(m + 1) % len(pieces)]
        if np.linalg.norm(b - a) < 1e-9:
            continue
        if p[0] == "line":
            edges.append(("line", a, b))
        else:
            mid = _project(Q[(p[1] + p[2]) // 2], p)
            if (p[2] - p[1]) % 2:
                mid = _project(0.5 * (Q[(p[1] + p[2]) // 2] + Q[(p[1] + p[2]) // 2 + 1]), p)
            edges.append(("arc", a, mid, b))
    return Loop(edges, _area(P))


def consolidate(loops, tol):
    """One value per design dimension: circle and pattern centres within tol become one centre,
    circle radii within tol one radius (a hole and a boss turned on one axis share it)."""
    cs, rs = [], []
    for L in loops:
        if L.edges and L.edges[0][0] == "circle" and len(L.edges) == 1:
            cs.append(L.edges[0][1])
            rs.append(L.edges[0][2])
        if L.pattern is not None:
            cs.append(L.pattern[0])

    def cluster(vals, eps):
        vals = [np.atleast_1d(np.asarray(v, float)) for v in vals]
        rep = []
        for v in vals:
            for g in rep:
                if np.linalg.norm(g[0] - v) < eps:
                    g[1].append(v)
                    break
            else:
                rep.append([v, [v]])
        return [(g[0], np.mean(g[1], axis=0)) for g in rep]
    cmap, rmap = cluster(cs, tol), cluster(rs, tol)

    def snap(v, table):
        v = np.atleast_1d(np.asarray(v, float))
        for key, mean in table:
            if np.linalg.norm(key - v) < tol:
                return mean
        return v
    out = []
    for L in loops:
        if L.edges and L.edges[0][0] == "circle" and len(L.edges) == 1:
            c = snap(L.edges[0][1], cmap)
            r = float(snap(L.edges[0][2], rmap)[0])
            out.append(Loop([("circle", c, r)], L.area))
        elif L.pattern is not None:
            c0, N, one = L.pattern
            c = snap(c0, cmap)
            d = c - c0                                  # move the whole outline with its centre
            one = [(e[0],) + tuple(np.asarray(q) + d for q in e[1:]) for e in one]
            edges = [(e[0],) + tuple(np.asarray(q) + d for q in e[1:]) for e in L.edges]
            out.append(Loop(edges, L.area, pattern=(c, N, one)))
        else:
            out.append(L)
    return out


def _same_section(A, B, tol):
    """Two raw sections (lists of polylines) describe the same region."""
    if len(A) != len(B):
        return False
    from scipy.spatial import cKDTree
    pa = np.vstack([p for _, p in A])
    pb = np.vstack([p for _, p in B])
    da, _ = cKDTree(pb).query(pa)
    db, _ = cKDTree(pa).query(pb)
    return max(da.max(), db.max()) < 2 * tol


# ----------------------------------------------------------------------------- solids
def _wire(loop, k, z):
    import cadquery as cq
    i, j = [m for m in range(3) if m != k]

    def V(p):
        q = [0.0, 0.0, 0.0]
        q[i], q[j], q[k] = float(p[0]), float(p[1]), float(z)
        return cq.Vector(*q)
    D = [0.0, 0.0, 0.0]
    D[k] = 1.0
    es = []
    for e in loop.edges:
        if e[0] == "circle":
            es.append(cq.Edge.makeCircle(e[2], V(e[1]), cq.Vector(*D)))
        elif e[0] == "line":
            es.append(cq.Edge.makeLine(V(e[1]), V(e[2])))
        else:
            es.append(cq.Edge.makeThreePointArc(V(e[1]), V(e[2]), V(e[3])))
    return cq.Wire.assembleEdges(es)


def _faces(loops, k, z):
    """Planar faces from loops: each outer loop (largest first) takes the loops it contains."""
    import cadquery as cq
    from matplotlib.path import Path
    order = sorted(range(len(loops)), key=lambda m: -abs(loops[m][1].area))
    pts = {m: loops[m][0] for m in order}
    depth = {}
    for m in order:
        depth[m] = sum(1 for o in order if o != m and abs(loops[o][1].area) > abs(loops[m][1].area)
                       and Path(pts[o]).contains_point(pts[m][0]))
    faces = []
    for m in order:
        if depth[m] % 2:
            continue
        inner = [o for o in order if depth[o] == depth[m] + 1 and Path(pts[m]).contains_point(pts[o][0])]
        outer_w = _wire(loops[m][1], k, z)
        inner_w = [_wire(loops[o][1], k, z) for o in inner]
        faces.append((cq.Face.makeFromWires(outer_w, inner_w), m, inner))
    return faces


def _distance_check(shape, mesh, n=4000, seed=0):
    """Max distances mesh -> model (a sample of mesh vertices, which lie on the true surfaces) and
    model -> mesh (a sample of points on the model surface; includes the chord sag of the mesh)."""
    import trimesh
    from .core import shape_to_mesh
    tm = shape_to_mesh(shape, tolerance=2e-5 * float(np.linalg.norm(mesh.extents)), angular=0.05)
    rng = np.random.default_rng(seed)
    V = mesh.vertices
    idx = rng.choice(len(V), min(n, len(V)), replace=False)
    _, d1, _ = trimesh.proximity.closest_point(tm, V[idx])
    S, _ = trimesh.sample.sample_surface(tm, n, seed=seed)
    _, d2, _ = trimesh.proximity.closest_point(mesh, S)
    return float(d1.max()), float(d2.max())


def reconstruct(mesh, tol, log=print):
    """Stacked-prism solid, or None when the mesh is not 2.5D or the fit fails the checks.
    A mesh of several closed bodies (an STL exported as overlapping shells) is rebuilt body by
    body, each checked against its own shell, and the solids are fused."""
    import time
    t0 = time.time()
    if mesh.body_count == 1:
        return _reconstruct_body(mesh, tol, log)
    bodies = mesh.split(only_watertight=False)
    log(f"  2.5D: {len(bodies)} closed shells; each is rebuilt and checked on its own")
    res = []
    for b in bodies:
        r = _reconstruct_body(b, tol, log)
        if r is None:
            return None
        res.append(r)
    solid = res[0].shape
    for r in res[1:]:
        solid = solid.fuse(r.shape)
    solid = solid.clean()
    vsum = sum(b.volume for b in bodies)
    log(f"  2.5D: fused {len(res)} shells into {len(solid.Solids())} solid(s); union volume "
        f"{solid.Volume():.1f} mm^3 vs {vsum:.1f} summed over the shells (overlap "
        f"{vsum - solid.Volume():.1f}) ({time.time() - t0:.1f}s)")
    if not solid.isValid():
        return None
    feats = []
    for r in res:
        feats += r.extras["prism25"]["features"]
    out = Result(shape=solid, mode="exact")
    out.extras = dict(prims=[], analytic_fraction=1.0, summary=f"2.5D: {len(feats)} extrusions",
                      volume_error=max(r.extras["volume_error"] for r in res), brep=None, valid=True,
                      passed=True, overlapping_shells=len(bodies),
                      prism25=dict(axis=res[0].extras["prism25"]["axis"], features=feats))
    return out


def _reconstruct_body(mesh, tol, log=print):
    try:
        r = _reconstruct_body_(mesh, tol, log)
        if r is None:
            r2 = _reconstruct_body_(mesh, tol, lambda *a: None, snap=False)
            if r2 is not None:
                log("  2.5D: passes without snapping; snapping skipped")
            return r2
        return r
    except Exception as e:                       # any construction failure: use the arrangement
        log(f"  2.5D: construction failed ({type(e).__name__}: {str(e)[:80]}); using the surface arrangement")
        return None


def _same_loop(A, B, tol):
    """Two section polylines trace the same closed curve: every vertex of each lies within 2 tol
    of the other polyline (nearest-vertex distance would mis-judge chords cut at other heights)."""
    import shapely
    if abs(abs(_area(A)) - abs(_area(B))) > 0.02 * abs(_area(A)) + tol ** 2:
        return False
    la = shapely.LinearRing(A)
    lb = shapely.LinearRing(B)
    da = np.max(shapely.distance(shapely.points(A), lb))
    db = np.max(shapely.distance(shapely.points(B), la))
    return max(da, db) < 2 * tol


def loop_tracks(slabs, secs, tol):
    """Follow each section loop through consecutive slabs. A track is a loop that stays the same
    over a height range; its nesting depth in each slab says whether it adds material (even) or
    removes it (odd). Returns tracks: dict(z0, z1, P, poly, op, order)."""
    from matplotlib.path import Path
    tracks, open_ = [], []
    for si, ((a, b), sec) in enumerate(zip(slabs, secs)):
        polys = [poly for _, poly in sec]
        areas = [abs(_area(p)) for p in polys]
        depth = []
        for m, p in enumerate(polys):
            depth.append(sum(1 for o, q in enumerate(polys) if o != m and areas[o] > areas[m]
                             and Path(q).contains_point(p[0])))
        nxt = []
        for m, (P, poly) in enumerate(sec):
            op = "add" if depth[m] % 2 == 0 else "cut"
            hit = None
            for t in open_:
                if abs(t["z1"] - a) < 1e-9 and t["op"] == op and _same_loop(t["poly"], poly, tol):
                    hit = t
                    break
            if hit is not None:
                open_ = [t for t in open_ if t is not hit]
                hit["z1"] = b
                hit["order"] = max(hit["order"], depth[m])
                nxt.append(hit)
            else:
                t = dict(z0=a, z1=b, P=P, poly=poly, op=op, order=depth[m])
                tracks.append(t)
                nxt.append(t)
        open_ = nxt
    return tracks


def _reconstruct_body_(mesh, tol, log=print, snap=True):
    import time
    import cadquery as cq
    t0 = time.time()
    k, frac = prism_axis(mesh)
    if k is None:
        log(f"  2.5D: no prism axis (best face fraction {frac:.3f})")
        return None
    lv = levels(mesh, k, tol)
    slabs = [(a, b) for a, b in zip(lv[:-1], lv[1:]) if b - a > 2 * tol]
    if not slabs:
        return None
    secs = [section_loops(mesh, k, 0.5 * (a + b), tol) for a, b in slabs]
    tracks = loop_tracks(slabs, secs, tol)
    # material is built in nesting order: outlines, then the holes in them, then islands in the
    # holes, and so on; at every height this reproduces the section
    tracks.sort(key=lambda t: (t["order"], t["op"] == "cut"))
    if not tracks or tracks[0]["op"] != "add":
        return None
    loops = consolidate([build_loop(t["P"], tol) for t in tracks], tol)
    if snap:
        from .snap2d import snap_loops
        loops = snap_loops(loops, [t["P"] for t in tracks], tol, log=log)
    feats, solid = [], None
    for t, loop in zip(tracks, loops):
        f = cq.Face.makeFromWires(_wire(loop, k, t["z0"]), [])
        D = [0.0, 0.0, 0.0]
        D[k] = t["z1"] - t["z0"]
        body = cq.Solid.extrudeLinear(f, cq.Vector(*D))
        if solid is None:
            solid = body
        elif t["op"] == "add":
            solid = solid.fuse(body)
        else:
            solid = solid.cut(body)
        feats.append(dict(z0=t["z0"], z1=t["z1"], op=t["op"], loop=loop))
    solid = solid.clean()
    n_edges = sum(len(f["loop"].edges) for f in feats)
    vol_err = abs(solid.Volume() / mesh.volume - 1)
    valid = solid.isValid() and len(solid.Solids()) == 1
    dmax, dback = _distance_check(solid, mesh) if valid else (np.inf, np.inf)
    n_add = sum(f["op"] == "add" for f in feats)
    log(f"  2.5D along {AX[k]}: {len(lv)} levels, {n_add} extrusion(s) and {len(feats) - n_add} cut(s), "
        f"{n_edges} profile edges, volume error {100 * vol_err:.3f} %, max distance mesh->model "
        f"{dmax:.4f} mm, model->mesh {dback:.4f} mm, valid={valid} ({time.time() - t0:.1f}s)")
    # the two-way distance is the real test: every mesh vertex within 1.5 tol of the model, the
    # model within 4 tol of the mesh (chord sag of the facets plus fitting error). The volume
    # differs by the sag of the facets (1.2 % on a thin gear ring), so it is only a gross check.
    if not (valid and vol_err < 0.03 and dmax < 1.5 * tol and dback < 4 * tol):
        log("  2.5D: checks failed; using the surface arrangement")
        return None
    res = Result(shape=solid, mode="exact")
    res.extras = dict(prims=[], analytic_fraction=1.0, summary=f"2.5D: {len(feats)} features",
                      volume_error=vol_err, brep=None, valid=True, passed=True,
                      prism25=dict(axis=k, features=feats))
    return res
