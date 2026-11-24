"""Design-intent snapping of fitted 2D profile loops (2.5D method).

The fitted lines and arcs are within tolerance of the mesh but carry fitting noise: a vertical
edge 0.03 mm off vertical, a length of 100.000001, an inner corner arc not quite tangent to its
lines. Snapping works on the geometry, not on vertices:

1. lines within tolerance of horizontal or vertical become exactly so;
2. the offsets of such lines and the coordinates of arc and circle centres are pooled per axis,
   clustered within tol/2 and rounded to the coarsest clean value within tol/4 (mm and inch
   grids); radii likewise;
3. an arc tangent to an axis-parallel line is moved to be exactly tangent;
4. corners are recomputed from the snapped geometry (intersections; tangent points);
5. a loop is kept only if every original section point lies within tol of the snapped loop.
Patterned loops (gear teeth) are left as fitted.
"""
import os

import numpy as np

DEBUG = os.environ.get("STL2CAD_SNAPDBG") == "1"

MM_GRIDS = [10.0, 5.0, 2.5, 1.0, 0.5, 0.25, 0.1, 0.05, 0.01]
IN_GRIDS = [25.4 * g for g in (1.0, 0.5, 0.25, 0.125, 0.1, 0.0625, 0.05, 0.03125, 0.025,
                               0.015625, 0.01, 0.005, 0.001)]
GRIDS = MM_GRIDS


def choose_grids(values, eps):
    """mm or inch grids for a part: the system in which more of its values lie within eps of a
    0.1 mm or 0.25 mm multiple, or of a 0.01 in multiple (ties go to mm)."""
    v = np.abs(np.asarray(values, float))
    v = v[v > eps]
    if not len(v):
        return MM_GRIDS
    on = lambda g: np.abs(v - np.round(v / g) * g) <= eps
    return IN_GRIDS if int(on(0.254).sum()) > int((on(0.1) | on(0.25)).sum()) else MM_GRIDS


def nice(v, eps, grids=None):
    """Coarsest clean value within eps of v (v itself if none)."""
    for g in (grids or GRIDS):
        r = round(v / g) * g
        if abs(r - v) <= eps:
            return float(r)
    return float(v)


def nice_any(v, eps):
    """Clean value within eps of v on the coarser of the mm and inch grids (53.975 = 2 1/8 in
    rather than 53.98 mm; 10 mm rather than 0.3937 in)."""
    cand = []
    for grids in (MM_GRIDS, IN_GRIDS):
        for g in grids:
            r = round(v / g) * g
            if abs(r - v) <= eps:
                cand.append((abs(r - v) > eps / 10, -g, float(r)))   # a near-exact hit first
                break
    return min(cand)[2] if cand else float(v)


def _circle3(p, q, r):
    ax, ay = p
    bx, by = q
    cx, cy = r
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-15:
        return None, None
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    c = np.array([ux, uy])
    return c, float(np.linalg.norm(np.asarray(p) - c))


class _Prim:
    def __init__(self, kind, **kw):
        self.kind = kind
        self.__dict__.update(kw)


def _prims(loop):
    out = []
    for e in loop.edges:
        if e[0] == "line":
            a, b = np.asarray(e[1], float), np.asarray(e[2], float)
            d = b - a
            L = float(np.linalg.norm(d))
            out.append(_Prim("line", a=a, b=b, d=d / L, L=L, axis=None, off=None))
        else:
            a, m, b = (np.asarray(q, float) for q in e[1:])
            c, r = _circle3(a, m, b)
            if c is None:
                return None
            out.append(_Prim("arc", a=a, m=m, b=b, c=c, r=r))
    return out


class _Pool:
    """Values of one kind (x, y or radius): clustered within eps, cluster mean rounded."""

    def __init__(self, eps, round_eps):
        self.vals, self.eps, self.round_eps = [], eps, round_eps

    def add(self, v):
        self.vals.append(float(v))

    def table(self, grids=None):
        out = []
        for v in sorted(self.vals):
            if out and v - out[-1][-1] <= self.eps:
                out[-1].append(v)
            else:
                out.append([v])
        return [(min(g), max(g), nice(float(np.mean(g)), self.round_eps, grids)) for g in out]

    @staticmethod
    def get(table, v):
        for lo, hi, rep in table:
            if lo - 1e-12 <= v <= hi + 1e-12:
                return rep
        return v


def _project_line(p, line):
    if line.axis == "v":
        return np.array([line.off, p[1]])
    if line.axis == "h":
        return np.array([p[0], line.off])
    return line.a + ((p - line.a) @ line.d) * line.d


def _line_point_dir(line):
    if line.axis == "v":
        return np.array([line.off, 0.0]), np.array([0.0, 1.0])
    if line.axis == "h":
        return np.array([0.0, line.off]), np.array([1.0, 0.0])
    return line.a, line.d


def _junction(p, q, near, tol):
    """Corner between primitives p (before) and q (after), near the old corner."""
    if p.kind == "line" and q.kind == "line":
        a0, d0 = _line_point_dir(p)
        a1, d1 = _line_point_dir(q)
        A = np.c_[d0, -d1]
        if abs(np.linalg.det(A)) < 1e-9:
            return None
        t = np.linalg.solve(A, a1 - a0)
        return a0 + t[0] * d0
    if p.kind == "arc" and q.kind == "arc":
        d = np.linalg.norm(q.c - p.c)
        if d < tol and abs(p.r - q.r) < tol:            # one circle split in two arcs
            v = near - p.c
            return p.c + p.r * v / np.linalg.norm(v)
        if d < 1e-12:
            return None
        u = (q.c - p.c) / d
        if abs(d - (p.r + q.r)) < 2 * tol or abs(d - abs(p.r - q.r)) < 2 * tol:   # tangent arcs
            s = 1.0 if abs(d - (p.r + q.r)) < 2 * tol else (1.0 if p.r > q.r else -1.0)
            return p.c + s * p.r * u
        a = (p.r ** 2 - q.r ** 2 + d ** 2) / (2 * d)
        h2 = p.r ** 2 - a ** 2
        if h2 < 0:
            return None
        w = np.array([-u[1], u[0]])
        cand = [p.c + a * u + s * np.sqrt(h2) * w for s in (1, -1)]
        return min(cand, key=lambda x: np.linalg.norm(x - near))
    ln, ar = (p, q) if p.kind == "line" else (q, p)
    a0, d0 = _line_point_dir(ln)
    foot = a0 + ((ar.c - a0) @ d0) * d0
    h = float(np.linalg.norm(ar.c - foot))
    if abs(h - ar.r) < 2 * tol:                      # tangent: the foot of the perpendicular
        return foot
    if h > ar.r:
        return None
    s = np.sqrt(ar.r ** 2 - h ** 2)
    cand = [foot + s * d0, foot - s * d0]
    return min(cand, key=lambda x: np.linalg.norm(x - near))


def _tangent_line(n0, h0, pts):
    """Line n.x = h tangent to each (point, signed radius) in pts (a radius 0 is a point on the
    line), nearest to the line n0.x = h0; None if there is none. pts: [(c1, s1 r1), (c2, s2 r2)]."""
    (c1, k1), (c2, k2) = pts
    D = c1 - c2
    dd = float(D @ D)
    k = k1 - k2
    if dd < 1e-18 or k * k > dd:
        return None
    w = np.array([-D[1], D[0]]) / np.sqrt(dd)
    base = (k / dd) * D
    t = np.sqrt(max(0.0, 1.0 - k * k / dd))
    best = None
    for n in (base + t * w, base - t * w):
        h = float(n @ c1 - k1)
        score = float(np.linalg.norm(n - n0)) + abs(h - h0)
        if best is None or score < best[0]:
            best = (score, n, h)
    return best[1], best[2]


def _tangent_slanted(pr, tol):
    """A slanted line next to an arc, tangent to it within 2 tol, is made exactly tangent: the
    common tangent of the two circles when arcs are on both sides, else the tangent from the line's
    far end (kept fixed) to the circle. Each circle keeps the side of the line it was on."""
    n = len(pr)
    for m, p in enumerate(pr):
        if p.kind != "line" or p.axis is not None:
            continue
        n0 = np.array([-p.d[1], p.d[0]])
        h0 = float(n0 @ p.a)
        prev, nxt = pr[m - 1], pr[(m + 1) % n]
        touch = []
        for q in (prev, nxt):
            if q.kind == "arc":
                dist = float(n0 @ q.c - h0)
                if abs(abs(dist) - q.r) < 2 * tol:
                    touch.append((q.c, np.sign(dist) * q.r))
                    continue
            touch.append(None)
        if touch == [None, None]:
            continue
        if touch[0] is None:
            touch[0] = (p.a, 0.0)                     # far end fixed
        if touch[1] is None:
            touch[1] = (p.b, 0.0)
        r = _tangent_line(n0, h0, touch)
        if r is None:
            continue
        nn, h = r
        d = np.array([nn[1], -nn[0]])
        if d @ p.d < 0:
            d = -d
        p.a = p.a - (nn @ p.a - h) * nn
        p.b = p.b - (nn @ p.b - h) * nn
        p.d = d


def _dist_to_loop(P, edges):
    """Max over points P of the exact distance to the loop (segments, arcs, circles)."""
    P = np.asarray(P, float)
    best = np.full(len(P), np.inf)
    for e in edges:
        if e[0] == "line":
            a, b = np.asarray(e[1], float), np.asarray(e[2], float)
            d = b - a
            t = np.clip((P - a) @ d / max(d @ d, 1e-30), 0, 1)
            dist = np.linalg.norm(P - (a + t[:, None] * d), axis=1)
        elif e[0] == "arc":
            c, r = _circle3(*e[1:])
            a0 = np.arctan2(*(np.asarray(e[1]) - c)[::-1])
            am = np.arctan2(*(np.asarray(e[2]) - c)[::-1])
            a1 = np.arctan2(*(np.asarray(e[3]) - c)[::-1])
            s1 = (a1 - a0) % (2 * np.pi)
            sm = (am - a0) % (2 * np.pi)
            ccw = sm < s1                                   # arc runs counter-clockwise from a0
            ang = np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0])
            rel = (ang - a0) % (2 * np.pi) if ccw else (a0 - ang) % (2 * np.pi)
            span = s1 if ccw else (2 * np.pi - s1)
            on = rel <= span
            radial = np.abs(np.linalg.norm(P - c, axis=1) - r)
            ends = np.minimum(np.linalg.norm(P - np.asarray(e[1]), axis=1),
                              np.linalg.norm(P - np.asarray(e[3]), axis=1))
            dist = np.where(on, radial, ends)
        elif e[0] == "ellipse":
            from .prism25 import _ellipse_dist
            dist = _ellipse_dist(P, np.asarray(e[1], float), e[2], e[3], e[4])
        else:
            dist = np.abs(np.linalg.norm(P - np.asarray(e[1]), axis=1) - e[2])
        best = np.minimum(best, dist)
    return float(best.max())


def snap_loops(loops, points, tol, log=print):
    """Snapped copies of the loops (prism25.Loop); points[i] are the section points of loop i."""
    from .prism25 import Loop
    eps, reps = 0.5 * tol, 0.25 * tol
    X, Y, R = _Pool(eps, reps), _Pool(eps, reps), _Pool(eps, reps)
    parsed = []
    for L in loops:
        if L.pattern is not None:
            parsed.append(None)
            continue
        if len(L.edges) == 1 and L.edges[0][0] == "circle":
            c, r = L.edges[0][1], L.edges[0][2]
            X.add(c[0]); Y.add(c[1]); R.add(r)
            parsed.append("circle")
            continue
        if len(L.edges) == 1 and L.edges[0][0] == "ellipse":
            _, c, a, b, _ = L.edges[0]
            X.add(c[0]); Y.add(c[1]); R.add(a); R.add(b)
            parsed.append("ellipse")
            continue
        pr = _prims(L)
        if pr is None:
            parsed.append(None)
            continue
        for p in pr:
            if p.kind == "line":
                dev_v = abs(p.b[0] - p.a[0])
                dev_h = abs(p.b[1] - p.a[1])
                if dev_v < tol and dev_v < np.tan(np.radians(2)) * p.L:
                    p.axis, p.off = "v", 0.5 * (p.a[0] + p.b[0])
                    X.add(p.off)
                elif dev_h < tol and dev_h < np.tan(np.radians(2)) * p.L:
                    p.axis, p.off = "h", 0.5 * (p.a[1] + p.b[1])
                    Y.add(p.off)
            else:
                X.add(p.c[0]); Y.add(p.c[1]); R.add(p.r)
        parsed.append(pr)
    grids = choose_grids(X.vals + Y.vals + R.vals, reps)
    tx, ty, tr = X.table(grids), Y.table(grids), R.table(grids)
    if DEBUG:
        print("    snap: grids", "inch" if grids is IN_GRIDS else "mm", "y table", ty)
    out, n_snapped = [], 0
    for L, pr, P in zip(loops, parsed, points):
        if pr is None:
            out.append(L)
            continue
        if pr == "ellipse":
            _, c, a, b, ang = L.edges[0]
            # axis direction: whole degrees, which include the axis-parallel case
            ang_s = np.radians(round(np.degrees(ang)))
            if abs(ang_s - ang) * a > tol / 4:
                ang_s = ang
            new = Loop([("ellipse", np.array([_Pool.get(tx, c[0]), _Pool.get(ty, c[1])]),
                         _Pool.get(tr, a), _Pool.get(tr, b), float(ang_s % np.pi))], L.area)
            if _dist_to_loop(P, new.edges) < tol:
                out.append(new)
                n_snapped += 1
            else:
                out.append(L)
            continue
        if pr == "circle":
            c, r = L.edges[0][1], L.edges[0][2]
            new = Loop([("circle", np.array([_Pool.get(tx, c[0]), _Pool.get(ty, c[1])]), _Pool.get(tr, r))], L.area)
            if _dist_to_loop(P, new.edges) < tol:
                out.append(new)
                n_snapped += 1
            else:
                out.append(L)
            continue
        n = len(pr)
        for p in pr:                                     # snapped values
            if p.kind == "line" and p.axis == "v":
                p.off = _Pool.get(tx, p.off)
            elif p.kind == "line" and p.axis == "h":
                p.off = _Pool.get(ty, p.off)
            elif p.kind == "arc":
                p.c = np.array([_Pool.get(tx, p.c[0]), _Pool.get(ty, p.c[1])])
                p.r = _Pool.get(tr, p.r)
        for m in range(n):                               # tangency to axis-parallel lines
            p = pr[m]
            if p.kind != "arc":
                continue
            for q in (pr[m - 1], pr[(m + 1) % n]):
                if q.kind != "line" or q.axis is None:
                    continue
                k = 0 if q.axis == "v" else 1
                side = np.sign(p.c[k] - q.off) or 1.0
                if abs(abs(p.c[k] - q.off) - p.r) < 2 * tol:
                    p.c[k] = q.off + side * p.r
        _tangent_slanted(pr, tol)                        # slanted lines tangent to arcs
        old = [e[1] for e in L.edges]                    # old corners: start of each edge
        V = []
        for m in range(n):
            v = _junction(pr[m - 1], pr[m], np.asarray(old[m], float), tol)
            # a corner at a tangent junction is poorly placed by the fit (it can slide along the
            # tangent), so it may move further; the point-to-loop check below guards the shape
            lim = 10 * tol
            for q in (pr[m - 1], pr[m]):                 # a tangent corner slides about sqrt(2 r tol)
                if q.kind == "arc" and {pr[m - 1].kind, pr[m].kind} == {"line", "arc"}:
                    lim = max(lim, 3 * np.sqrt(2 * q.r * tol))
            if v is None or np.linalg.norm(v - old[m]) > lim:
                V = None
                break
            V.append(v)
        if V is None:
            if DEBUG:
                print("    snap: corner not found or moved too far", m, pr[m - 1].kind, pr[m].kind, old[m], v)
            out.append(L)
            continue
        edges = []
        for m, p in enumerate(pr):
            a, b = V[m], V[(m + 1) % n]
            if p.kind == "line":
                edges.append(("line", a, b))
            else:
                ang = np.arctan2(*(p.m - p.c)[::-1])     # middle point on the snapped circle
                mid = p.c + p.r * np.array([np.cos(ang), np.sin(ang)])
                edges.append(("arc", a, mid, b))
        dd = _dist_to_loop(P, edges)
        if DEBUG:
            print(f"    snap: loop of {n} edges, max point distance {dd:.4f} (tol {tol:.4f})")
            for e in edges:
                print("      ", e[0], [np.round(np.asarray(q, float), 4).tolist() for q in e[1:]])
        if dd < tol:
            out.append(Loop(edges, L.area))
            n_snapped += 1
        else:
            out.append(L)
    log(f"  snapping: {n_snapped} of {len(loops)} loop(s) snapped to clean geometry")
    return out
