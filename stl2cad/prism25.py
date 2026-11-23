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


def _wall_or_cap(mesh, k, deg=1.0):
    c = np.abs(mesh.face_normals[:, k])
    return (c > np.cos(np.radians(deg))) | (c < np.sin(np.radians(deg)))


def secondary_cylinders(mesh, k, tol, max_fraction=0.1, log=print, ignore=None):
    """Bosses and holes along directions other than the prism axis k. The faces that are neither
    walls nor caps for k must form connected patches, each a cylinder (normals perpendicular to
    one axis, vertices on one circle about it). Returns a list of dict(op, c, d, r, t0, t1, ext)
    with c a point on the axis, d the unit axis (from the attached end to the free end of a boss),
    the feature spanning c + t d for t in [t0, t1], and ext its extent along k; or None if the
    off-axis faces are not all such cylinders."""
    import trimesh
    off = ~_wall_or_cap(mesh, k)
    if ignore is not None:                       # faces of revolved features about axes along k
        off &= ~ignore
    a = mesh.area_faces
    if not off.any():
        return []
    if a[off].sum() > max_fraction * a.sum():
        return None
    idx = np.nonzero(off)[0]
    adj = mesh.face_adjacency
    keep = off[adj[:, 0]] & off[adj[:, 1]]
    comps = trimesh.graph.connected_components(adj[keep], nodes=idx, min_len=1)
    fits = []                                    # (point on axis, unit axis, radius) per patch
    for comp in comps:
        comp = np.asarray(comp)
        n, w = mesh.face_normals[comp], a[comp]
        M = (n * w[:, None]).T @ n
        _, vec = np.linalg.eigh(M)
        d = vec[:, 0]
        if np.abs(n @ d).max() > np.sin(np.radians(2.0)):
            log(f"  2.5D: off-axis patch of {len(comp)} faces is not a cylinder")
            return None
        V = mesh.vertices[np.unique(mesh.faces[comp])]
        u, v = _perp(d)
        c2, r, err = _circle_fit(np.c_[V @ u, V @ v])
        if err > tol:
            log(f"  2.5D: off-axis patch of {len(comp)} faces is not circular (error {err:.3f} mm)")
            return None
        c = c2[0] * u + c2[1] * v
        for f in fits:                           # same cylinder as an earlier patch
            if (np.linalg.norm(np.cross(f[1], d)) < 1e-3 and abs(f[2] - r) < tol
                    and np.linalg.norm(np.cross(c - f[0], f[1])) < tol):
                break
        else:
            fits.append((c, d, float(r)))
    out = []
    for c, d, r in fits:
        if np.abs(d).max() > np.cos(np.radians(0.1)):   # along a model axis
            m_ = int(np.argmax(np.abs(d)))
            d = np.sign(d[m_]) * np.eye(3)[m_] + 0.0
            c = c - (c @ d) * d
            from .snap2d import nice
            c = np.array([nice(float(x), tol / 4) for x in c]) + 0.0
            r = nice(float(r), tol / 4)
        # every facet on this cylinder, including those whose normals happen to be along or across
        # the prism axis (a cross hole's wall facets facing +-x or +-z when the prism axis is y)
        u, v = _perp(d)
        T = mesh.triangles
        rad = T - c - (T @ d)[..., None] * d
        on = (np.abs(np.linalg.norm(rad, axis=2) - r) < tol).all(axis=1)
        n = mesh.face_normals
        on &= np.abs(n @ d) < np.sin(np.radians(2.0))
        rc = rad.mean(axis=1)
        rc /= np.linalg.norm(rc, axis=1)[:, None] + 1e-30
        on &= np.abs(np.einsum("ij,ij->i", rc, n)) > np.cos(np.radians(10.0))
        comp_all = np.nonzero(on)[0]
        # pieces of one cylinder separated along its axis (a stepped cross hole: the wide bore
        # through both lugs, a narrow one between them) are separate features
        tt = (mesh.triangles[comp_all] - c) @ d
        o = np.argsort(tt.min(axis=1))
        pieces, hi = [], -np.inf
        for f_ in o:
            if tt[f_].min() > hi + max(4 * tol, 0.5 * r):
                pieces.append([])
            pieces[-1].append(comp_all[f_])
            hi = max(hi, tt[f_].max())
        for comp in map(np.asarray, pieces):
            w = a[comp]
            V = mesh.vertices[np.unique(mesh.faces[comp])]
            err = float(np.abs(np.linalg.norm(V - c - np.outer(V @ d, d), axis=1) - r).max())
            t = V @ d
            t0, t1 = float(t.min()), float(t.max())
            # outward normals point away from the axis on a boss, towards it in a hole
            F = mesh.triangles_center[comp]
            radial = F - c - np.outer(F @ d, d)
            op = "add" if float(np.sum(w * np.einsum("ij,ij->i", radial, n[comp]))) > 0 else "cut"
            # an end is free when the point just beyond it on the axis lies outside the material.
            # d is oriented so that t0 is the end a designer sketches on: the attached end of a boss,
            # the open end of a hole
            e = max(2 * tol, 0.05 * r)
            inside = mesh.contains([c + (t0 - e) * d, c + (t1 + e) * d])
            if op == "cut":                            # an end is closed if any cap facet lies on its disk
                inside = np.array([inside[0] or _cap_at(mesh, c, d, r, t0, tol),
                                   inside[1] or _cap_at(mesh, c, d, r, t1, tol)])
            through_boss = False
            if op == "add":
                if inside[0] and inside[1]:
                    log("  2.5D: off-axis boss is closed at both ends")
                    return None
                through_boss = not inside[0] and not inside[1]   # passes through the prism body
                flip = bool(inside[1])
            else:
                if inside[0] and inside[1]:
                    log("  2.5D: off-axis hole is closed at both ends")
                    return None
                flip = bool(inside[0])
                if not inside[0] and not inside[1]:      # open at both ends: sketch on a flat end face
                    f0 = _face_point(mesh, c + t0 * d, d, u, v, r, tol)
                    flip = f0 is None and _face_point(mesh, c + t1 * d, d, u, v, r, tol) is not None
            if flip:
                d, t0, t1 = -d, -t1, -t0
                inside = inside[::-1]
            tf = t0
            if op == "cut":                            # cut past the open end(s)
                t0 -= e
                t1 += 0 if inside[1] else e
            face_pt = _face_point(mesh, c + tf * d, d, u, v, r, tol)
            if through_boss:                           # sketched on a face of the body it crosses
                face_pt = None
                for tc in _cap_heights(mesh, c, d, r, t0, t1, tol):
                    face_pt = _face_point(mesh, c + tc * d, d, u, v, r, tol)
                    if face_pt is not None:
                        tf = tc
                        break
            ext = (float(V[:, k].min()), float(V[:, k].max()))
            out.append(dict(op=op, c=c, d=d, r=float(r), t0=t0, t1=t1, tf=tf, face_pt=face_pt, spans=through_boss,
                            through=bool(op == "cut" and not inside[1]), ext=ext, err=float(err),
                            faces=comp))
    return out


def _cap_heights(mesh, c, d, r, t0, t1, tol):
    """Heights t in (t0, t1) of planar faces normal to d that meet the circle of radius r about
    the axis from outside (faces of a body the cylinder passes through), largest area first."""
    F = mesh.triangles_center
    cap = np.abs(mesh.face_normals @ d) > np.cos(np.radians(1.0))
    t = (F - c) @ d
    rad = np.linalg.norm(F - c - np.outer(t, d), axis=1)
    sel = cap & (t > t0 + tol) & (t < t1 - tol) & (rad > r) & (rad < 4 * r)
    hs = {}                                      # grid key -> [area, sum of area * t]
    for tt, w in zip(t[sel], mesh.area_faces[sel]):
        g = hs.setdefault(int(round(tt / tol)), [0.0, 0.0])
        g[0] += w
        g[1] += w * tt
    return [g[1] / g[0] for g in sorted(hs.values(), key=lambda g: -g[0])]


def _cap_at(mesh, c, d, r, t, tol):
    """Some facet with normal along d lies on the disk of radius r about the axis at c + t d."""
    F = mesh.triangles_center
    cap = np.abs(mesh.face_normals @ d) > np.cos(np.radians(1.0))
    on = np.abs((F - c) @ d - t) < tol
    rad = np.linalg.norm(F - c - np.outer((F - c) @ d, d), axis=1) < r - tol
    return bool((cap & on & rad).any())


def _perp(d):
    u = np.cross(d, [1.0, 0, 0] if abs(d[0]) < 0.9 else [0, 1.0, 0])
    u /= np.linalg.norm(u)
    return u, np.cross(d, u)


def _face_point(mesh, p, d, u, v, r, tol):
    """A point on the face the feature starts from, just outside its circle, on a mesh facet whose
    normal is along d (the face a designer sketches on); None if there is no such facet."""
    for ang in np.radians(np.arange(0, 360, 30)):
        q = p + (r + max(4 * tol, 0.1 * r)) * (np.cos(ang) * u + np.sin(ang) * v)
        cp, dist, tri = mesh.nearest.on_surface([q])
        if dist[0] < tol and abs(mesh.face_normals[tri[0]] @ d) > np.cos(np.radians(1.0)):
            return cp[0]
    return None


def _free_height(a, b, feats, tol):
    """Height in (a, b) away from the k-extent of every secondary feature (the section there is
    that of the prism alone), at the middle of the longest free interval; None if none is free."""
    cuts = sorted((max(a, lo - tol), min(b, hi + tol)) for lo, hi in (f["ext"] for f in feats)
                  if hi + tol > a and lo - tol < b)
    best, cur = None, a
    for lo, hi in cuts + [(b, b)]:
        if lo - cur > 2 * tol and (best is None or lo - cur > best[1] - best[0]):
            best = (cur, lo)
        cur = max(cur, hi)
    return None if best is None else 0.5 * (best[0] + best[1])


def _open_section(loops, span, k, z, tol):
    """Section loops with the slivers left by boss chords removed: the section polygon is opened
    by the widest chord, and loop points away from the opened polygon are dropped."""
    import shapely
    w = 0.0
    for c in span:
        dk = abs(z - c["ext"][0]), abs(c["ext"][1] - z)
        sag = min(dk)
        w = max(w, np.sqrt(max(0.0, 2 * c["r"] * sag)))
    rad = w + 0.25 * tol
    out = []
    for pts, poly in loops:
        pg = shapely.Polygon(poly).buffer(0)
        op = pg.buffer(-rad, join_style="mitre").buffer(rad, join_style="mitre")
        if op.is_empty:
            continue
        keep = shapely.distance(shapely.points(pts), op.boundary) < tol
        q = pts[keep]
        if len(q) >= 3:
            out.append((q, np.asarray(op.exterior.coords) if op.geom_type == "Polygon" else poly))
    return out


def cylinder_solid(cy):
    import cadquery as cq
    c, d = np.asarray(cy["c"]), np.asarray(cy["d"])
    base = c + cy["t0"] * d
    return cq.Solid.makeCylinder(cy["r"], cy["t1"] - cy["t0"], cq.Vector(*base), cq.Vector(*d))


def levels(mesh, k, tol, skip=()):
    """Heights of the cap faces along axis k, clustered within tol; faces in skip (facets of
    off-axis cylinders that happen to face along k) are not caps."""
    cap = np.abs(mesh.face_normals[:, k]) > np.cos(np.radians(1.0))
    cap[np.asarray(list(skip), int)] = False
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
    from .snap2d import nice
    return [nice(v, 1e-6) for v in out]          # 3e-17 is 0


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
    # a wall facet's vertices project onto the section curve; a tilted facet's (a badly shaped
    # triangle where a cross hole meets a wall) belong to other heights: its crossing point is used
    tilted = np.abs(mesh.face_normals[fidx, k]) > 1e-3
    out = []
    for order in loops:
        poly = np.array([S[s, e] for s, e in order])
        pts, par = [], []
        for n_, (s, e) in enumerate(order):
            a, b = S[s, e], S[s, 1 - e]
            d = b - a
            L2 = float(d @ d) + 1e-30
            if tilted[s]:
                pts.append(a)
                par.append(float(n_))
                continue
            for v in V[F[fidx[s]]]:
                t = np.clip((v - a) @ d / L2, 0.0, 1.0)
                if np.linalg.norm(a + t * d - v) < 0.05 * np.sqrt(L2) + 4 * tol:
                    pts.append(v)
                    par.append(n_ + t)
        pts, par = np.array(pts), np.array(par)
        o = np.argsort(par, kind="stable")
        pts = pts[o]
        # each vertex once, at its first place: on small fillets a facet's vertex is also picked
        # up from the neighbouring segment, after the next vertex, which would fold the loop back
        keep, seen = [], set()
        for m in range(len(pts)):
            kk = key(pts[m] / 10)
            if kk not in seen:
                seen.add(kk)
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


def ellipse_points(c, a, b, ang, n):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    u = np.array([np.cos(ang), np.sin(ang)])
    v = np.array([-u[1], u[0]])
    return np.asarray(c) + np.outer(a * np.cos(t), u) + np.outer(b * np.sin(t), v)


def _ellipse_dist(P, c, a, b, ang):
    """Distance from points to an ellipse: nearest of dense samples, refined by Newton steps on
    the parameter (exact to well below the fit tolerance)."""
    u = np.array([np.cos(ang), np.sin(ang)])
    v = np.array([-u[1], u[0]])
    X = (np.asarray(P) - c) @ u
    Y = (np.asarray(P) - c) @ v
    t0 = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    t = t0[np.argmin((X[:, None] - a * np.cos(t0)) ** 2 + (Y[:, None] - b * np.sin(t0)) ** 2, axis=1)]
    for _ in range(8):                           # minimise |(X, Y) - (a cos t, b sin t)|^2
        ex, ey = a * np.cos(t), b * np.sin(t)
        dx, dy = -a * np.sin(t), b * np.cos(t)
        g = (ex - X) * dx + (ey - Y) * dy
        h = dx * dx + dy * dy + (ex - X) * (-a * np.cos(t)) + (ey - Y) * (-b * np.sin(t))
        t = t - g / np.where(np.abs(h) > 1e-12, h, 1e-12)
    return np.hypot(X - a * np.cos(t), Y - b * np.sin(t))


def _full_ellipse(P, tol):
    """Ellipse through a closed loop (centre, semi-axes a >= b, angle of the major axis), or None.
    A general conic is fitted, it must be an ellipse that is not a circle, every point must lie
    within tol of it, and the loop must go once around its centre without large gaps."""
    if len(P) < 8:
        return None
    m = P.mean(0)
    s = np.abs(P - m).max() or 1.0
    x, y = ((P - m) / s).T
    A = np.c_[x * x, x * y, y * y, x, y, np.ones(len(x))]
    _, _, vt = np.linalg.svd(A)
    a_, b_, c_, d_, e_, f_ = vt[-1]
    if b_ * b_ - 4 * a_ * c_ >= 0:
        return None
    x0, y0 = np.linalg.solve([[2 * a_, b_], [b_, 2 * c_]], [-d_, -e_])
    f0 = a_ * x0 * x0 + b_ * x0 * y0 + c_ * y0 * y0 + d_ * x0 + e_ * y0 + f_
    lam, vec = np.linalg.eigh([[a_, b_ / 2], [b_ / 2, c_]])
    ax2 = -f0 / lam
    if (ax2 <= 0).any():
        return None
    semi = np.sqrt(ax2) * s
    i = int(np.argmax(semi))
    a, b = float(semi[i]), float(semi[1 - i])
    ang = float(np.arctan2(vec[1, i], vec[0, i]))
    c = m + s * np.array([x0, y0])
    if a < 1.01 * b:                             # a circle (handled by the circle test)
        return None
    if _ellipse_dist(P, c, a, b, ang).max() > tol:
        return None
    Q = np.vstack([P, P[:1]])
    w = np.unwrap(np.arctan2(Q[:, 1] - c[1], Q[:, 0] - c[0]))
    if abs(abs(w[-1] - w[0]) - 2 * np.pi) > 0.2:
        return None
    if np.diff(np.sort(np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]))).max(initial=0) > np.radians(30):
        return None
    return dict(c=c, a=a, b=b, ang=ang % np.pi)


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
    if cosang[s] > np.cos(np.radians(30.0)):      # no sharp corner (every corner filleted):
        s = int(np.argmax(np.linalg.norm(b, axis=1)))   # start of the longest step, a tangent point
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


def _mirror_axes(P, tol):
    """Sketch axes through the origin (0: the line a = 0, 1: the line b = 0) about which the closed
    point loop P is mirror-symmetric within tol."""
    import shapely
    ring = shapely.LinearRing(P)
    out = []
    for ax in (0, 1):
        M = np.array([1.0, 1.0])
        M[ax] = -1.0
        if np.max(shapely.distance(shapely.points(P * M), ring)) < tol:
            out.append(ax)
    return out


def _half_chain(P, ax):
    """The mesh points of closed loop P on the side a >= 0 (ax 0) or b >= 0 (ax 1), in loop order
    between the loop's two crossings of the axis (no interpolated points: a crossing interpolated
    along a facet chord lies off a curved edge); None unless the loop crosses the axis twice."""
    c = P[:, ax]
    n = len(P)
    pos = c >= 0
    ups = [i for i in range(n) if not pos[i] and pos[(i + 1) % n]]
    downs = [i for i in range(n) if pos[i] and not pos[(i + 1) % n]]
    if len(ups) != 1 or len(downs) != 1:
        return None
    i, j = ups[0], downs[0]
    idx = [(i + 1 + k) % n for k in range((j - i) % n)]
    Q = P[idx].copy()
    Q[np.abs(Q[:, ax]) < 1e-9, ax] = 0.0
    return Q if len(Q) >= 3 else None


def _chain_edges(Q, ax_start, ax_end, tol):
    """Lines and arcs through the open chain Q whose first point lies on axis ax_start and last on
    ax_end (the lines a = 0 or b = 0). The chain is extended past both axes by its own mirror
    image, so that a curve crossing an axis is fitted with support on both sides; the fitted
    pieces are then cut where they cross the axes."""
    def mir(X, ax):
        M = np.array([1.0, 1.0])
        M[ax] = -1.0
        return X * M
    k = min(len(Q) - 2, max(4, len(Q) // 3))
    # across the axis the loop continues with the mirror images of the chain's own points; a first
    # (last) point on the axis is its own image and is not repeated
    on0, on1 = Q[0][ax_start] == 0.0, Q[-1][ax_end] == 0.0
    pre = mir(Q[1:k + 1] if on0 else Q[:k], ax_start)[::-1]
    post = mir(Q[-k - 1:-1] if on1 else Q[-k:], ax_end)[::-1]
    Qe = np.vstack([pre, Q, post])
    s0, e0 = len(pre), len(pre) + len(Q) - 1
    pieces = fit_chain(Qe, tol)
    if not pieces:
        return None

    def axis_piece(ax):
        d = np.zeros(2)
        d[1 - ax] = 1.0
        return ("line", 0, 0, dict(c=np.zeros(2), d=d))

    def cut(piece, near, ax):
        gap = abs(near[ax])                       # the last mesh point may lie well off the axis
        near = np.array(near, float)
        near[ax] = 0.0
        x = _intersect(piece, axis_piece(ax), near)
        if x is None or np.linalg.norm(x - near) > 3 * tol + gap:
            x = _project(near, piece).copy()
        x = np.array(x, float)
        x[ax] = 0.0
        return x
    # the axes cross the chain between points s0 - 1 and s0, and between e0 and e0 + 1 (or at
    # s0 / e0 themselves when those lie on the axis)
    sa = s0 if Q[0][ax_start] == 0.0 else s0 - 1
    eb = e0 if Q[-1][ax_end] == 0.0 else e0 + 1
    first = next(m for m, p in enumerate(pieces) if p[2] > sa)      # piece reaching past the start axis
    last = next(m for m, p in enumerate(pieces) if p[2] >= eb)       # piece reaching the end axis
    sel = pieces[first:last + 1]
    V = [cut(sel[0], Qe[s0], ax_start)]
    for m in range(1, len(sel)):
        V.append(_junction(sel[m - 1], sel[m], Qe[sel[m][1]], tol))
    V.append(cut(sel[-1], Qe[e0], ax_end))
    edges = []
    for m, p in enumerate(sel):
        a, b = V[m], V[m + 1]
        if np.linalg.norm(b - a) < 1e-9:
            continue
        if p[0] == "line":
            edges.append(("line", a, b))
        else:
            lo, hi = max(p[1], s0), min(p[2], e0)
            mid = _project(Qe[(lo + hi) // 2], p)
            if np.linalg.norm(mid - a) < 1e-9 or np.linalg.norm(mid - b) < 1e-9:
                mid = _project(0.5 * (a + b), p)
            edges.append(("arc", a, mid, b))
    return edges


def _mirror_edges(edges, ax):
    """The mirror image of an edge chain about axis ax, traversed backwards (it continues the
    chain on the other side)."""
    M = np.array([1.0, 1.0])
    M[ax] = -1.0
    return [(e[0],) + tuple(np.asarray(q) * M for q in e[1:][::-1]) for e in edges[::-1]]


def _join(e0, e1, tol=1e-6):
    """One edge for two edges meeting at a point on a symmetry axis, if they are one line or one
    arc (the meeting point is then the arc's middle); else None."""
    if e0[0] == "line" and e1[0] == "line":
        u, v = e0[2] - e0[1], e1[2] - e1[1]
        # mirror images meet at the axis; nearly collinear (within 2 deg) means one line across it
        if abs(u[0] * v[1] - u[1] * v[0]) < np.sin(np.radians(2.0)) * np.linalg.norm(u) * np.linalg.norm(v) and u @ v > 0:
            return ("line", e0[1], e1[2])
    if e0[0] == "arc" and e1[0] == "arc":
        c0, c1 = _centre3(*e0[1:]), _centre3(*e1[1:])
        if c0 is None or c1 is None:
            return None
        r0, r1 = np.linalg.norm(e0[1] - c0), np.linalg.norm(e1[1] - c1)
        if np.linalg.norm(c0 - c1) < tol and abs(r0 - r1) < tol:
            return ("arc", e0[1], e0[3], e1[3])      # through the two ends and the axis point
    return None


def _centre3(p, q, r):
    """Centre of the circle through three points, or None if they are collinear."""
    p, q, r = (np.asarray(v, float) for v in (p, q, r))
    den = 2 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1]))
    if abs(den) < 1e-12:
        return None
    return np.array([((p @ p) * (q[1] - r[1]) + (q @ q) * (r[1] - p[1]) + (r @ r) * (p[1] - q[1])) / den,
                     ((p @ p) * (r[0] - q[0]) + (q @ q) * (p[0] - r[0]) + (r @ r) * (q[0] - p[0])) / den])


def _glue(A, B, tol=1e-6):
    """Chain A followed by chain B (B starts where A ends, on an axis), joining the two edges at
    the junction when they are one curve."""
    j = _join(A[-1], B[0], tol) if A and B else None
    return A[:-1] + [j] + B[1:] if j is not None else A + B


def mirrored_loop(P, tol):
    """A loop symmetric about one or both sketch axes through the origin, fitted on one half (or
    one quarter) and mirrored, so that mirrored edges are exact mirror images; None if P is not
    symmetric or the result is not within tol of P."""
    axes = _mirror_axes(P, tol)
    if not axes:
        return None
    from .snap2d import _dist_to_loop
    a0 = axes[0]
    H = _half_chain(P, a0)
    if H is None:
        return None
    half = None
    if len(axes) == 2:
        a1 = axes[1]
        c = H[:, a1]
        # the half chain runs from one end on axis a0 to the other, crossing axis a1 once
        k = [m for m in range(len(H) - 1) if (c[m] < 0) != (c[m + 1] < 0)]
        if len(k) == 1:
            Q = H[k[0] + 1:]
            if len(Q) >= 3:
                quarter = _chain_edges(Q, a1, a0, tol)
                if quarter:
                    half = _glue(_mirror_edges(quarter, a1), quarter, tol)
    if half is None:
        half = _chain_edges(H, a0, a0, tol)
    if not half:
        return None
    loop = _glue(half, _mirror_edges(half, a0), tol)
    j = _join(loop[-1], loop[0], tol)                # the closing junction, on axis a0
    if j is not None:
        loop = [j] + loop[1:-1]
    for m in range(len(loop)):                       # identical shared vertices
        nxt = loop[(m + 1) % len(loop)]
        loop[m] = loop[m][:-1] + (np.asarray(nxt[1]),)
    if len(loop) < 2 or _dist_to_loop(P, loop) > tol:
        return None
    return Loop(loop, _area(P))


MIN_PATTERN = 6                       # rotational patterns of fewer repeats are fitted as plain loops


class Loop:
    """A fitted closed loop: list of edges ("line", a, b) / ("arc", a, mid, b) / ("circle", c, r).
    pattern: (centre, N, edges of one period) when the loop repeats N times about a centre."""

    def __init__(self, edges, area, pattern=None, linear=None):
        self.edges, self.area, self.pattern = edges, area, pattern
        self.linear = linear          # (p, N, step): edges[:p N] are N copies of edges[:p] stepped


def resymmetrize(L, P, tol):
    """After snapping: a loop whose corners and arc middles still pair up with their mirror
    images about a sketch axis (within tol) is made exactly symmetric by averaging each point
    with its partner's mirror image (snapping fixes junctions one at a time, from data that is
    not exactly symmetric). Kept only if the loop stays within tol of P."""
    from .snap2d import _dist_to_loop
    E = L.edges
    if L.pattern is not None or len(E) < 2 or any(e[0] not in ("line", "arc") for e in E):
        return L
    pts = [[np.asarray(q, float) for q in e[1:]] for e in E]
    kinds = [e[0] for e in E]
    n = len(E)
    changed = False
    for ax in (0, 1):
        M = np.array([1.0, 1.0])
        M[ax] = -1.0
        C = np.array([q[0] for q in pts])
        D = np.linalg.norm(C[None, :, :] - (C * M)[:, None, :], axis=2)
        pc = D.argmin(axis=1)
        if D[np.arange(n), pc].max() > tol or len(set(pc.tolist())) != n:
            continue
        # edge i (corner i -> i+1) mirrors to the edge between the partners, traversed backwards
        pe = [(pc[(i + 1) % n]) % n for i in range(n)]
        if any(kinds[i] != kinds[pe[i]] for i in range(n)):
            continue
        newC = [(C[i] + C[pc[i]] * M) / 2 for i in range(n)]
        newM = {i: (pts[i][1] + pts[pe[i]][1] * M) / 2 for i in range(n) if kinds[i] == "arc"}
        for i in range(n):
            pts[i][0] = newC[i]
            pts[i][-1] = newC[(i + 1) % n]
            if i in newM:
                pts[i][1] = newM[i]
        changed = True
    if not changed:
        return L
    edges = [(kinds[i],) + tuple(pts[i]) for i in range(n)]
    if _dist_to_loop(P, edges) > tol:
        return L
    return Loop(edges, L.area)


MIN_LINEAR = 3                        # linear patterns of fewer copies are written out


def linear_pattern(L, tol):
    """The loop with its longest run of edges that repeats by a translation found and made exact:
    the edges are rotated so the run comes first, as N copies of a period of p edges, each the
    previous one moved by a step vector; all copies are made identical to the first. Returns the
    loop unchanged if no run has MIN_LINEAR or more copies."""
    E = L.edges
    n = len(E)
    if L.pattern is not None or n < 2 * MIN_LINEAR or any(e[0] not in ("line", "arc") for e in E):
        return L
    A = [np.array([np.asarray(q, float) for q in e[1:]]) for e in E]
    kinds = [e[0] for e in E]
    eps = 0.05 * tol

    def same(a, b, v):
        return kinds[a] == kinds[b] and np.abs(A[b] - A[a] - v).max() < eps
    best = None                                   # (coverage, -p, i, p, N, v)
    for p in range(1, n // MIN_LINEAR + 1):
        for i in range(n):
            v = A[(i + p) % n][0] - A[i][0]
            if np.linalg.norm(v) < 2 * tol:
                continue
            run = 0
            while run < n - p and same((i + run) % n, (i + run + p) % n, v):
                run += 1
            N = run // p + 1
            if N >= MIN_LINEAR and N * p < n and (best is None or (N * p, -p) > best[:2]):
                best = (N * p, -p, i, p, N, v)
    if best is None:
        return L
    _, _, i, p, N, v = best
    E = E[i:] + E[:i]
    A = A[i:] + A[:i]
    # least-squares step (to 1e-6 mm, so a round pitch stays round); the first copy, with its
    # snapped values, is the period
    v = np.mean([A[m * p + j][0] - A[(m - 1) * p + j][0] for m in range(1, N) for j in range(p)], axis=0)
    v = np.round(v, 6)
    period = [A[j].copy() for j in range(p)]
    for j in range(1, p):                         # shared vertices within the period
        period[j][0] = period[j - 1][-1]
    period[-1][-1] = period[0][0] + v
    K = [kinds[(i + j) % n] for j in range(n)]          # kinds in the rotated order
    # a part of a period after the last copy (the rack's last tooth, followed by the long flat
    # instead of a gap) and before the first copy are copies too, made exact like the others
    rest_n = n - N * p
    after = 0
    while after < min(p - 1, rest_n - 1) and K[N * p + after] == K[after] and \
            np.abs(A[N * p + after] - period[after] - N * v).max() < eps:
        after += 1
    before = 0
    while before < min(p - 1, rest_n - after - 1) and K[n - 1 - before] == K[p - 1 - before] and \
            np.abs(A[n - 1 - before] - period[p - 1 - before] + v).max() < eps:
        before += 1
    edges = [(K[j],) + tuple(q + m * v for q in period[j]) for m in range(N) for j in range(p)]
    edges += [(K[j],) + tuple(q + N * v for q in period[j]) for j in range(after)]
    rest = [list(e) for e in E[N * p + after:n - before]]
    tail = [(K[p - before + j],) + tuple(q - v for q in period[p - before + j]) for j in range(before)]
    start_rest = (period[after][0] if after < p else period[0][0] + v) + N * v
    rest[0][1] = start_rest
    rest[-1][-1] = tail[0][1] if tail else period[0][0]
    edges += [tuple(e) for e in rest] + tail
    return Loop(edges, L.area, linear=(p, N, v))


def build_loop(P, tol):
    if _full_circle(P, tol) is None:
        el = _full_ellipse(P, tol)
        if el is not None:
            return Loop([("ellipse", el["c"], el["a"], el["b"], el["ang"])], _area(P))
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
    if pat is None:
        mir = mirrored_loop(P, tol)
        if mir is not None:
            return mir
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


def exact_arcs(L, eps=1e-6):
    """Arc middle points put exactly on the arc's circle when its centre and radius are clean
    values within eps (a snapped arc keeps a fitted middle point a few 1e-9 off, and three
    rounded points then give a radius off by 1e-7: OCC fails between coaxial cylinders then)."""
    from .snap2d import nice
    E = list(L.edges)
    changed = False
    for m, e in enumerate(E):
        if e[0] != "arc":
            continue
        a, mid, b = (np.asarray(q, float) for q in e[1:])
        c = _centre3(a, mid, b)
        if c is None:
            continue
        c2 = np.array([nice(float(c[0]), eps), nice(float(c[1]), eps)])
        r = nice(float(np.linalg.norm(a - c)), eps)
        if abs(np.linalg.norm(a - c2) - r) > eps or abs(np.linalg.norm(b - c2) - r) > eps:
            continue
        u = mid - c2
        E[m] = ("arc", e[1], c2 + r * u / np.linalg.norm(u), e[3])
        changed = True
    if not changed:
        return L
    out = Loop(E, L.area, pattern=L.pattern, linear=L.linear)
    return out


def absorb_chords(L, tol):
    """Lines that are chords of a neighbouring arc's circle (both ends on it within tol / 4, the
    line within 1.5 tol of it: no data can tell them from the arc) and neighbouring arcs on one
    circle are merged into one arc. Where a section's points are sparse (a few mm apart on a
    large radius), chords fit within tol and would leave a faceted wall."""
    if L.pattern is not None or L.linear is not None or len(L.edges) < 3:
        return L
    if any(e[0] not in ("line", "arc") for e in L.edges):
        return L
    E = [tuple(np.asarray(q, float) if m else q for m, q in enumerate(e)) for e in L.edges]
    cr = lambda u, v: u[0] * v[1] - u[1] * v[0]

    def circle(e):
        c = _centre3(*e[1:])
        if c is None:
            return None
        return c, float(np.linalg.norm(e[1] - c)), np.sign(cr(e[1] - c, e[2] - c) + cr(e[2] - c, e[3] - c))

    def fits(e, c, R, s):
        a, b = e[1], e[-1]
        if abs(np.linalg.norm(a - c) - R) > tol / 4 or abs(np.linalg.norm(b - c) - R) > tol / 4:
            return False
        if e[0] == "arc":
            g = circle(e)
            return g is not None and np.linalg.norm(g[0] - c) < tol / 10 and abs(g[1] - R) < tol / 10 and g[2] == s
        d = b - a
        dist = abs(cr(d, c - a)) / (np.linalg.norm(d) + 1e-30)
        return -tol / 4 <= R - dist < 1.5 * tol and np.sign(cr(a - c, b - c)) == s

    def arc(a, b, c, R, s):
        ta, tb = np.arctan2(*(a - c)[::-1]), np.arctan2(*(b - c)[::-1])
        sw = (tb - ta) % (2 * np.pi) if s > 0 else -((ta - tb) % (2 * np.pi))
        if abs(sw) > 2 * np.pi - 0.05:
            return None
        tm = ta + sw / 2
        return ("arc", a, c + R * np.array([np.cos(tm), np.sin(tm)]), b)
    changed = True
    while changed and len(E) > 2:
        changed = False
        n = len(E)
        for m in range(n):
            if E[m][0] != "arc":
                continue
            g = circle(E[m])
            if g is None:
                continue
            c, R, s = g
            nx, pv = (m + 1) % n, (m - 1) % n
            if fits(E[nx], c, R, s):
                new = arc(E[m][1], E[nx][-1], c, R, s)
                if new is not None:
                    E[m] = new
                    del E[nx]
                    changed = True
                    break
            if fits(E[pv], c, R, s):
                new = arc(E[pv][1], E[m][-1], c, R, s)
                if new is not None:
                    E[m] = new
                    del E[pv]
                    changed = True
                    break
    if len(E) == len(L.edges):
        return L
    return Loop(E, L.area)


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
        elif e[0] == "ellipse":
            c = np.asarray(e[1], float)
            xdir = V(c + np.array([np.cos(e[4]), np.sin(e[4])])) - V(c)
            es.append(cq.Edge.makeEllipse(e[2], e[3], V(c), cq.Vector(*D), xdir))
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
    from .core import closest_distance
    d1 = closest_distance(tm, V[idx])
    S, _ = trimesh.sample.sample_surface(tm, n, seed=seed)
    d2 = closest_distance(mesh, S)
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
    feats, cyls, revs = [], [], []
    for r in res:
        feats += r.extras["prism25"]["features"]
        cyls += r.extras["prism25"].get("cylinders", [])
        revs += r.extras["prism25"].get("revolved", [])
    out = Result(shape=solid, mode="exact")
    out.extras = dict(prims=[], analytic_fraction=1.0, summary=f"2.5D: {len(feats)} extrusions",
                      volume_error=max(r.extras["volume_error"] for r in res), brep=None, valid=True,
                      passed=True, overlapping_shells=len(bodies),
                      prism25=dict(axis=res[0].extras["prism25"]["axis"], features=feats,
                                   cylinders=cyls, revolved=revs))
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
    from . import revfeat
    k, frac = prism_axis(mesh, min_fraction=0.0)
    groups, member = [], None
    if frac < 0.999:                             # surfaces of revolution about axes along k
        groups, member = revfeat.find_groups(mesh, k, tol, log=log)
    cyls = [] if frac >= 0.999 else secondary_cylinders(mesh, k, tol, log=log, ignore=member)
    if cyls is None and groups:
        log(f"  2.5D along {AX[k]}: {len(groups)} revolved group(s), but the other off-axis faces "
            "are not cylinders")
    if cyls is None:                             # another axis, whose off-axis faces are cylinders
        groups = []
        quiet = lambda *a, **kw: None
        for k2 in range(3):
            f2 = float(mesh.area_faces[_wall_or_cap(mesh, k2)].sum() / mesh.area_faces.sum())
            if k2 == k or f2 < 0.7:
                continue
            c2 = secondary_cylinders(mesh, k2, tol, max_fraction=0.3, log=quiet)
            if c2 is not None:
                k, frac, cyls = k2, f2, c2
                break
    if cyls is None:
        log(f"  2.5D: no prism axis (best face fraction {frac:.3f})")
        return None
    if cyls:
        log(f"  2.5D along {AX[k]} ({100 * frac:.1f} % of the area) with {len(cyls)} off-axis "
            f"cylinder(s): " + ", ".join(f"{'boss' if c['op'] == 'add' else 'hole'} d={2 * c['r']:.3f}"
                                         for c in cyls))
    if groups:
        log(f"  2.5D along {AX[k]} with {len(groups)} revolved group(s) about axes along {AX[k]}")
    skip = [f for c in cyls for f in c["faces"]] + list(revfeat.cap_skip(mesh, k, tol, groups))
    lv = levels(mesh, k, tol, skip=skip)
    slabs = [(a, b) for a, b in zip(lv[:-1], lv[1:]) if b - a > 2 * tol]
    if not slabs:
        return None
    revfeat.set_modes(groups, slabs, tol)
    local = [dict(ext=iv) for g in groups if g["mode"] == "local" for iv in g["spans"]]

    def disks(sec):
        out = revfeat.disk_sections(sec, groups, k, tol)
        if out is None:
            raise ValueError("a section loop crosses the circle of a revolved disk")
        return out
    zs, secs = [], []
    for a, b in slabs:
        z = _free_height(a, b, cyls + local, tol)
        if z is not None:
            zs.append(z)
            secs.append(disks(section_loops(mesh, k, z, tol)))
            continue
        if local and _free_height(a, b, local, tol) is None:
            log("  2.5D: a revolved feature spans a whole slab; no clean section")
            return None
        # bosses that cross the whole slab: section next to the slab end where their chords
        # are narrowest, and drop the slivers the chords leave (the prism under a boss is
        # unknown; any profile within the boss is the same part)
        span = [c for c in cyls if c["op"] == "add" and c["ext"][0] <= a + tol and c["ext"][1] >= b - tol]
        rest = [c for c in cyls if all(c is not s_ for s_ in span)]
        z = None
        for zc in (b - 0.01 * tol, a + 0.01 * tol):
            if span and all(not (lo - tol < zc < hi + tol) for lo, hi in (c["ext"] for c in rest)):
                z = zc
                break
        if z is None:
            log("  2.5D: an off-axis feature spans a whole slab; no clean section")
            return None
        zs.append(z)
        secs.append(disks(_open_section(section_loops(mesh, k, z, tol), span, k, z, tol)))
    tracks = loop_tracks(slabs, secs, tol)
    # material is built in nesting order: outlines, then the holes in them, then islands in the
    # holes, and so on; at every height this reproduces the section
    tracks.sort(key=lambda t: (t["order"], t["op"] == "cut"))
    if not tracks or tracks[0]["op"] != "add":
        return None
    loops = [build_loop(t["P"], tol) for t in tracks]
    # section points are exact mesh vertices: a fit that cuts corners by more than tol / 2 has
    # absorbed small fillets into lines; refit finer and keep that if it is close and not much
    # longer
    from .snap2d import _dist_to_loop
    for m, L in enumerate(loops):
        if L.pattern is not None or len(L.edges) == 1:
            continue
        P = tracks[m]["P"]
        if _dist_to_loop(P, L.edges) <= tol / 2:
            continue
        Lf = build_loop(P, tol / 4)
        if len(Lf.edges) <= 1.75 * len(L.edges) and _dist_to_loop(P, Lf.edges) <= tol / 4:
            log(f"  loop {m + 1}: refitted at tol/4 ({len(L.edges)} -> {len(Lf.edges)} edges)")
            loops[m] = Lf
    loops = [absorb_chords(L, tol) for L in loops]
    loops = consolidate(loops, tol)
    if snap:
        from .snap2d import snap_loops
        loops = snap_loops(loops, [t["P"] for t in tracks], tol, log=log)
        loops = [resymmetrize(L, t["P"], tol) for L, t in zip(loops, tracks)]
    loops = [linear_pattern(exact_arcs(L), tol) for L in loops]
    loops = revfeat.exact_circles(loops, tracks, groups)
    for L in loops:
        if L.linear is not None:
            p_, N_, v_ = L.linear
            log(f"  linear pattern: {N_} copies of {p_} edge(s), step ({v_[0]:.4f}, {v_[1]:.4f}) mm")
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
    revs = []
    if groups:
        revs = revfeat.features(mesh, k, tol, groups, slabs, secs, snap=snap, log=log, prism=solid)
        if revs is None:
            return None
        for f in revs:                           # a disk boss only touches the prism on its cylinder
            if f["op"] == "cut":
                solid = solid.cut(f["solid"])
            else:
                solid = solid.fuse(f["solid"])
    for cy in cyls:                              # off-axis bosses first, then the holes
        if cy["op"] == "add":
            solid = solid.fuse(cylinder_solid(cy))
    for cy in cyls:
        if cy["op"] == "cut":
            solid = solid.cut(cylinder_solid(cy))
    solid = solid.clean()
    n_edges = sum(len(f["loop"].edges) for f in feats)
    from .core import mass_props
    vol_err = abs(mass_props(solid)[0] / mesh.volume - 1)
    valid = solid.isValid() and len(solid.Solids()) == 1
    dmax, dback = _distance_check(solid, mesh) if valid else (np.inf, np.inf)
    n_add = sum(f["op"] == "add" for f in feats)
    log(f"  2.5D along {AX[k]}: {len(lv)} levels, {n_add} extrusion(s) and {len(feats) - n_add} cut(s), "
        f"{len(revs)} revolved feature(s), "
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
                      prism25=dict(axis=k, features=feats, cylinders=_order_cyls(cyls),
                                   revolved=[{q: v for q, v in f.items() if q != "solid"} for f in revs]))
    return res


def _order_cyls(cyls):
    return [dict(op=c["op"], c=c["c"], d=c["d"], r=c["r"], t0=c["t0"], t1=c["t1"], tf=c["tf"],
                 face_pt=c["face_pt"], through=c["through"], spans=c.get("spans", False))
            for c in sorted(cyls, key=lambda c: c["op"] == "cut")]
