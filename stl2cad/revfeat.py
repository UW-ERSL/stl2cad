"""Revolved features on a 2.5D part: surfaces of revolution about axes parallel to the prism axis.

A plate or ring whose walls and caps are those of a stack of prisms may also carry chamfers,
countersinks, rounded edges about a hole or a boss, and formed dishes: faces that are neither walls
nor caps, but whose normals all lie in the planes through one axis parallel to the prism axis. Such
faces are grouped by axis, and each group becomes one or more revolved features:

- local group (its faces leave part of every slab free): the prism sections are taken at a free
  height, and the difference between the prism and the part in a half plane through the axis is
  revolved as cuts (chamfers, countersinks) and bosses;
- disk group (its faces span a whole slab, as a formed dish): every section loses the disk of
  radius R about the axis, R being the largest radius of the group's faces, and the part inside
  that cylinder is one revolved boss.

The half plane is the one, of the four along the model axes, whose profile agrees best with the
other three (a gap, a slot or a cross hole changes the profile at one of them).
"""
import numpy as np

from .prism25 import _wall_or_cap

AX = "xyz"
ANG = np.sin(np.radians(5.0))       # coarse rings: a centroid lies off its facet's meridian


def _unit(k):
    e = np.zeros(3)
    e[k] = 1.0
    return e


def _residual(n, p, a3, e):
    """Sine of the angle between each normal and the meridian plane through its point, and the
    distance of the point from the axis (a3, e)."""
    q = p - a3
    q = q - np.outer(q @ e, e)
    rr = np.linalg.norm(q, axis=1)
    t = np.cross(e, q) / np.maximum(rr, 1e-30)[:, None]
    return np.abs((n * t).sum(1)), rr


def _fit_axis(n, p, w, e, i, j):
    W = np.cross(n, e)                              # n . (e x (p - a)) = (p - a) . (n x e)
    A = np.c_[W[:, i], W[:, j]] * np.sqrt(w)[:, None]
    b = (p * W).sum(1) * np.sqrt(w)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    a3 = np.zeros(3)
    a3[i], a3[j] = sol
    return a3


def _on_rings(mesh, comp, a3, e, tol):
    """Faces of comp whose vertices all lie on rings: at least 3 vertices of the patch share
    their (radius, height) about the axis. A tessellated surface of revolution is made of such
    rings; a cross hole's facets meet the test by their normals at most by a mirror image."""
    from scipy.spatial import cKDTree
    F = mesh.faces[comp]
    vid = np.unique(F)
    q = mesh.vertices[vid] - a3
    rh = np.c_[np.linalg.norm(q - np.outer(q @ e, e), axis=1), q @ e]
    cnt = np.array([len(x) for x in cKDTree(rh).query_ball_point(rh, tol / 4)])
    ok = np.zeros(int(vid.max()) + 1, bool)
    ok[vid] = cnt >= 3
    return ok[F].all(axis=1)


def _spans(H, tol):
    """Union of the height intervals of the facets (rows of H: vertex heights), merged across
    gaps below tol."""
    iv = sorted(zip(H.min(axis=1), H.max(axis=1)))
    out = [list(iv[0])]
    for a, b in iv[1:]:
        if a <= out[-1][1] + tol:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(float(a), float(b)) for a, b in out]


def _coverage(p, a3, i, j):
    """Degrees of angle about the axis covered by the points (360 less the widest gap)."""
    t = np.sort(np.degrees(np.arctan2(p[:, j] - a3[j], p[:, i] - a3[i])))
    gaps = np.diff(np.r_[t, t[0] + 360.0])
    return 360.0 - float(gaps.max())


def find_groups(mesh, k, tol, log=print):
    """Groups of off-axis faces that are surfaces of revolution about one axis parallel to k.
    Returns (groups, member): groups as dicts (a3: axis point at height 0, faces, rmin, rmax,
    ext: (kmin, kmax)), member a face mask."""
    import trimesh
    e = _unit(k)
    i, j = [m for m in range(3) if m != k]
    off = ~_wall_or_cap(mesh, k)
    member = np.zeros(len(mesh.faces), bool)
    if not off.any():
        return [], member
    N, C, A = mesh.face_normals, mesh.triangles_center, mesh.area_faces
    idx = np.nonzero(off)[0]
    adj = mesh.face_adjacency
    keep = off[adj[:, 0]] & off[adj[:, 1]]
    comps = trimesh.graph.connected_components(adj[keep], nodes=idx, min_len=1)
    found = []                                      # (axis point, faces) per patch
    for comp in sorted(comps, key=lambda c: -A[c].sum()):
        comp = np.asarray(comp)
        n, p, w = N[comp], C[comp], A[comp]
        sel = np.ones(len(comp), bool)
        a3, claimed = None, comp[:0]
        for _ in range(8):
            if sel.sum() < 6:
                break
            a3 = _fit_axis(n[sel], p[sel], w[sel], e, i, j)
            res, rr = _residual(n, p, a3, e)
            new = (res < ANG) & (rr > 2 * tol)
            if (new == sel).all():
                break
            sel = new
        if a3 is not None and sel.sum() >= 6:
            sel &= _on_rings(mesh, comp, a3, e, tol)
            # facets with one radius or one height are walls or caps of the prism, tilted by
            # poor triangle shapes where a cross hole pierces them: claimed (not cross-hole
            # facets), but not part of the revolved surface
            T = mesh.triangles[comp] - a3
            rT = np.linalg.norm(T - (T @ e)[..., None] * e, axis=2)
            hT = T @ e
            flat = (np.ptp(rT, axis=1) <= tol / 4) | (np.ptp(hT, axis=1) <= tol / 4)
            claimed = comp[sel & flat]
            sel &= ~flat
        if a3 is None or sel.sum() < 6 or w[sel].sum() < 0.3 * w.sum():
            continue
        if _coverage(p[sel], a3, i, j) < 45.0:     # a short cross cylinder seen from afar
            continue
        member[claimed] = True
        for g in found:
            if np.linalg.norm(a3 - g[0]) < tol:
                g[1].extend(comp[sel])
                break
        else:
            found.append((a3, list(comp[sel])))
    # facets left over at one radius about a group's axis or at one height are tilted walls and
    # caps of the prism (badly shaped triangles where a cross hole meets them)
    rest = np.nonzero(off & ~member)[0]
    if len(rest):
        T = mesh.triangles[rest]
        flat = np.ptp(T[..., k], axis=1) <= tol / 4
        for a3, _ in found:
            q = T - a3
            rT = np.linalg.norm(q - (q @ e)[..., None] * e, axis=2)
            flat |= np.ptp(rT, axis=1) <= tol / 4
        member[rest[flat]] = True
    from .snap2d import nice
    groups = []
    for a3, f in found:
        f = np.asarray(f)
        a3 = a3.copy()
        a3[i], a3[j] = nice(float(a3[i]), tol / 4), nice(float(a3[j]), tol / 4)
        member[f] = True
        V = mesh.vertices[np.unique(mesh.faces[f])]
        q = V - a3
        r = np.linalg.norm(q - np.outer(q @ e, e), axis=1)
        groups.append(dict(a3=a3, faces=f, rmin=float(r.min()), rmax=float(r.max()),
                           ext=(float(V[:, k].min()), float(V[:, k].max())),
                           spans=_spans(mesh.triangles[f][..., k], tol)))
    return groups, member


# ----------------------------------------------------------------------------- slabs and sections
def cap_skip(mesh, k, tol, groups):
    """Cap facets that belong to a group's surface (nearly flat facets of a dish or a round),
    inside its annulus and strictly between its end heights; they are not prism levels."""
    e = _unit(k)
    cap = np.abs(mesh.face_normals[:, k]) > np.cos(np.radians(1.0))
    out = np.zeros(len(mesh.faces), bool)
    T = mesh.triangles
    for g in groups:
        q = T - g["a3"]
        r = np.linalg.norm(q - (q @ e)[..., None] * e, axis=2)
        h = T[..., k].mean(axis=1)
        out |= (cap & (r > g["rmin"] - tol).all(1) & (r < g["rmax"] + tol).all(1)
                & (h > g["ext"][0] + tol) & (h < g["ext"][1] - tol))
    return np.nonzero(out)[0]


def set_modes(groups, slabs, tol):
    """A group whose faces span a whole slab is a disk group; the others are local."""
    from .snap2d import nice
    for g in groups:
        g["mode"] = "disk" if any(lo <= a + tol and hi >= b - tol for a, b in slabs
                                  for lo, hi in g["spans"]) else "local"
        g["R"] = nice(g["rmax"], tol / 4)
    return groups


def disk_sections(sec, groups, k, tol):
    """Section loops less the disks of the disk groups: loops inside a disk are dropped, and the
    disk's circle is added as a hole where it lies in material. None if a loop crosses a circle."""
    from matplotlib.path import Path
    i, j = [m for m in range(3) if m != k]
    out = list(sec)
    for g in groups:
        if g["mode"] != "disk":
            continue
        c2, R = np.array([g["a3"][i], g["a3"][j]]), g["R"]
        keep = []
        for pts, poly in out:
            d = np.linalg.norm(poly - c2, axis=1)
            if d.max() < R + tol:
                continue
            if d.min() < R - tol:
                return None
            keep.append((pts, poly))
            import shapely                       # nearest outline beyond R (its edges, not vertices)
            g["clear"] = min(g.get("clear", np.inf),
                             float(shapely.LinearRing(poly).distance(shapely.Point(*c2))))
        probe = c2 + [R + 2 * tol, 0.0]
        depth = sum(Path(poly).contains_point(probe) for _, poly in keep)
        if depth % 2:                            # the cylinder's wall: an exact circle (see exact_circles)
            t = np.linspace(0, 2 * np.pi, 180, endpoint=False)
            P = c2 + R * np.c_[np.cos(t), np.sin(t)]
            g.setdefault("circles", []).append(P)
            keep.append((P, P))
        out = keep
    return out


def overlap_radius(g, tol):
    """Outer radius of a disk group's boss profile: past the prism's opening (radius R) into the
    material around it, by half the clearance to the next outline (at most R / 10), so that the
    boss overlaps the prism instead of touching it on one cylinder (a merge that SolidWorks lost on
    a rebuild). R itself when there is no room for more than 4 tol."""
    from .snap2d import nice
    w = min(0.5 * (g.get("clear", g["R"]) - g["R"]), 0.1 * g["R"])
    return g["R"] if w < 4 * tol else nice(g["R"] + w, w / 4)


def exact_circles(loops, tracks, groups):
    """Loops fitted to a disk group's circle are replaced by the exact circle, so that the prism's
    opening and the revolved boss meet on one cylinder."""
    from .prism25 import Loop
    for g in groups:
        for P in g.get("circles", []):
            for m, t in enumerate(tracks):
                if t["P"] is P:
                    c2 = P.mean(axis=0)
                    loops[m] = Loop([("circle", c2, g["R"])], np.pi * g["R"] ** 2)
    return loops


# ----------------------------------------------------------------------------- profiles
def _material(polys):
    """Even-odd union of closed rings."""
    import shapely
    geo = shapely.Polygon()
    for P in polys:
        if len(P) >= 3:
            geo = geo.symmetric_difference(shapely.Polygon(P).buffer(0))
    return geo


def _half_section(mesh, k, a3, e1):
    """The part's section in the plane through the axis (a3, k) and e1, as a shapely geometry in
    (u, v) = (distance along e1, height)."""
    e = _unit(k)
    sec = mesh.section(plane_origin=a3, plane_normal=np.cross(e, e1))
    if sec is None:
        return None
    rings = [np.c_[(P - a3) @ e1, P[:, k]] for P in sec.discrete]
    return _material(rings)


def _clean(geo, tol):
    """Polygons of geo without slivers and strips thinner than tol (the mesh section's chords
    against the exact prism)."""
    w = tol / 2                                  # opening: strips thinner than tol go, corners stay
    geo = geo.buffer(-w, join_style="mitre").buffer(w, join_style="mitre")
    out = []
    for p in getattr(geo, "geoms", [geo]):
        if p.geom_type != "Polygon" or p.is_empty:
            continue
        if 2 * p.area > tol * p.length:
            out.append(p)
    return out


def _prism_half(prism, a3, e1, k, r0, r1, lo, hi):
    """Section of the prism solid by the half-plane rectangle [r0, r1] x [lo, hi] through the axis,
    as a shapely geometry in (u, v); its edges are lines (walls and caps cut by the plane)."""
    import cadquery as cq
    import shapely
    V = lambda u, v: cq.Vector(*to3((u, v), a3, e1, k))
    face = cq.Face.makeFromWires(cq.Wire.makePolygon([V(r0, lo), V(r1, lo), V(r1, hi), V(r0, hi)], close=True))
    uv = lambda p: (float((np.array(p.toTuple()) - a3) @ e1), float(p.toTuple()[k]))

    def ring(w):
        segs = [shapely.LineString([uv(e.startPoint()), uv(e.endPoint())]) for e in w.Edges()]
        return shapely.union_all(list(shapely.polygonize(segs).geoms))
    out = []
    for f in prism.intersect(face).Faces():
        p = ring(f.outerWire())
        for w in f.innerWires():
            p = p.difference(ring(w))
        out.append(p)
    return shapely.union_all(out) if out else shapely.Polygon()


def profile(mesh, k, tol, g, slabs, secs, prism=None):
    """Revolved regions of group g: (e1, adds, cuts), adds and cuts as shapely polygons in (u, v).
    secs are the prism's section loops per slab (after the disks are removed); prism, if given, is
    the extruded solid, whose half section is exact (the section loops' points are chords)."""
    import shapely
    i, j = [m for m in range(3) if m != k]
    lo_k = min(a for a, _ in slabs) - 1.0
    hi_k = max(b for _, b in slabs) + 1.0
    lo_k, hi_k = min(lo_k, g["ext"][0] - 1.0), max(hi_k, g["ext"][1] + 1.0)
    cands = []
    for m in (i, j):
        for s in (1.0, -1.0):
            e1 = s * _unit(m)
            S = _half_section(mesh, k, g["a3"], e1)
            if S is None:
                continue
            if g["mode"] == "disk":
                band = shapely.box(0.0, lo_k, overlap_radius(g, tol), hi_k)
                adds, cuts = _clean(S.intersection(band), tol), []
            else:
                r0, r1 = max(0.0, g["rmin"] - 2 * tol), g["rmax"] + 2 * tol
                band = shapely.box(r0, lo_k, r1, hi_k)
                c2 = np.array([g["a3"][i], g["a3"][j]])
                d2 = np.array([e1[i], e1[j]])
                if prism is not None:
                    Ph = _prism_half(prism, g["a3"], e1, k, r0, r1, lo_k, hi_k)
                else:
                    rects = []
                    for (a, b), sec in zip(slabs, secs):
                        M = _material([poly for _, poly in sec])
                        ray = shapely.LineString([c2 + r0 * d2, c2 + r1 * d2])
                        hit = M.intersection(ray)
                        for seg in getattr(hit, "geoms", [hit]):
                            if seg.geom_type != "LineString" or seg.is_empty:
                                continue
                            u = [(np.asarray(q) - c2) @ d2 for q in seg.coords]
                            rects.append(shapely.box(min(u), a, max(u), b))
                    Ph = shapely.union_all(rects) if rects else shapely.Polygon()
                Sb = S.intersection(band)
                cuts = _clean(Ph.difference(Sb), tol)
                adds = _clean(Sb.difference(Ph), tol)
            cands.append((e1, adds, cuts))
    if not cands:
        return None
    def geo(c):
        return shapely.union_all(c[1] + c[2]) if c[1] or c[2] else shapely.Polygon()
    G = [geo(c) for c in cands]
    score = [sum(G[a].symmetric_difference(G[b]).area for b in range(len(G)) if b != a) for a in range(len(G))]
    return cands[int(np.argmin(score))]


def _exact_ring(P, UV, tol):
    """The ring P with its points on curved surfaces replaced by the mesh vertices UV (the vertex
    rings of the group's faces as (radius, height)): points where the section plane cuts a
    facet edge lie inside the surface by the chord sag. Ring points on straight axis-aligned
    pieces (prism faces, the band's cuts) are exact and kept."""
    from .revolve import _exact_points
    if UV is None or not len(UV):
        return P
    n = len(P)
    keep = []
    for m in range(n):
        a, b, c = P[m - 1], P[m], P[(m + 1) % n]
        aligned = lambda p, q: abs(p[0] - q[0]) < 1e-6 or abs(p[1] - q[1]) < 1e-6
        if aligned(a, b) and aligned(b, c):
            keep.append(b)
    cand = np.vstack([UV] + ([np.array(keep)] if keep else []))
    X = _exact_points(np.vstack([P, P[:1]]), cand, tol)
    if len(X) < 3:
        return P
    # a straight axis-aligned piece (a prism face, the band's cut) has points only at its ends;
    # points along it keep the fitter from folding it into a neighbouring arc
    out = []
    for m in range(len(X)):
        a, b = X[m], X[(m + 1) % len(X)]
        out.append(a)
        L = float(np.linalg.norm(b - a))
        if (abs(a[0] - b[0]) < 1e-9 or abs(a[1] - b[1]) < 1e-9) and L > 4 * tol:
            k = min(int(L / tol), 8)
            out += [a + (b - a) * t for t in np.arange(1, k) / k]
    return np.array(out)


def fit_region(poly, tol, UV=None):
    """Fitted loops of one profile polygon: [outer, holes...], with the point sets they fit."""
    from .prism25 import build_loop
    rings = [np.asarray(poly.exterior.coords)[:-1]] + [np.asarray(r.coords)[:-1] for r in poly.interiors]
    rings = [_exact_ring(P, UV, tol) for P in rings]
    from .snap2d import _dist_to_loop
    loops = []
    for P in rings:                              # the vertex rings are exact: fit them closely
        L = build_loop(P, tol / 4)
        if L is None or _dist_to_loop(P, L.edges) > tol / 2:
            L = build_loop(P, tol)
        if L is None:
            return None
        loops.append(L)
    return loops, rings


def to3(q, a3, e1, k):
    p = a3 + float(q[0]) * e1
    p[k] = float(q[1])
    return p


def region_solid(loops, a3, e1, k):
    """Revolve one profile region (outer loop, then holes) 360 degrees about the axis."""
    import cadquery as cq

    def V(q):
        return cq.Vector(*to3(q, a3, e1, k))
    wires = []
    for L in loops:
        es = []
        for e in L.edges:
            if e[0] == "line":
                es.append(cq.Edge.makeLine(V(e[1]), V(e[2])))
            elif e[0] == "arc":
                es.append(cq.Edge.makeThreePointArc(V(e[1]), V(e[2]), V(e[3])))
            elif e[0] == "circle":
                c3 = to3(e[1], a3, e1, k)
                es.append(cq.Edge.makeCircle(e[2], cq.Vector(*c3), cq.Vector(*np.cross(_unit(k), e1))))
            else:
                return None
        wires.append(cq.Wire.assembleEdges(es))
    face = cq.Face.makeFromWires(wires[0], wires[1:])
    ax0 = a3.copy()
    ax1 = a3 + _unit(k)
    return cq.Solid.revolve(face, 360.0, cq.Vector(*ax0), cq.Vector(*ax1))


def features(mesh, k, tol, groups, slabs, secs, snap=True, log=print, prism=None):
    """Fitted revolved features of all groups: list of dict(op, a3, e1, k, loops, mode, R, solid)
    in build order (disk bosses, local bosses, local cuts); None if a group cannot be fitted."""
    from .snap2d import snap_loops
    out = []
    for gi, g in enumerate(groups):
        pr = profile(mesh, k, tol, g, slabs, secs, prism)
        if pr is None:
            log(f"  revolved group {gi + 1}: no section through its axis")
            return None
        e1, adds, cuts = pr
        e = _unit(k)
        V = mesh.vertices[np.unique(mesh.faces[g["faces"]])] - g["a3"]
        UV = np.unique(np.round(np.c_[np.linalg.norm(V - np.outer(V @ e, e), axis=1), V @ e], 6), axis=0)
        n_feat = 0
        for op, polys in (("add", adds), ("cut", cuts)):
            for poly in polys:
                fr = fit_region(poly, tol, UV)
                if fr is None:
                    return None
                loops, rings = fr
                if snap:
                    loops = snap_loops(loops, rings, tol, log=lambda *a: None)
                from .prism25 import exact_arcs
                loops = [exact_arcs(L) for L in loops]
                solid = region_solid(loops, g["a3"], e1, k)
                if solid is None or not solid.isValid():
                    log(f"  revolved group {gi + 1}: invalid revolved {op}")
                    return None
                out.append(dict(op=op, a3=g["a3"], e1=e1, k=k, loops=loops, mode=g["mode"], R=g["R"],
                                solid=solid, group=gi))
                n_feat += 1
        log(f"  revolved group {gi + 1} ({g['mode']}, axis through ({g['a3'][0]:.3f}, {g['a3'][1]:.3f}, "
            f"{g['a3'][2]:.3f}), r {g['rmin']:.3f}..{g['rmax']:.3f}, half plane {'+-'[int(e1.sum() < 0)]}"
            f"{AX[int(np.argmax(np.abs(e1)))]}): {len(adds)} boss(es), {len(cuts)} cut(s), "
            f"{sum(len(L.edges) for f in out[len(out) - n_feat:] for L in f['loops'])} profile edges")
    rank = {("disk", "add"): 0, ("local", "add"): 1, ("local", "cut"): 2, ("disk", "cut"): 2}
    return sorted(out, key=lambda f: rank[(f["mode"], f["op"])])
