"""Straight cuts across a 2.5D part along a second model axis (a chamfer, a slanted notch).

The stacked prism is built from sections that ignore the slanted planar faces, so it holds the
material those faces cut away. Each slanted plane whose normal is perpendicular to a model axis
j bounds a cut that is a prism along j: its outline is the difference between the prism and the
part in a plane across j, fitted with lines and arcs, and it runs over the plane's extent along
j (through the part where the plane reaches the part's ends).
"""
import numpy as np


def plane_frame(j):
    """In-plane axes (w, h) across axis j, with w x h = e_j."""
    E = np.eye(3)
    return E[(j + 1) % 3], E[(j + 2) % 3]


def evenodd(loops):
    """Section loops (2D point arrays) as one shapely geometry by even-odd nesting."""
    import shapely
    geo = shapely.Polygon()
    for P in loops:
        if len(P) >= 3:
            geo = geo.symmetric_difference(shapely.Polygon(P).buffer(0))
    return geo


def loops_of(geo):
    """Boundary loops (2D point arrays) of a shapely geometry."""
    out = []
    for p in getattr(geo, "geoms", [geo]):
        if p.is_empty or p.geom_type != "Polygon":
            continue
        out.append(np.asarray(p.exterior.coords)[:-1])
        out += [np.asarray(r.coords)[:-1] for r in p.interiors]
    return out


def shadow(sections, tol):
    """Union of several sections of one slab (each a list of (exact points, polyline) loops):
    the prism that holds every one of them. Each boundary loop keeps the exact points of the
    sections that lie on it, in order along it."""
    import shapely
    geo = shapely.Polygon()
    for sec in sections:
        geo = geo.union(evenodd([poly for _, poly in sec]))
    src = np.vstack([pts for sec in sections for pts, _ in sec])
    out = []
    for R in loops_of(geo):
        ring = shapely.LinearRing(R)
        near = shapely.distance(shapely.points(src), ring) < 0.25 * tol
        Q = src[near]
        if len(Q) < 3:
            out.append((R, R))
            continue
        # ring corners that no section point is close to (where a hole's fill meets the outline)
        far = shapely.distance(shapely.points(R), shapely.MultiPoint(Q)) > 0.25 * tol
        Q = np.vstack([Q, R[far]])
        t = shapely.line_locate_point(ring, shapely.points(Q))
        o = np.argsort(t, kind="stable")
        Q, t = Q[o], t[o]
        keep = np.r_[True, np.diff(t) > 1e-9]
        out.append((Q[keep], R))
    return out


def brep_section(solid, org, E, w, h, tol):
    """Section of a B-rep solid by the plane through org normal to E, as a shapely geometry in
    (w, h) coordinates (exact faces, discretized finely; no tessellation of the solid)."""
    import cadquery as cq
    import shapely
    R = 4.0 * solid.BoundingBox().DiagonalLength
    pl = cq.Face.makePlane(R, R, cq.Vector(*org), cq.Vector(*E))
    sec = solid.intersect(pl)
    geo = shapely.Polygon()

    for f in sec.Faces():
        poly = shapely.Polygon()
        for k_, wi in enumerate([f.outerWire()] + f.innerWires()):
            pts = []
            for e in _ordered_edges(wi):                 # in the order and direction traversed
                pts += e
            Q = np.array(pts) - org
            ring2 = shapely.Polygon(np.c_[Q @ w, Q @ h]).buffer(0)
            poly = ring2 if k_ == 0 else poly.difference(ring2)
        geo = geo.union(poly)
    return geo


def _ordered_edges(wire, n_arc=64):
    """Points of a wire's edges in the order and direction they are traversed."""
    from OCP.BRepTools import BRepTools_WireExplorer
    import cadquery as cq
    ex = BRepTools_WireExplorer(wire.wrapped)
    out = []
    while ex.More():
        e = cq.Edge(ex.Current())
        n = 1 if e.geomType() == "LINE" else n_arc
        ts = np.linspace(0, 1, n + 1)[:-1]
        P = [e.positionAt(t).toTuple() for t in ts]
        from OCP.TopAbs import TopAbs_REVERSED
        if ex.Current().Orientation() == TopAbs_REVERSED:
            P = [e.positionAt(t).toTuple() for t in (1 - np.linspace(0, 1, n + 1))[:-1]]
        out.append(P)
        ex.Next()
    return out


def cross_cuts(mesh, solid, planes, tol, log=print):
    """Cuts along second axes for the slanted planes; returns (solid, cuts) or None when a plane
    has no cut region or lies at a compound angle. A cut is dict(j, lo, hi, loops) with loops in
    the (w, h) coordinates of plane_frame(j)."""
    import cadquery as cq
    from .core import shape_to_mesh
    from .prism25 import build_loop, Loop
    from .revfeat import _clean
    from .revolve import _plane_section
    from .snap2d import snap_loops, nice
    if not planes:
        return solid, []
    if any(p["j"] is None for p in planes):
        log("  2.5D: a slanted face at a compound angle")
        return None
    cuts = []
    covered = [False] * len(planes)
    lo_all, hi_all = mesh.bounds
    for n_, pl in enumerate(planes):
        if covered[n_]:
            continue
        j = pl["j"]
        w, h = plane_frame(j)
        E = np.eye(3)[j]
        V = mesh.vertices[np.unique(mesh.faces[pl["faces"]])]
        org = V.mean(axis=0)
        # off the middle: a symmetric part's tessellation has seams in its plane of symmetry
        org[j] = V[:, j].min() + 0.382 * np.ptp(V[:, j])
        S = _plane_section(mesh, org, E, w, h)
        P = brep_section(solid, org, E, w, h, tol)
        region = _clean(P.difference(S).buffer(-2 * tol, join_style="mitre").buffer(2 * tol, join_style="mitre"), tol)
        import shapely
        seg = shapely.MultiPoint(np.c_[(V - org) @ w, (V - org) @ h])
        region = [r for r in region if r.distance(seg) < 2 * tol and 2 * r.area > 4 * tol * r.length]
        if not region:                           # a tilted facet of a cap or wall, not a cut:
            log(f"  2.5D: no cut region for the slanted face with normal {np.round(pl['n'], 3)}; "
                "left to the distance check")
            covered[n_] = True
            continue
        # every slanted face across j that bounds this region belongs to the same cut
        a, b = float(V[:, j].min()), float(V[:, j].max())
        for m_, q in enumerate(planes):
            if q["j"] != j or covered[m_]:
                continue
            Vq = mesh.vertices[np.unique(mesh.faces[q["faces"]])]
            sq = shapely.MultiPoint(np.c_[(Vq - org) @ w, (Vq - org) @ h])
            if any(r.distance(sq) < 2 * tol for r in region):
                covered[m_] = True
                a, b = min(a, float(Vq[:, j].min())), max(b, float(Vq[:, j].max()))
        # outline in the plane through the origin: offset by the section position along j
        off = np.array([org @ w, org @ h])
        loops, pts = [], []
        for r in region:
            Q = np.asarray(r.exterior.coords)[:-1] + off
            L = build_loop(Q, tol)
            if L is None:
                return None
            loops.append(L)
            pts.append(Q)
        loops = snap_loops(loops, pts, tol, log=lambda *a_: None)
        for m_, L in enumerate(loops):                 # corners of line outlines to clean values
            if all(e[0] == "line" for e in L.edges):
                r_ = lambda q: np.array([nice(float(q[0]), 1e-5), nice(float(q[1]), 1e-5)]) + 0.0
                loops[m_] = Loop([("line", r_(e[1]), r_(e[2])) for e in L.edges], L.area)
        a, b = nice(a, tol / 4), nice(b, tol / 4)
        open_lo = a <= lo_all[j] + tol
        open_hi = b >= hi_all[j] - tol
        reach = 0.05 * float(np.linalg.norm(mesh.extents)) + 1.0
        z0 = a - (reach if open_lo else 0.0)
        z1 = b + (reach if open_hi else 0.0)
        cut = dict(j=j, lo=z0, hi=z1, a=a, b=b, open_lo=bool(open_lo), open_hi=bool(open_hi), loops=loops)
        solid = solid.cut(cut_solid(cut)).clean()
        cuts.append(cut)
        log(f"  2.5D: cut along {'xyz'[j]} over [{a:.4f}, {b:.4f}]"
            f"{' (through)' if open_lo and open_hi else ''}, {sum(len(L.edges) for L in loops)} edges")
    return solid, cuts


def _wire(L, j, z):
    import cadquery as cq
    w, h = plane_frame(j)
    E = np.eye(3)[j]
    V = lambda q: cq.Vector(*(float(q[0]) * w + float(q[1]) * h + z * E))
    es = []
    for e in L.edges:
        if e[0] == "line":
            es.append(cq.Edge.makeLine(V(e[1]), V(e[2])))
        elif e[0] == "arc":
            es.append(cq.Edge.makeThreePointArc(V(e[1]), V(e[2]), V(e[3])))
        else:
            es.append(cq.Edge.makeCircle(e[2], V(e[1]), cq.Vector(*E)))
    return cq.Wire.assembleEdges(es)


def cut_solid(cut):
    import cadquery as cq
    j = cut["j"]
    body = None
    for L in cut["loops"]:
        f = cq.Face.makeFromWires(_wire(L, j, cut["lo"]), [])
        D = [0.0, 0.0, 0.0]
        D[j] = cut["hi"] - cut["lo"]
        piece = cq.Solid.extrudeLinear(f, cq.Vector(*D))
        body = piece if body is None else body.fuse(piece)
    return body
