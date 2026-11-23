"""Bodies of revolution: a profile in a plane through the axis, revolved 360 degrees.

A mesh whose two moments of inertia about axes perpendicular to one candidate axis are equal (a
turned part: pawn, shaft, knob) is cut by the plane through the axis and one of its own mesh
vertices, so the section passes through the vertex rings of the tessellation. The section loop is
mirror-symmetric about the axis; one half is fitted with lines and arcs (prism25's mirrored chain
fitting, cut exactly at the axis), closed along the axis, snapped, revolved, and checked against
the mesh both ways like a 2.5D model. Loops that do not reach the axis (grooves, inner voids) are
fitted on their own and revolved as holes in the profile.

The profile lives in (a, b) coordinates of the model plane through the axis spanned by two model
axes, so the CadQuery and SolidWorks generators can draw it like a 2.5D sketch.
"""
import time

import numpy as np

from .core import Result

AX = "xyz"


def _axis_candidates(mesh):
    """Axis directions along which the mesh may be a body of revolution: model axes whose two
    perpendicular principal moments agree within 1 %, through the centre of mass."""
    c = mesh.center_mass
    I = mesh.moment_inertia
    out = []
    for k in range(3):
        i, j = [m for m in range(3) if m != k]
        off = max(abs(I[i, j]), abs(I[i, k]), abs(I[j, k]))
        if abs(I[i, i] - I[j, j]) < 0.01 * max(I[i, i], I[j, j]) and off < 0.01 * max(I[i, i], I[j, j]):
            d = np.zeros(3)
            d[k] = 1.0
            out.append((k, c - (c @ d) * d, d))       # the axis point at level 0 along the axis
    return out


def _section(mesh, c, d, e1):
    """Closed section loops of the mesh in the plane through the axis (c, d) spanned by d and e1,
    as (u, v) points: u along e1, v along d (relative to c)."""
    n = np.cross(d, e1)
    sec = mesh.section(plane_origin=c, plane_normal=n)
    if sec is None:
        return []
    loops = []
    for P3 in sec.discrete:
        Q = np.column_stack([(P3 - c) @ e1, (P3 - c) @ d])
        if len(Q) > 3 and np.linalg.norm(Q[0] - Q[-1]) < 1e-9:
            Q = Q[:-1]
        # drop repeated points
        keep = [0]
        for m in range(1, len(Q)):
            if np.linalg.norm(Q[m] - Q[keep[-1]]) > 1e-9:
                keep.append(m)
        Q = Q[keep]
        if len(Q) >= 3:
            loops.append(Q)
    return loops


def _exact_points(G, UV, tol, near=None):
    """The mesh vertices (UV: (radius, height) of every vertex, all on the profile curve) near the
    section polyline G, ordered by where they project onto it. The section itself mixes vertices
    with points where the plane cuts facet edges, which lie on chords inside the surface; the
    vertices of all the rings are exact and dense."""
    seg_a, seg_b = G[:-1], G[1:]
    d = seg_b - seg_a
    L2 = np.maximum((d * d).sum(1), 1e-30)
    cum = np.r_[0.0, np.cumsum(np.sqrt(L2))]
    best_s = np.full(len(UV), np.inf)
    best_d = np.full(len(UV), np.inf)
    for m in range(len(d)):                              # nearest segment of the polyline
        t = np.clip(((UV - seg_a[m]) @ d[m]) / L2[m], 0.0, 1.0)
        q = seg_a[m] + t[:, None] * d[m]
        dist = np.linalg.norm(UV - q, axis=1)
        better = dist < best_d
        best_d[better] = dist[better]
        best_s[better] = cum[m] + t[better] * np.sqrt(L2[m])
    # chords of the section lie inside the surface by their sag; vertices lie on it
    near = best_d < (near if near is not None else max(0.02 * float(np.ptp(G, axis=0).max()), 10 * tol))
    P, sp = UV[near], best_s[near]
    o = np.argsort(sp, kind="stable")
    P = P[o]
    keep = [0]                       # one point per vertex ring (a ring's vertices agree to
    for m in range(1, len(P)):       # rounding, and would jitter the order)
        if np.linalg.norm(P[m] - P[keep[-1]]) > 0.05 * tol:
            keep.append(m)
    return P[keep]


def _profile(loops, tol, log, UV=None):
    """Profile loops (u >= 0 half plane) as prism25 Loops in (u, v): the loop crossing the axis is
    fitted on its u >= 0 half, cut at the axis and closed along it; loops entirely on u > 0 are
    fitted whole; loops on u < 0 are mirror images and dropped."""
    from .prism25 import Loop, _half_chain, _chain_edges, build_loop, _area
    out, pts = [], []
    for Q in loops:
        if Q[:, 0].min() > tol:                       # off the axis: a groove or void ring
            L = build_loop(Q, tol)
            if L is None:
                return None
            out.append(L)
            pts.append(Q)
            continue
        if Q[:, 0].max() < -tol:
            continue
        # turn (u, v) into the convention of the mirrored fitting (axis = line a = 0)
        H = _half_chain(Q, 0)
        if H is None:
            log("  revolve: the section does not cross the axis twice")
            return None
        if UV is not None:
            G = np.vstack([[0.0, H[0][1]], H, [0.0, H[-1][1]]]) if False else H
            X = _exact_points(G, UV[UV[:, 0] >= 0], tol)
            if len(X) >= 3:
                H = X
        E = None
        if UV is not None:                            # vertex rings are exact: fit them closely
            from .snap2d import _dist_to_loop
            E = _chain_edges(H, 0, 0, tol / 4)
            if E and _dist_to_loop(H, [e for e in E]) > tol / 2:
                E = None
        if not E:
            E = _chain_edges(H, 0, 0, tol)
        if not E:
            return None
        a, b = np.asarray(E[-1][-1], float), np.asarray(E[0][1], float)
        if np.linalg.norm(a - b) > 1e-9:              # close along the axis
            E = list(E) + [("line", a, b)]
        Ph = np.vstack([H, [[0.0, H[-1][1]], [0.0, H[0][1]]]])
        poly = np.array([np.asarray(e[1], float) for e in E])
        out.append(Loop([(e[0],) + tuple(np.asarray(q, float) for q in e[1:]) for e in E], _area(poly)))
        pts.append(Ph)
    return out, pts


def _solid(loops, c, d, e1):
    """Revolve the profile loops (u, v) about the axis (c, d): the outer loop(s) less the inner."""
    import cadquery as cq
    from .prism25 import _wire

    def to3(q):
        q = np.asarray(q, float)
        return c + q[0] * e1 + q[1] * d
    wires = []
    for L in loops:
        es = []
        for e in L.edges:
            if e[0] == "line":
                es.append(cq.Edge.makeLine(cq.Vector(*to3(e[1])), cq.Vector(*to3(e[2]))))
            elif e[0] == "arc":
                es.append(cq.Edge.makeThreePointArc(*(cq.Vector(*to3(q)) for q in e[1:])))
            else:
                return None
        wires.append(cq.Wire.assembleEdges(es))
    order = sorted(range(len(loops)), key=lambda m: -abs(loops[m].area))
    outer = wires[order[0]]
    inner = [wires[m] for m in order[1:]]
    face = cq.Face.makeFromWires(outer, inner)
    return cq.Solid.revolve(face, 360.0, cq.Vector(*c), cq.Vector(*(c + d)))


def reconstruct(mesh, tol, log=print):
    """Revolved solid, or None when the mesh is not a body of revolution or the fit fails."""
    from .prism25 import _distance_check
    from .snap2d import snap_loops
    t0 = time.time()
    cands = _axis_candidates(mesh)
    if not cands:
        log("  revolve: no axis with equal perpendicular moments")
        return None
    for k, c, d in cands:
        i, j = [m for m in range(3) if m != k]
        V = mesh.vertices - c
        r = np.linalg.norm(V - np.outer(V @ d, d), axis=1)
        # the meridian through a vertex of the widest ring (that ring's vertices lie on the plane);
        # the model axis closest to it spans the sketch plane with the axis
        v0 = V[int(np.argmax(r))]
        e1 = v0 - (v0 @ d) * d
        e1 /= np.linalg.norm(e1)
        # prefer a model axis direction when a vertex lies on it (a CAD tessellation usually has a
        # vertex at angle 0)
        for m in (i, j):
            for s in (1.0, -1.0):
                em = np.zeros(3)
                em[m] = s
                ang = np.abs(np.arctan2(np.cross(np.outer(np.ones(len(V)), em), V) @ d, V @ em))
                if np.any((ang < 1e-6) & (r > 0.5 * r.max())):
                    e1 = em
                    break
            else:
                continue
            break
        loops = _section(mesh, c, d, e1)
        if not loops:
            continue
        UV = np.unique(np.round(np.column_stack([r, V @ d]), 5), axis=0)
        prof = _profile(loops, tol, log, UV)
        if prof is None:
            continue
        L, pts = prof
        L = snap_loops(L, pts, tol, log=log)
        try:
            solid = _solid(L, c, d, e1)
        except Exception as ex:
            log(f"  revolve about {AX[k]}: solid failed ({type(ex).__name__})")
            continue
        if solid is None:
            continue
        valid = solid.isValid()
        vol_err = abs(solid.Volume() / mesh.volume - 1)
        dmax, dback = _distance_check(solid, mesh) if valid else (np.inf, np.inf)
        n_edges = sum(len(x.edges) for x in L)
        log(f"  revolve about {AX[k]}: {len(L)} profile loop(s), {n_edges} edges, volume error "
            f"{100 * vol_err:.3f} %, max distance mesh->model {dmax:.4f} mm, model->mesh {dback:.4f} mm, "
            f"valid={valid} ({time.time() - t0:.1f}s)")
        if not (valid and vol_err < 0.03 and dmax < 1.5 * tol and dback < 4 * tol):
            log("  revolve: checks failed")
            continue
        res = Result(shape=solid, mode="exact")
        res.extras = dict(prims=[], analytic_fraction=1.0, summary=f"revolve about {AX[k]}",
                          volume_error=vol_err, brep=None, valid=True, passed=True,
                          revolve=dict(axis=k, c=c, d=d, e1=e1, loops=L))
        return res
    return None


# ----------------------------------------------------------------------------- revolve + through cuts
def _half_area(mesh, c, d, e1):
    """Area of the part's section in the half plane u >= 0 through the axis (c, d) towards e1."""
    import shapely
    loops = _section(mesh, c, d, e1)
    geo = shapely.Polygon()
    for Q in loops:
        geo = geo.symmetric_difference(shapely.Polygon(Q).buffer(0))
    return geo.intersection(shapely.box(0, -1e6, 1e6, 1e6)).area


def _profile_solid(mesh, k, c, d, e1, tol, log, UV):
    from .snap2d import snap_loops
    loops = _section(mesh, c, d, e1)
    if not loops:
        return None
    prof = _profile(loops, tol, log, UV)
    if prof is None:
        return None
    L, pts = prof
    L = snap_loops(L, pts, tol, log=lambda *a: None)
    from .prism25 import exact_arcs
    L = [exact_arcs(x) for x in L]
    try:
        solid = _solid(L, c, d, e1)
    except Exception:
        return None
    return (L, solid) if solid is not None and solid.isValid() else None


def _cut_directions(mesh, c, d, tol):
    """Directions D (perpendicular to the axis) of straight cuts: the planar faces whose normals
    leave the meridian planes (notch walls) have normals perpendicular to D. Returns a list of
    (D, faces)."""
    from .revfeat import _residual
    N, C = mesh.face_normals, mesh.triangles_center
    res, rr = _residual(N, C, c, d)
    off = (res > np.sin(np.radians(5.0))) & (rr > 2 * tol)
    nh = N - np.outer(N @ d, d)
    ln = np.linalg.norm(nh, axis=1)
    off &= ln > 0.3
    dirs = []
    for f in np.nonzero(off)[0]:
        D = np.cross(d, nh[f] / ln[f])
        for g in dirs:
            if abs(g[0] @ D) > np.cos(np.radians(1.0)):
                g[1].append(f)
                break
        else:
            dirs.append((D, [f]))
    A = mesh.area_faces
    tot = A[off].sum()
    return [(D, np.array(f)) for D, f in dirs if A[f].sum() > 0.02 * tot]


def _plane_section(shape_or_mesh, origin, normal, w, e):
    """Section by a plane as a shapely geometry in (w, e) coordinates."""
    import shapely
    import trimesh
    if not isinstance(shape_or_mesh, trimesh.Trimesh):
        from .core import shape_to_mesh
        shape_or_mesh = shape_to_mesh(shape_or_mesh, tolerance=1e-3, angular=0.02)
    sec = shape_or_mesh.section(plane_origin=origin, plane_normal=normal)
    geo = shapely.Polygon()
    if sec is None:
        return geo
    for P in sec.discrete:
        Q = np.c_[(P - origin) @ w, (P - origin) @ e]
        if len(Q) >= 3:
            geo = geo.symmetric_difference(shapely.Polygon(Q).buffer(0))
    return geo


def reconstruct_with_cuts(mesh, tol, log=print):
    """A body of revolution with straight cuts across it (a chess rook's crown notches): the
    profile is taken in the half plane of largest section (clear of the cuts), revolved, and each
    cut, a prism along a direction D perpendicular to the axis, is the difference between the
    revolved body and the part in a plane across D near the rim (where no other cut reaches),
    extruded through the body (both sides), else from the axis outward."""
    import shapely
    import cadquery as cq
    from .prism25 import _distance_check, build_loop
    from .revfeat import _clean
    t0 = time.time()
    from .snap2d import nice
    for k, c, d in _axis_candidates(mesh):
        i, j = [m for m in range(3) if m != k]
        c = np.array([nice(float(x), tol / 4) for x in c]) + 0.0
        cuts = _cut_directions(mesh, c, d, tol)
        if not cuts:
            continue
        # half plane: the largest section (between the cuts); a model axis when it is as large
        cand = []
        for ang in np.radians(np.arange(0, 360, 5.0)):
            e1 = np.cos(ang) * np.eye(3)[i] + np.sin(ang) * np.eye(3)[j]
            cand.append((_half_area(mesh, c, d, e1), e1))
        amax = max(a for a, _ in cand)
        axis_dirs = [s * np.eye(3)[m] for m in (i, j) for s in (1.0, -1.0)]
        good = [(a, e) for a, e in cand if a > amax * (1 - 1e-4)]
        e1 = next((e for a, e in good if any(np.allclose(e, x) for x in axis_dirs)), good[0][1])
        e1 = np.round(e1, 12) + 0.0
        # exact profile points: vertices of facets on the surface of revolution
        from .revfeat import _residual
        res, _ = _residual(mesh.face_normals, mesh.triangles_center, c, d)
        V = mesh.vertices[np.unique(mesh.faces[res < np.sin(np.radians(5.0))])] - c
        UV = np.unique(np.round(np.c_[np.linalg.norm(V - np.outer(V @ d, d), axis=1), V @ d], 5), axis=0)
        fit = _profile_solid(mesh, k, c, d, e1, tol, log, UV)
        if fit is None:
            continue
        L, rev = fit
        R = float(np.linalg.norm(mesh.vertices - c - np.outer((mesh.vertices - c) @ d, d), axis=1).max())
        out, solid = [], rev
        from .core import shape_to_mesh
        rev_mesh = shape_to_mesh(rev, tolerance=tol / 4, angular=0.1)   # sections of the revolve
        for D, faces in cuts:
            # the cut direction's angle about the axis, to a whole degree when within 0.05 deg
            ang = np.degrees(np.arctan2(D @ np.eye(3)[j], D @ np.eye(3)[i]))
            if abs(ang - round(ang)) < 0.05:
                a_ = np.radians(round(ang))
                D = np.round(np.cos(a_) * np.eye(3)[i] + np.sin(a_) * np.eye(3)[j], 12) + 0.0
            w = np.cross(D, d)
            t = (mesh.triangles_center[faces] - c) @ D
            side = 1.0 if (t > 0).sum() >= (t < 0).sum() else -1.0
            tt = float(np.median(np.abs(t[t * side > 0])))
            org = c + side * tt * D
            S = _plane_section(mesh, org, D, w, d)
            Rv = _plane_section(rev_mesh, org, D, w, d)
            # strips along curved walls (the tessellated revolve against the mesh) are not cuts
            region = _clean(Rv.difference(S).buffer(-2 * tol, join_style="mitre").buffer(2 * tol, join_style="mitre"), tol)
            if region:
                big = max(p.area for p in region)
                region = [p for p in region if p.area > 0.01 * big and 2 * p.area > 4 * tol * p.length]
            if not region:
                log(f"  revolve+cuts: no cut region across {np.round(D, 3)}")
                break
            loops, pts = [], []
            for poly in region:
                P = np.asarray(poly.exterior.coords)[:-1]
                Lc = build_loop(P, tol)
                if Lc is None:
                    break
                loops.append(Lc)
                pts.append(P)
            from .snap2d import snap_loops, nice
            loops = snap_loops(loops, pts, tol, log=lambda *a: None)
            # outlines of lines: corners to clean values within 1e-5 mm (corners on slanted lines
            # are recomputed intersections, 2.5000005 for 2.5); shared corners round alike
            from .prism25 import Loop
            for m_, Lc in enumerate(loops):
                if all(e[0] == "line" for e in Lc.edges):
                    r_ = lambda q: np.array([nice(float(q[0]), 1e-5), nice(float(q[1]), 1e-5)]) + 0.0
                    loops[m_] = Loop([("line", r_(e[1]), r_(e[2])) for e in Lc.edges], Lc.area)
            # through the body when the cut's walls lie on both sides of the axis
            both = bool((t > tol).any() and (t < -tol).any())
            body = None
            for Lc in loops:
                f = cq.Face.makeFromWires(_cut_wire(Lc, c, w, d), [])
                if both:
                    piece = cq.Solid.extrudeLinear(f.translate(cq.Vector(*(-1.2 * R * D))), cq.Vector(*(2.4 * R * D)))
                else:
                    piece = cq.Solid.extrudeLinear(f, cq.Vector(*(1.2 * R * side * D)))
                body = piece if body is None else body.fuse(piece)
            solid = solid.cut(body).clean()
            out.append(dict(D=D, w=w, loops=loops, both=both, side=side))
            log(f"  revolve+cuts: cut across ({D[0]:.3f}, {D[1]:.3f}, {D[2]:.3f}), "
                f"{sum(len(x.edges) for x in loops)} edges, {'through' if both else 'one side'}")
        else:
            from .core import mass_props
            valid = solid.isValid() and len(solid.Solids()) == 1
            vol_err = abs(mass_props(solid)[0] / mesh.volume - 1)
            dmax, dback = _distance_check(solid, mesh) if valid else (np.inf, np.inf)
            log(f"  revolve about {AX[k]} with {len(out)} cut(s): {sum(len(x.edges) for x in L)} profile edges, "
                f"volume error {100 * vol_err:.3f} %, max distance mesh->model {dmax:.4f} mm, model->mesh "
                f"{dback:.4f} mm, valid={valid} ({time.time() - t0:.1f}s)")
            if not (valid and vol_err < 0.03 and dmax < 1.5 * tol and dback < 4 * tol):
                log("  revolve+cuts: checks failed")
                continue
            r = Result(shape=solid, mode="exact")
            r.extras = dict(prims=[], analytic_fraction=1.0, summary=f"revolve about {AX[k]} with {len(out)} cut(s)",
                            volume_error=vol_err, brep=None, valid=True, passed=True,
                            revolve=dict(axis=k, c=c, d=d, e1=e1, loops=L, cuts=out))
            return r
    return None


def _cut_wire(Lc, c, w, d):
    """A cut's outline (w, h) in the plane through the axis point c spanned by w and d."""
    import cadquery as cq
    V = lambda q: cq.Vector(*(c + float(q[0]) * w + float(q[1]) * d))
    es = []
    for e in Lc.edges:
        if e[0] == "line":
            es.append(cq.Edge.makeLine(V(e[1]), V(e[2])))
        elif e[0] == "arc":
            es.append(cq.Edge.makeThreePointArc(V(e[1]), V(e[2]), V(e[3])))
        else:
            es.append(cq.Edge.makeCircle(e[2], V(e[1]), cq.Vector(*np.cross(w, d))))
    return cq.Wire.assembleEdges(es)
