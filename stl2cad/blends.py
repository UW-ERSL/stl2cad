"""Edge blends (fillets and chamfers) on a 2.5D part.

A blend is a strip of off-axis faces along an edge of the sharp part: a partial cylinder (a
fillet along a straight edge), a strip of planes and cones (a chamfer along an outline), a cone
around a cross hole (a countersink). The stacked prism is built from sections that ignore the
blend faces (taken away from them), so it has the sharp edge; the blend is then made on that
edge with the radius or distance measured from the mesh, as a designer adds a Fillet or a
Chamfer feature.

For an edge where two faces meet at the angle b (the wedge that is rounded or cut, 90 degrees at
a box corner), with legs L (the distance from the edge to where the blend meets each face):
a chamfer has distance L and its middle lies L cos(b/2) from the edge; a fillet has radius
r = L tan(b/2) and its middle lies r / sin(b/2) - r from the edge.
"""
import numpy as np


def blend_patches(mesh, k, tol, off, log=print):
    """Connected groups of off-axis faces (mask off) that are blends: a group that is not one
    cylinder, or whose cylinder (with the other groups about the same axis) covers less than 150
    degrees. Returns (patches, member) with each patch an array of face indices."""
    import trimesh
    from .prism25 import _perp, _circle_fit
    member = np.zeros(len(mesh.faces), bool)
    if not off.any():
        return [], member, []
    adj = mesh.face_adjacency
    keep = off[adj[:, 0]] & off[adj[:, 1]] & (mesh.face_adjacency_angles < np.radians(30))
    comps = [np.asarray(c) for c in trimesh.graph.connected_components(adj[keep], nodes=np.nonzero(off)[0], min_len=1)]
    fits = []
    for comp in comps:
        n, w = mesh.face_normals[comp], mesh.area_faces[comp]
        _, vec = np.linalg.eigh((n * w[:, None]).T @ n)
        d = vec[:, 0]
        fit = None
        if np.abs(n @ d).max() < np.sin(np.radians(2.0)) and len(comp) > 2:
            V = mesh.vertices[np.unique(mesh.faces[comp])]
            u, v = _perp(d)
            c2, r, err = _circle_fit(np.c_[V @ u, V @ v])
            if err <= tol:
                fit = (c2[0] * u + c2[1] * v, d, float(r), u, v)
        fits.append(fit)
    out, cones = [], []
    for comp, fit in zip(comps, fits):
        cn = None if fit is not None else _cone(mesh, comp, tol)
        if cn is not None:
            cones.append(cn)
            member[comp] = True
            continue
        if fit is not None:
            c, d, r, u, v = fit
            # every group about the same axis with the same radius: the cylinder's coverage
            ang = []
            for comp2, f2 in zip(comps, fits):
                if f2 is None or np.linalg.norm(np.cross(f2[1], d)) > 1e-3 or abs(f2[2] - r) > tol \
                        or np.linalg.norm(np.cross(f2[0] - c, d)) > tol:
                    continue
                V = mesh.vertices[np.unique(mesh.faces[comp2])] - c
                ang.append(np.arctan2(V @ v, V @ u))
            a = np.sort(np.mod(np.concatenate(ang), 2 * np.pi))
            span = 2 * np.pi - np.diff(np.r_[a, a[0] + 2 * np.pi]).max()
            if np.degrees(span) >= 150:
                continue
        out.append(comp)
        member[comp] = True
    return out, member, cones


def _cone(mesh, comp, tol):
    """A cone through the faces (normals at one angle to an axis): dict(faces, d, apex, s) with
    s = n . d for the outward normals n; None if the faces are not on one cone."""
    if len(comp) < 4:
        return None
    n, w = mesh.face_normals[comp], mesh.area_faces[comp]
    m = (n * w[:, None]).sum(0) / w.sum()
    C = ((n - m) * w[:, None]).T @ (n - m)
    _, vec = np.linalg.eigh(C)
    d = vec[:, 0]
    if m @ d < 0:
        d = -d
    s_ = n @ d
    if np.abs(s_ - s_.mean()).max() > np.sin(np.radians(1.0)) or not 0.05 < abs(s_.mean()) < 0.98:
        return None
    X = mesh.triangles_center[comp]
    A, *_ = np.linalg.lstsq(n, np.einsum("ij,ij->i", n, X), rcond=None)
    V = mesh.vertices[np.unique(mesh.faces[comp])]
    # every vertex on the cone: distance from the axis proportional to the distance from the apex
    t = (V - A) @ d
    rad = np.linalg.norm(V - A - np.outer(t, d), axis=1)
    ca = np.sqrt(1 - s_.mean() ** 2)
    if np.abs(rad * ca - np.abs(t) * abs(s_.mean())).max() > tol:
        return None
    return dict(faces=comp, d=d, apex=A, s=float(s_.mean()))


def attach_tips(cones, cyls, tol, log=print):
    """Drill points: a cone coaxial with a blind hole, its base the hole's circle at the closed
    end, its apex beyond. Sets cy["tip"] (apex) and returns the cones left over."""
    rest = []
    for cn in cones:
        for cy in cyls:
            if cy["op"] != "cut":
                continue
            c, d = np.asarray(cy["c"]), np.asarray(cy["d"])
            if np.linalg.norm(np.cross(cn["d"], d)) > 1e-3 or np.linalg.norm(np.cross(cn["apex"] - c, d)) > tol:
                continue
            ta = float((cn["apex"] - c) @ d)
            h = ta - cy["traw"][1]
            ang = np.degrees(np.arctan2(cy["r"], h)) if h > 0 else 0.0
            if h <= 0:
                continue
            from .snap2d import nice
            ang = nice(2 * ang, 0.25) / 2                  # included angle to a clean value
            h = cy["r"] / np.tan(np.radians(ang))
            cy["tip"] = dict(h=float(h), angle=float(2 * ang))
            cy["t1"], cy["through"] = cy["traw"][1], False     # the hole ends where the point starts
            log(f"  drill point: {2 * ang:g} deg on the hole d={2 * cy['r']:g}")
            break
        else:
            rest.append(cn)
    return rest


def _edge_samples(e, n=7):
    return np.array([e.positionAt(t).toTuple() for t in np.linspace(0.15, 0.85, n)])


def _poly(e, n=60):
    return np.array([e.positionAt(t).toTuple() for t in np.linspace(0.0, 1.0, n)])


def _curve_distance(e, X, poly):
    """Distance of points X from the edge's curve extended past its ends (a line, a circle),
    else from its polyline."""
    t = e.geomType()
    if t == "LINE":
        p0 = np.array(e.startPoint().toTuple())
        d = np.array(e.endPoint().toTuple()) - p0
        d /= np.linalg.norm(d)
        Q = X - p0
        return np.linalg.norm(Q - np.outer(Q @ d, d), axis=1)
    if t == "CIRCLE":
        c = np.array(e.arcCenter().toTuple())
        r = e.radius()
        n = np.cross(poly[len(poly) // 3] - c, poly[2 * len(poly) // 3] - c)
        if np.linalg.norm(n) < 1e-12:
            n = np.cross(poly[1] - c, poly[len(poly) // 2] - c)
        n /= np.linalg.norm(n)
        Q = X - c
        h = Q @ n
        rho = np.linalg.norm(Q - np.outer(h, n), axis=1)
        return np.sqrt(h ** 2 + (rho - r) ** 2)
    return np.linalg.norm(X[:, None, :] - poly[None, :, :], axis=2).min(axis=1)


def _on_face(f, X, eps):
    """Points of X on the surface of face f (within eps)."""
    import cadquery as cq
    if f.geomType() == "PLANE":
        o = np.array(f.Center().toTuple())
        n = np.array(f.normalAt(f.Center()).toTuple())
        return np.abs((X - o) @ n) < eps
    sel = np.zeros(len(X), bool)
    idx = np.arange(len(X))
    if len(idx) > 300:
        idx = idx[np.linspace(0, len(idx) - 1, 300).astype(int)]
    surf = f._geomAdaptor()
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.gp import gp_Pnt
    for i_ in idx:
        pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*map(float, X[i_])), surf)
        sel[i_] = pr.NbPoints() > 0 and pr.LowerDistance() < eps
    return sel


def find_blends(mesh, solid, patches, tol, log=print):
    """For each blend patch, the edges of the sharp solid it rounds or cuts and its kind and size.
    Returns a list of dict(kind, size, points) (points: edge midpoints, to select the edges), or
    None when a patch has no edge or no consistent size."""
    from trimesh.proximity import closest_point
    from .snap2d import nice
    edges = solid.Edges()
    S = [_edge_samples(e) for e in edges]
    allS = np.vstack(S)
    cp, dist, tri = closest_point(mesh, allS)
    nS = len(S[0])
    out = []
    for comp in patches:
        inpatch = np.isin(tri, comp)
        sel = [m for m in range(len(edges)) if (dist[m * nS:(m + 1) * nS] > tol / 2).all()
               and inpatch[m * nS:(m + 1) * nS].all()]
        if not sel:
            log(f"  blends: no sharp edge under a blend of {len(comp)} faces")
            return None
        V = mesh.vertices[np.unique(mesh.faces[comp])]
        polys = [_poly(edges[m]) for m in sel]
        # each patch vertex to its nearest edge; the leg is the largest distance from the edge's
        # curve extended past its ends (a line, a full circle), so that the vertices of mitred
        # corners, beyond the edge ends, count at their distance across the edge
        best = np.full(len(V), np.inf)
        own = np.full(len(V), -1)
        for q, P in enumerate(polys):
            dm = np.linalg.norm(V[:, None, :] - P[None, :, :], axis=2).min(axis=1)
            better = dm < best
            best[better], own[better] = dm[better], q
        found = []
        for q, m in enumerate(sel):
            mine = own == q
            if not mine.any():
                continue
            e = edges[m]
            mid = e.positionAt(0.5)
            fs = e.ancestors(solid, "Face").Faces()
            if len(fs) != 2:
                return None
            # the leg on each face: the farthest patch vertex lying on that face (where the blend
            # meets it), measured across the edge's extended curve
            Vm = V[mine]
            dc = _curve_distance(e, Vm, polys[q])
            legs = []
            for f in fs:
                on = _on_face(f, Vm, tol / 4)
                if on.any():
                    legs.append(float(dc[on].max()))
            Lg = float(max(legs)) if legs else float(dc.max())
            import os
            if os.environ.get("BLEND_DEBUG"):
                print("   edge", e.geomType(), np.round(mid.toTuple(), 2), "legs", np.round(legs, 4), "n", len(Vm))
            n1, n2 = (np.array(f.normalAt(mid).toTuple()) for f in fs)
            beta = np.pi - np.arccos(np.clip(n1 @ n2, -1.0, 1.0))
            dmid = float(np.median(dist[m * nS:(m + 1) * nS]))
            ch = np.cos(beta / 2)
            fi = np.tan(beta / 2) * (1 / np.sin(beta / 2) - 1)
            kind = "chamfer" if abs(dmid / Lg - ch) < abs(dmid / Lg - fi) else "fillet"
            size = Lg if kind == "chamfer" else Lg * np.tan(beta / 2)
            found.append((kind, size, tuple(np.array(mid.toTuple()))))
        if not found:
            return None
        kinds = {f[0] for f in found}
        if len(kinds) > 1:
            log(f"  blends: a blend of {len(comp)} faces reads as both fillet and chamfer")
            return None
        sz = float(np.median([f[1] for f in found]))
        from .snap2d import nice_any
        size = nice_any(sz, tol / 4)
        # every edge under the patch, also one whose vertices all went to a neighbour
        pts = [tuple(np.array(edges[m].positionAt(0.5).toTuple())) for m in sel]
        out.append(dict(kind=found[0][0], size=size, points=pts))
        log(f"  blends: {found[0][0]} {size:g} mm on {len(pts)} edge(s) (measured {sz:.4f})")
    return out


def edges_at(solid, points):
    """The edges of solid nearest to the points (edge midpoints)."""
    import cadquery as cq
    es = solid.Edges()
    sel = []
    for p in points:
        v = cq.Vertex.makeVertex(*p)
        sel.append(min(es, key=lambda e: e.distance(v)))
    return sel


def apply_blends(solid, blends):
    for b in blends:
        es = edges_at(solid, b["points"])
        if b["kind"] == "fillet":
            solid = solid.fillet(b["size"], es)
        else:
            solid = solid.chamfer(b["size"], None, es)
    return solid
