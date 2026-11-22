"""Blend (fillet) handling for exact mode.

Tangent contacts between surfaces leave gaps in the cell arrangement. A blend surface is therefore
removed before the arrangement is built: the two surfaces it connects are extended until they meet
in a sharp edge, the sharp solid is built, and the blend is re-created on that edge with the OCC
fillet operation and the fitted radius, which is how such parts are modelled in CAD.

A blend is a cylinder or torus region that spans less than a half turn and whose spine (the
cylinder axis, the torus tube-centre circle) lies at distance r from two neighbouring surfaces
along its whole length, i.e. the region is tangent to both.
"""
from dataclasses import dataclass

import numpy as np

from .primitives import perp_basis


@dataclass
class Blend:
    index: int          # position in the primitive list before removal
    kind: str
    r: float
    nbrs: tuple         # the two tangent neighbours (indices before removal)
    region: np.ndarray  # region points


def _angular_span(ang):
    """Angular extent (rad) of a set of angles: 2*pi minus the largest gap."""
    a = np.sort(np.mod(ang, 2 * np.pi))
    gaps = np.diff(np.r_[a, a[0] + 2 * np.pi])
    return 2 * np.pi - gaps.max()


def _spine(pr):
    """Sample points of the spine of a cylinder or torus region, and the region's angular span."""
    X, p = pr.pts, pr.p
    a = p["a"]
    if pr.kind == "cylinder":
        h = (X - p["c"]) @ a
        S = p["c"] + np.outer(np.linspace(h.min(), h.max(), 7), a)
        u, v = perp_basis(a)
        W = X - p["c"]
        return S, _angular_span(np.arctan2(W @ v, W @ u))
    # torus: tube-centre circle points in the radial directions of the region points
    W = X - p["c"]
    z = W @ a
    radial = W - np.outer(z, a)
    rho = np.linalg.norm(radial, axis=1)
    e = radial / np.maximum(rho, 1e-12)[:, None]
    idx = np.linspace(0, len(X) - 1, min(len(X), 24)).astype(int)
    S = p["c"] + p["R"] * e[idx]
    span = _angular_span(np.arctan2(z, rho - p["R"]))      # angle around the tube
    return S, span


def _region_adjacency(prims, mesh):
    own = np.full(len(mesh.faces), -1)
    for i, p in enumerate(prims):
        own[p.faces] = i
    adj = mesh.face_adjacency
    a, b = own[adj[:, 0]], own[adj[:, 1]]
    m = (a >= 0) & (b >= 0) & (a != b)
    count = {}
    for x, y in zip(a[m], b[m]):
        k = (min(x, y), max(x, y))
        count[k] = count.get(k, 0) + 1
    nb = {}
    for (x, y), c in count.items():
        nb.setdefault(x, {})[y] = c
        nb.setdefault(y, {})[x] = c
    return nb


def detect_blends(prims, mesh, tol, exclude=(), max_span_deg=175.0, log=print):
    nb = _region_adjacency(prims, mesh)
    blends = []
    for i, pr in enumerate(prims):
        if pr.kind not in ("cylinder", "torus") or i in exclude or i not in nb:
            continue
        if getattr(pr, "is_hole", False) and pr.kind == "cylinder":
            continue
        S, span = _spine(pr)
        if np.degrees(span) > max_span_deg:
            continue
        r = float(pr.p["r"])
        tang = []
        for j, cnt in sorted(nb[i].items(), key=lambda t: -t[1]):
            q = prims[j]
            if q.kind not in ("plane", "cylinder", "cone", "sphere"):
                continue
            d = np.abs(q.distance(S))
            if np.abs(d - r).max() < 3 * tol:
                tang.append(j)
        if len(tang) >= 2:
            blends.append(Blend(i, pr.kind, r, (tang[0], tang[1]), pr.pts))
    if blends:
        log(f"  blends: {len(blends)} fillet surface(s) recognised, radii "
            + ", ".join(sorted({f'{b.r:.3g}' for b in blends})) + " mm")
    return blends


def remove_blends(prims, blends, pairs):
    """Drop blend regions; extend each neighbour's tool extent over the blend region plus r.
    Returns the reduced primitive list, remapped tangent pairs and blends with remapped neighbours."""
    drop = {b.index for b in blends}
    for b in blends:
        for j in b.nbrs:
            q = prims[j]
            X = b.region
            P = X - q.distance(X)[:, None] * q.normal(X)
            q.extent_pts = np.vstack([getattr(q, "extent_pts", q.pts), P])
            q.extra_margin = max(getattr(q, "extra_margin", 0.0), b.r)
    new_index = {}
    out = []
    for i, p in enumerate(prims):
        if i in drop:
            continue
        new_index[i] = len(out)
        out.append(p)
    pairs = [(new_index[a], new_index[b]) for a, b in pairs if a in new_index and b in new_index]
    kept = []
    for b in blends:
        if all(j in new_index for j in b.nbrs):
            kept.append(Blend(b.index, b.kind, b.r, tuple(new_index[j] for j in b.nbrs), b.region))
    return out, pairs, kept


# ----------------------------------------------------------------------------------- OCC side
def _face_on(face, pr, tol):
    """True if an OCC face lies on the surface of primitive pr (tested at points of the face)."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepTools import BRepTools
    s = BRepAdaptor_Surface(face)
    u0, u1, v0, v1 = BRepTools.UVBounds_s(face)
    P = []
    for fu in (0.3, 0.5, 0.7):
        for fv in (0.3, 0.5, 0.7):
            P.append(s.Value(u0 + fu * (u1 - u0), v0 + fv * (v1 - v0)).Coord())
    return np.abs(pr.distance(np.array(P))).max() < 3 * tol


def _edge_points(edge, n=5):
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    c = BRepAdaptor_Curve(edge)
    t0, t1 = c.FirstParameter(), c.LastParameter()
    return np.array([c.Value(t0 + f * (t1 - t0)).Coord() for f in np.linspace(0, 1, n)])


def blend_edges(shape, prims, blends, tol):
    """For each blend, the edges of the sharp solid between faces on its two neighbours, near it."""
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
    from OCP.TopoDS import TopoDS
    from scipy.spatial import cKDTree
    emap = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape.wrapped, TopAbs_EDGE, TopAbs_FACE, emap)
    on_cache = {}

    def on(face, j):
        key = (face.__hash__(), j)
        if key not in on_cache:
            on_cache[key] = _face_on(face, prims[j], tol)
        return on_cache[key]

    out = []
    for b in blends:
        tree = cKDTree(b.region)
        A, B = b.nbrs
        edges = []
        for k in range(1, emap.Extent() + 1):
            faces = list(emap.FindFromIndex(k))
            if len(faces) != 2:
                continue
            f0, f1 = TopoDS.Face_s(faces[0]), TopoDS.Face_s(faces[1])
            if not ((on(f0, A) and on(f1, B)) or (on(f0, B) and on(f1, A))):
                continue
            e = TopoDS.Edge_s(emap.FindKey(k))
            d, _ = tree.query(_edge_points(e))
            if np.median(d) < 1.5 * b.r + 3 * tol:
                edges.append(e)
        out.append(edges)
    return out


def _edge_mid(edge):
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    c = BRepAdaptor_Curve(edge)
    return np.array(c.Value(0.5 * (c.FirstParameter() + c.LastParameter())).Coord())


def apply_fillets(shape, prims, blends, tol, log=print, record=None):
    """Fillet the sharp solid. All blends at once first; on failure, one blend at a time, keeping
    the ones that succeed. record (dict) receives the steps for script generation: "steps", a list
    of fillet operations, each a list of (radius, edge midpoint); "tori", the exact tori that
    replace B-spline fillet faces."""
    if record is not None:
        record.setdefault("steps", [])
        record.setdefault("tori", [])
    import cadquery as cq
    from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet

    groups = blend_edges(shape, prims, blends, tol)
    todo = [(b, es) for b, es in zip(blends, groups) if es]
    missing = sum(1 for es in groups if not es)
    if missing:
        log(f"  fillets: no sharp edge found for {missing} blend(s)")
    if not todo:
        return shape, 0

    def run(cur, items):
        mk = BRepFilletAPI_MakeFillet(cur.wrapped)
        seen = set()
        for b, es in items:
            for e in es:
                h = e.__hash__()
                if h in seen:
                    continue
                seen.add(h)
                mk.Add(b.r, e)
        mk.Build()
        if not mk.IsDone():
            return None
        res = cq.Shape.cast(mk.Shape())
        if not res.isValid() or abs(res.Volume() - cur.Volume()) > 0.2 * abs(cur.Volume()):
            return None
        return res

    try:
        res = run(shape, todo)
    except Exception:
        res = None
    if res is not None:
        log(f"  fillets: {len(todo)} blend(s) applied on {sum(len(es) for _, es in todo)} edge(s)")
        if record is not None:
            record["steps"].append([(float(b.r), _edge_mid(e).tolist()) for b, es in todo for e in es])
        return canonical_tori(res, tol, log, record), len(todo)
    # one at a time: edges must be re-found on the current shape
    cur, done = shape, 0
    for b in blends:
        es = blend_edges(cur, prims, [b], tol)[0]
        if not es:
            continue
        try:
            nxt = run(cur, [(b, es)])
        except Exception:
            nxt = None
        if nxt is not None:
            if record is not None:
                record["steps"].append([(float(b.r), _edge_mid(e).tolist()) for e in es])
            cur, done = nxt, done + 1
    log(f"  fillets: {done} of {len(todo)} blend(s) applied one by one")
    return (canonical_tori(cur, tol, log, record) if done else cur), done


def _face_samples(face, n=7):
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepTools import BRepTools
    from OCP.BRepLProp import BRepLProp_SLProps
    s = BRepAdaptor_Surface(face)
    u0, u1, v0, v1 = BRepTools.UVBounds_s(face)
    P, N = [], []
    for fu in np.linspace(0.05, 0.95, n):
        for fv in np.linspace(0.05, 0.95, n):
            pr = BRepLProp_SLProps(s, u0 + fu * (u1 - u0), v0 + fv * (v1 - v0), 1, 1e-9)
            if not pr.IsNormalDefined():
                continue
            P.append(pr.Value().Coord())
            N.append(pr.Normal().Coord())
    return np.array(P), np.array(N)


def canonical_tori(shape, tol, log=print, record=None):
    """Rebuild B-spline fillet faces that lie on a torus on the exact toroidal surface. OCC returns
    fillets along closed circular edges as periodic B-splines, which do not survive a STEP round
    trip (the face reads back as unorientable); a toroidal face does."""
    import cadquery as cq
    from OCP.Geom import Geom_ToroidalSurface
    from OCP.gp import gp_Ax3, gp_Pnt, gp_Dir
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.BRepTools import BRepTools, BRepTools_ReShape
    from OCP.ShapeFix import ShapeFix_Face, ShapeFix_Shape
    from OCP.TopAbs import TopAbs_REVERSED
    from .primitives import fit, torus_init
    reshape = BRepTools_ReShape()
    n = 0
    tori = []
    for f in shape.Faces():
        if f.geomType() != "BSPLINE":
            continue
        P, N = _face_samples(f.wrapped)
        if len(P) < 12:
            continue
        init = torus_init(P, N)
        if init is None:
            continue
        q = fit("torus", P, init=init)
        if q is None or not np.isfinite(q.rms) or np.abs(q.distance(P)).max() > 1e-3 * tol:
            continue
        c, a = q.p["c"], q.p["a"]
        surf = Geom_ToroidalSurface(gp_Ax3(gp_Pnt(*map(float, c)), gp_Dir(*map(float, a))),
                                    float(q.p["R"]), float(q.p["r"]))
        wires = f.Wires()
        outer = f.outerWire()
        mk = BRepBuilderAPI_MakeFace(surf, outer.wrapped, True)
        for w in wires:
            if not w.wrapped.IsSame(outer.wrapped):
                mk.Add(w.wrapped)
        if not mk.IsDone():
            continue
        sf = ShapeFix_Face(mk.Face())
        sf.Perform()
        new = sf.Face()
        # keep the orientation of the old face: compare normals at a sample point
        nf = cq.Face(new)
        p0 = P[len(P) // 2]
        n_old = np.array(f.normalAt(cq.Vector(*p0)).toTuple())
        n_new = np.array(nf.normalAt(cq.Vector(*p0)).toTuple())
        if n_old @ n_new < 0:
            new = new.Reversed()
        reshape.Replace(f.wrapped, new)
        n += 1
        tori.append(dict(c=[float(x) for x in c], a=[float(x) for x in a], R=float(q.p["R"]), r=float(q.p["r"])))
    if not n:
        return shape
    res = cq.Shape.cast(reshape.Apply(shape.wrapped))
    sfs = ShapeFix_Shape(res.wrapped)
    sfs.Perform()
    res = cq.Shape.cast(sfs.Shape())
    BRepTools.Clean_s(res.wrapped, True)
    if res.isValid() and abs(res.Volume() - shape.Volume()) < 1e-5 * abs(shape.Volume()):
        log(f"  fillets: {n} B-spline fillet face(s) rebuilt as exact tori")
        if record is not None:
            record.setdefault("tori", []).extend(tori)
        return res
    log(f"  fillets: torus rebuild of {n} face(s) rejected (valid={res.isValid()})")
    return shape


def network_cylinders(prims, mesh, blends, tol, candidates):
    """Edge-fillet cylinders that continue a corner torus fillet (same tube radius, axis
    perpendicular to the torus axis, regions adjacent): the straight runs of one fillet network,
    which the fillet operation re-creates together with the corners."""
    nb = _region_adjacency(prims, mesh)
    tori = [b.index for b in blends if b.kind == "torus"]
    out = set()
    for i in candidates:
        p = prims[i]
        if p.kind != "cylinder" or i not in nb:
            continue
        for t in tori:
            T = prims[t]
            if t in nb[i] and abs(p.p["r"] - T.p["r"]) < 2 * tol and abs(p.p["a"] @ T.p["a"]) < 0.02:
                out.add(i)
    return out
