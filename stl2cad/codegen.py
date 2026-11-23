"""CadQuery script generation for exact-mode results.

Two forms, tried in this order; a script is written only after it has been run and its solid
compared with the reconstruction:

1. Feature script: when the solid is a prism (two parallel planar caps; every other face a plane
   or cylinder parallel to the extrusion direction and spanning the full thickness), a parametric
   script in the free-function style of cadquery.func: the cap outline as lines and arcs, extruded
   by `thickness`, minus round holes listed with their diameters and inner cut-outs.
2. Construction script: for every other result, the fitted surfaces as an editable table, the
   bounding box, one interior point per kept cell and the fillet steps, with verbatim copies of the
   construction helpers, so the script runs with CadQuery alone.
"""
import inspect
import os
import subprocess
import sys
import tempfile
import textwrap

import math

import numpy as np

from . import brep_helpers as BH

# ----------------------------------------------------------------------------- formatting
def _f(x, nd=4):
    v = round(float(x), nd)
    if v == 0:
        v = 0.0
    s = f"{v:.{nd}f}".rstrip("0")
    return s + "0" if s.endswith(".") else s


def _exact(x):
    """Full-precision float literal: tangent and coincident surfaces depend on exact values, and
    rounding (even to 1e-6 mm) changes the arrangement."""
    v = float(x)
    return repr(0.0 if v == 0 else v)


def _lit(v, exact=True):
    f = _exact if exact else _f
    if isinstance(v, dict):
        return "dict(" + ", ".join(f"{k}={_lit(x, exact)}" for k, x in v.items()) + ")"
    if isinstance(v, (str, bool, np.bool_)) or v is None:
        return repr(bool(v)) if isinstance(v, np.bool_) else repr(v)
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return f(v)
    a = np.asarray(v, dtype=object) if not isinstance(v, np.ndarray) else v
    if isinstance(v, (list, tuple)) and v and isinstance(v[0], str):
        return "(" + ", ".join(_lit(x, exact) for x in v) + ")"
    if isinstance(v, np.ndarray) and v.dtype.kind in "fiu":
        if v.ndim == 0:
            return f(v)
        return "(" + ", ".join(_lit(x, exact) for x in v) + ("," if len(v) == 1 else "") + ")"
    return "(" + ", ".join(_lit(x, exact) for x in v) + ("," if len(v) == 1 else "") + ")"


# ----------------------------------------------------------------------------- verification
def verify(script, ref_shape, tol, timeout=900):
    """Run the script in a fresh interpreter and compare its solid with ref_shape.
    Returns (ok, message)."""
    import cadquery as cq
    import trimesh
    from .core import shape_to_mesh
    d = tempfile.mkdtemp(prefix="stl2cad_script_")
    py, out = os.path.join(d, "part.py"), os.path.join(d, "part.step")
    with open(py, "w") as fh:
        fh.write(script)
    try:
        r = subprocess.run([sys.executable, py, out], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, "script timed out"
    if r.returncode != 0 or not os.path.exists(out):
        tail = (r.stderr or r.stdout).strip().splitlines()[-1:] or ["no output"]
        return False, "script failed: " + tail[0][:200]
    s = cq.importers.importStep(out).val()
    BH.clean_mesh(s)
    if not s.isValid():
        return False, "script solid is invalid"
    from .core import mass_props
    v0, v1 = abs(mass_props(ref_shape)[0]), abs(mass_props(s)[0])
    dv = abs(v1 - v0) / v0
    if dv > 1e-4:
        return False, f"volume differs by {100 * dv:.3f} %"
    # tessellation vertices lie on the exact surfaces: measure them against the other B-rep
    da = _max_distance(ref_shape, s, tol)
    db = _max_distance(s, ref_shape, tol)
    dmax = max(da, db)
    if dmax > 0.5 * tol:
        return False, f"surfaces differ by up to {dmax:.4f} mm"
    return True, f"verified: volume {100 * dv:.4f} %, max distance {dmax:.4f} mm"


def _max_distance(own, other, tol, n=800, seed=0):
    """Largest distance from area-weighted samples of `own` to `other`. Samples are drawn from a
    tessellation of `own` and projected onto its exact surface, so chord sag does not count, and
    zero-area slivers carry no weight."""
    import trimesh
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.gp import gp_Pnt
    from .core import shape_to_mesh
    import cadquery as cq
    m = shape_to_mesh(own, tolerance=0.25 * tol)
    P, _ = trimesh.sample.sample_surface_even(m, n, seed=seed)
    # distances to the boundary faces: to a solid, a point inside it is at distance 0
    own_s = cq.Compound.makeCompound(own.Faces()).wrapped
    other_s = cq.Compound.makeCompound(other.Faces()).wrapped
    # the exact distance can overestimate on a face merged by clean() with an irregular boundary
    # (lever: one R16 cylinder face over two height ranges gave 1.7 mm for a point on it), so the
    # smaller of it and the distance to a fine tessellation (deflection 0.05 tol) is taken
    other_tm = shape_to_mesh(other, tolerance=0.05 * tol)
    Q, D = [], []
    for p in P:
        v = BRepBuilderAPI_MakeVertex(gp_Pnt(*map(float, p))).Vertex()
        d0 = BRepExtrema_DistShapeShape(v, own_s)
        if not d0.IsDone():
            continue
        q = d0.PointOnShape2(1)
        d1 = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(q).Vertex(), other_s)
        if d1.IsDone():
            Q.append((q.X(), q.Y(), q.Z()))
            D.append(d1.Value())
    if not Q:
        return 0.0
    from .core import closest_distance
    dt = closest_distance(other_tm, np.array(Q))
    return float(np.max(np.minimum(np.array(D), dt)))


# ----------------------------------------------------------------------------- feature script
def _face_frame(face):
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    pl = BRepAdaptor_Surface(face.wrapped).Plane()
    n = np.array(pl.Axis().Direction().Coord())
    return n, float(n @ np.array(pl.Location().Coord()))


def _ordered_edges(wire):
    """Edges of a wire in traversal order with their start and end points."""
    from OCP.BRepTools import BRepTools_WireExplorer
    from OCP.BRep import BRep_Tool
    from OCP.TopoDS import TopoDS
    import cadquery as cq
    ex = BRepTools_WireExplorer(wire.wrapped)
    out = []
    while ex.More():
        e = cq.Edge(ex.Current())
        v0 = np.array(BRep_Tool.Pnt_s(TopoDS.Vertex_s(ex.CurrentVertex())).Coord())
        p0, p1 = np.array(e.startPoint().toTuple()), np.array(e.endPoint().toTuple())
        a, b = (p0, p1) if np.linalg.norm(p0 - v0) <= np.linalg.norm(p1 - v0) else (p1, p0)
        out.append((e, a, b))
        ex.Next()
    return out


def prism_script(shape, name, tol):
    """Feature script for a prismatic solid, or None."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    if len(shape.Solids()) != 1:
        return None
    faces = shape.Faces()
    planes = [f for f in faces if f.geomType() == "PLANE"]
    eps = 1e-7
    for cap in planes:
        d, _ = _face_frame(cap)
        k = int(np.argmax(np.abs(d)))
        if abs(abs(d[k]) - 1) > eps:
            continue                                   # world-aligned extrusion directions only
        d = np.zeros(3)
        d[k] = 1.0
        levels, caps, ok = [], [], True
        for f in faces:
            t = f.geomType()
            if t == "PLANE":
                n, off = _face_frame(f)
                if abs(abs(n @ d) - 1) < eps:
                    caps.append((f, off * np.sign(n @ d)))
                elif abs(n @ d) > eps:
                    ok = False
            elif t == "CYLINDER":
                ax = np.array(BRepAdaptor_Surface(f.wrapped).Cylinder().Axis().Direction().Coord())
                if abs(abs(ax @ d) - 1) > eps:
                    ok = False
            else:
                ok = False
            if not ok:
                break
        if not ok:
            continue
        offs = sorted({round(o, 6) for _, o in caps})
        if len(offs) != 2:
            continue
        z0, z1 = offs
        h = z1 - z0
        # every side face spans the full thickness
        for f in faces:
            if f.geomType() == "PLANE" and any(f is c for c, _ in caps):
                continue
            zs = [v.toTuple()[k] for v in f.Vertices()]   # bounding boxes of curved faces are padded
            lo, hi = min(zs), max(zs)
            if abs(lo - z0) > 1e-6 * max(1, abs(h)) or abs(hi - z1) > 1e-6 * max(1, abs(h)):
                ok = False
                break
        bottom = [c for c, o in caps if abs(o - z0) < 1e-6]
        if not ok or len(bottom) != 1:
            continue
        code = _prism_code(bottom[0], k, z0, h, name)
        if code is not None:
            return code
    return None


def _prism_code(face, k, z0, h, name):
    axes = "xyz"
    i, j = [m for m in range(3) if m != k]
    names = (axes[i], axes[j])

    def P(p):
        return f"P({_f(p[i])}, {_f(p[j])})"

    def edge_code(e, a, b):
        t = e.geomType()
        if t == "LINE":
            return f"cf.segment({P(a)}, {P(b)})"
        if t == "CIRCLE":
            if e.IsClosed():
                from OCP.BRepAdaptor import BRepAdaptor_Curve
                circ = BRepAdaptor_Curve(e.wrapped).Circle()
                c = np.array(circ.Location().Coord())
                return f"Edge.makeCircle({_f(circ.Radius())}, {P(c)}, DIR)"
            m = np.array(e.positionAt(0.5).toTuple())
            return f"Edge.makeThreePointArc({P(a)}, {P(m)}, {P(b)})"
        return None

    outer = face.outerWire()
    out_edges = []
    for e, a, b in _ordered_edges(outer):
        c = edge_code(e, a, b)
        if c is None:
            return None
        out_edges.append(c)
    holes, cutouts = [], []
    for w in face.innerWires():
        es = w.Edges()
        if len(es) == 1 and es[0].geomType() == "CIRCLE" and es[0].IsClosed():
            from OCP.BRepAdaptor import BRepAdaptor_Curve
            circ = BRepAdaptor_Curve(es[0].wrapped).Circle()
            c = np.array(circ.Location().Coord())
            holes.append((c[i], c[j], 2 * circ.Radius()))
            continue
        cs = []
        for e, a, b in _ordered_edges(w):
            cc = edge_code(e, a, b)
            if cc is None:
                return None
            cs.append(cc)
        cutouts.append(cs)
    diam = sorted({round(dd, 4) for _, _, dd in holes})
    dname = {dd: (f"hole_d{n + 1}" if len(diam) > 1 else "hole_d") for n, dd in enumerate(diam)}
    D = ["0", "0", "0"]
    D[k] = "1"
    pt = ["0", "0", "0"]
    pt[i], pt[j], pt[k] = "a", "b", "z"
    L = []
    L.append(f'"""{name}: parametric CadQuery model generated by stl2cad (exact mode).\n')
    L.append(f"The part is a prism: the outline below, in the {names[0]}{names[1]} plane, extruded by")
    L.append(f"`thickness` along +{axes[k]}, minus the holes and cut-outs. Units: mm.")
    L.append(f"Run:  python {name}.py [out.step]\n\"\"\"\n")
    L.append("import sys\n")
    L.append("from cadquery import Edge")
    L.append("from cadquery import func as cf\n")
    L.append("# ---- parameters ----------------------------------------------------------")
    L.append(f"thickness = {_f(h)}          # extrusion length along +{axes[k]}")
    L.append(f"base = {_f(z0)}               # {axes[k]} of the bottom face")
    for dd in diam:
        L.append(f"{dname[dd]} = {_f(dd)}")
    L.append("")
    L.append(f"DIR = ({', '.join(D)})")
    L.append("")
    L.append("")
    L.append(f"def P(a, b, z=None):")
    L.append(f'    """Point ({names[0]}, {names[1]}) on the bottom face."""')
    L.append(f"    z = base if z is None else z")
    L.append(f"    return ({', '.join(pt)})")
    L.append("")
    L.append("")
    L.append("# ---- outline (lines and arcs) ---------------------------------------------")
    L.append("outline = cf.wire(")
    for c in out_edges:
        L.append(f"    {c},")
    L.append(")")
    L.append("part = cf.extrude(cf.face(outline), [thickness * c for c in DIR])")
    if holes:
        L.append("")
        L.append(f"# ---- round through-holes: ({names[0]}, {names[1]}, diameter) ----------------------")
        L.append("holes = [")
        for a, b, dd in holes:
            L.append(f"    ({_f(a)}, {_f(b)}, {dname[round(dd, 4)]}),")
        L.append("]")
        L.append("for a, b, d in holes:")
        L.append("    disk = cf.face(Edge.makeCircle(d / 2, P(a, b, base - 1), DIR))")
        L.append("    part = part - cf.extrude(disk, [(thickness + 2) * c for c in DIR])")
    if cutouts:
        L.append("")
        L.append("# ---- inner cut-outs (lines and arcs) ---------------------------------------")
        for n, cs in enumerate(cutouts):
            L.append(f"cutout_{n + 1} = cf.wire(")
            for c in cs:
                L.append(f"    {c},")
            L.append(")")
            L.append(f"part = part - cf.extrude(cf.face(cutout_{n + 1}), [thickness * c for c in DIR])")
    L.append("")
    L.append('if __name__ == "__main__":')
    L.append(f'    part.exportStep(sys.argv[1] if len(sys.argv) > 1 else "{name}.step")')
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------- construction script
_HELPERS = [BH._pnt, BH._ax3, BH.plane_face, BH._ax2, BH._lateral, BH.cylinder_face, BH.cone_face,
            BH.sphere_face, BH.torus_face, BH.extrusion_face, BH.revolution_face, BH.profile_face,
            BH.facet_faces, BH._split, BH.sew_cells, BH._drop_dangling, BH.clean_mesh, BH.strip_internal,
            BH.fix_shape, BH.fuse_cells]

_RUNTIME = '''
def make_tool(sp):
    k = sp["kind"]
    if k == "plane":
        return plane_face(sp["n"], sp["o"], sp["u"], sp["rect"])
    if k == "cylinder":
        return cylinder_face(sp["c"], sp["a"], sp["r"], sp["h"])
    if k == "cone":
        return cone_face(sp["v"], sp["a"], sp["t"], sp["h"])
    if k == "sphere":
        return sphere_face(sp["c"], sp["r"])
    if k == "torus":
        return torus_face(sp["c"], sp["a"], sp["R"], sp["r"])
    if k == "extrusion":
        return extrusion_face(sp["points"], sp["vec"], sp["closed"])
    if k == "profile":
        return profile_face(sp["edges"], sp["vec"])
    if k == "revolution":
        return revolution_face(sp["points"], sp["c"], sp["a"])
    return facet_faces(sp["tris"])


def split_box(lo, hi, tools):
    """Cut the box by the tool faces; retry with small fuzzy values as exact mode does."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    size, vbox = float(np.linalg.norm(hi - lo)), float(np.prod(hi - lo))
    best = (np.inf, None)
    for fz in (0.0, 1e-8 * size, 1e-6 * size, 1e-5 * size):
        cells = _split(Solid.makeBox(*(hi - lo), Vector(*lo)), tools, fz, True)
        uncut = len(cells) == 1 and len(cells[0].Faces()) == 6
        err = abs(sum(c.Volume() for c in cells) - vbox) / vbox if cells else np.inf
        if cells and not uncut and err < 1e-4:
            return cells
        if not uncut and err < best[0]:
            best = (err, cells)
    if best[0] < 1e-3:
        return best[1]
    raise RuntimeError("splitting failed")


def keep_cells(cells, points):
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN
    keep = []
    for c in cells:
        cl = BRepClass3d_SolidClassifier(c.wrapped)
        for p in points:
            cl.Perform(_pnt(p), 1e-7)
            if cl.State() == TopAbs_IN:
                keep.append(c)
                break
    return keep


def safe_clean(shape):
    """Face unification (stl2cad runs it in a memory-limited child process)."""
    try:
        return shape.clean()
    except Exception:
        return None


def _edge_mid(e):
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    c = BRepAdaptor_Curve(e.wrapped)
    return np.array(c.Value(0.5 * (c.FirstParameter() + c.LastParameter())).Coord())


def apply_fillets(part, steps):
    """Each step fillets the edges nearest the recorded midpoints with their radii."""
    from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
    for step in steps:
        edges = part.Edges()
        mids = np.array([_edge_mid(e) for e in edges])
        mk = BRepFilletAPI_MakeFillet(part.wrapped)
        for r, m in step:
            e = edges[int(np.argmin(np.linalg.norm(mids - np.asarray(m), axis=1)))]
            mk.Add(float(r), e.wrapped)
        mk.Build()
        part = Shape.cast(mk.Shape())
    return part


def exact_tori(part, tori):
    """Rebuild B-spline fillet faces that lie on a listed torus on the exact toroidal surface."""
    if not tori:
        return part
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.BRepTools import BRepTools, BRepTools_ReShape
    from OCP.Geom import Geom_ToroidalSurface
    from OCP.ShapeFix import ShapeFix_Face, ShapeFix_Shape
    from cadquery import Face
    rs = BRepTools_ReShape()
    for f in part.Faces():
        if f.geomType() != "BSPLINE":
            continue
        s = BRepAdaptor_Surface(f.wrapped)
        u0, u1, v0, v1 = BRepTools.UVBounds_s(f.wrapped)
        P = np.array([s.Value(u0 + a * (u1 - u0), v0 + b * (v1 - v0)).Coord()
                      for a in (0.2, 0.5, 0.8) for b in (0.2, 0.5, 0.8)])
        for t in tori:
            c, a = np.asarray(t["c"]), np.asarray(t["a"]) / np.linalg.norm(t["a"])
            W = P - c
            z = W @ a
            rho = np.linalg.norm(W - np.outer(z, a), axis=1)
            if np.abs(np.hypot(rho - t["R"], z) - t["r"]).max() > 1e-4:
                continue
            surf = Geom_ToroidalSurface(_ax3(c, a), float(t["R"]), float(t["r"]))
            mk = BRepBuilderAPI_MakeFace(surf, f.outerWire().wrapped, True)
            for w in f.innerWires():
                mk.Add(w.wrapped)
            fx = ShapeFix_Face(mk.Face())
            fx.Perform()
            new = fx.Face()
            p0 = P[4]
            if np.dot(f.normalAt(Vector(*p0)).toTuple(), Face(new).normalAt(Vector(*p0)).toTuple()) < 0:
                new = new.Reversed()
            rs.Replace(f.wrapped, new)
            break
    out = Shape.cast(rs.Apply(part.wrapped))
    fx = ShapeFix_Shape(out.wrapped)
    fx.Perform()
    return Shape.cast(fx.Shape())


def build():
    tools = [make_tool(sp) for sp in SURFACES]
    cells = split_box(BOX[0], BOX[1], tools)
    part = fuse_cells(keep_cells(cells, INSIDE))
    sols = part.Solids()
    if len(sols) > 1:                    # drop zero-volume slivers, as exact mode does
        big = [x for x in sols if x.Volume() > 1e-6 * max(y.Volume() for y in sols)]
        part = big[0] if len(big) == 1 else cf.compound(big)
    part = apply_fillets(part, FILLETS)
    return exact_tori(part, TORI)


if __name__ == "__main__":
    build().exportStep(sys.argv[1] if len(sys.argv) > 1 else NAME + ".step")
'''


def _spec_line(sp):
    """One surface as a dict literal, one field per line."""
    sp = {k: v for k, v in sp.items() if k not in ("h_region", "separator")}
    if sp["kind"] == "cylinder":
        sp["hole"] = bool(sp.get("hole", False))
    if sp["kind"] == "profile":
        body = ["        edges=("] + [f"            {_lit(e)}," for e in sp["edges"]] + ["        ),"]
        body.append(f"        vec={_lit(np.asarray(sp['vec']))},")
    else:
        body = [f"        {k}={_lit(v)}," for k, v in sp.items() if k != "kind"]
    return "\n".join([f"    dict(kind={sp['kind']!r},"] + body + ["    ),"])


def construction_script(info, name):
    specs = info["specs"]
    L = []
    L.append(f'"""{name}: CadQuery construction script generated by stl2cad (exact mode).\n')
    L.append("The part is rebuilt from the surfaces fitted to the mesh. Each surface in SURFACES is")
    L.append("extended past its region; together they cut the box BOX into cells, and the cells that")
    L.append("contain a point of INSIDE form the part. FILLETS then re-applies the fillets (radius and")
    L.append("edge midpoint) and TORI makes the circular fillets exact tori. Units: mm.")
    L.append("Edit a parameter (a radius r, a plane offset, a cylinder axis) and re-run;")
    L.append("keep each edit small enough that every INSIDE point stays inside the part.")
    L.append(f"Run:  python {name}.py [out.step]\n\"\"\"\n")
    L.append("import sys\n")
    L.append("import numpy as np")
    L.append("from cadquery import Compound, Shape, Solid, Vector")
    L.append("from cadquery import func as cf")
    src = inspect.getsource(BH)
    imports = [ln for ln in src.splitlines() if ln.startswith("from OCP")]
    # multi-line OCP imports: take the block between the first OCP import and TWO_PI
    i0 = src.index("from OCP.BOPAlgo")
    i1 = src.index("TWO_PI")
    L.append(src[i0:i1].rstrip())
    L.append("")
    L.append(f"NAME = {name!r}")
    L.append("TWO_PI = 2 * np.pi\n")
    kinds = {}
    for sp in specs:
        kinds[sp["kind"]] = kinds.get(sp["kind"], 0) + 1
    L.append("# ---- fitted surfaces (" + ", ".join(f"{v} {k}" for k, v in kinds.items()) + ") " + "-" * 10)
    L.append("# plane: normal n, point o, in-plane axis u, rect = (u0, u1, v0, v1)")
    L.append("# cylinder: axis point c, direction a, radius r, axial range h; cone: apex v, half-angle t")
    L.append("# profile: line/arc chain (tangent planes and cylinders) swept along vec")
    L.append("SURFACES = [")
    for sp in specs:
        L.append(_spec_line(sp))
    L.append("]\n")
    lo, hi = info["box"]
    L.append(f"BOX = ({_lit(np.asarray(lo))}, {_lit(np.asarray(hi))})\n")
    L.append("# one interior point of each cell that belongs to the part")
    L.append("INSIDE = [")
    for p in info["points"]:
        L.append(f"    {_lit(np.asarray(p))},")
    L.append("]\n")
    steps = info.get("fillets", [])
    L.append("# fillet operations: (radius, edge midpoint) for each edge")
    L.append("FILLETS = [")
    for st in steps:
        L.append("    [" + ", ".join(f"({_exact(r)}, {_lit(np.asarray(m))})" for r, m in st) + "],")
    L.append("]\n")
    L.append("TORI = [")
    for t in info.get("tori", []):
        L.append(f"    {_lit(t)},")
    L.append("]\n")
    L.append("\n# ---- construction helpers (copied from stl2cad.brep_helpers) " + "-" * 14)
    for fn in _HELPERS:
        L.append("")
        L.append(textwrap.dedent(inspect.getsource(fn)).rstrip())
        L.append("")
    L.append(_RUNTIME)
    return "\n".join(L)


# ----------------------------------------------------------------------------- stacked prisms
def revolve_script(rv, name, nd=4):
    """Feature script for a body of revolution: the profile (lines and arcs in the half plane
    through the axis: u from the axis, v along it) revolved 360 degrees."""
    f = lambda x: _f(x, nd)
    axes = "xyz"
    c, d, e1 = (np.asarray(rv[k], float) for k in ("c", "d", "e1"))
    vec = lambda v: "(" + ", ".join(f(x) for x in v) + ")"
    loops = sorted(rv["loops"], key=lambda L: -abs(L.area))

    def P(q):
        return f"P({f(q[0])}, {f(q[1])})"

    def edges(L):
        out = []
        for e in L.edges:
            if e[0] == "line":
                out.append(f"Edge.makeLine({P(e[1])}, {P(e[2])})")
            else:
                out.append(f"Edge.makeThreePointArc({P(e[1])}, {P(e[2])}, {P(e[3])})")
        return out
    L = [f'"""{name}: parametric CadQuery model generated by stl2cad (exact mode, revolve).\n',
         f"The part is one profile revolved 360 degrees about the {axes[rv['axis']]} axis "
         f"through {vec(c)}. The profile is drawn",
         "with u = distance from the axis and v = position along it. Units: mm.",
         f"Run:  python {name}.py [out.step]\n\"\"\"\n",
         "import sys\n", "from cadquery import Edge, Face, Solid, Vector, Wire\n",
         "# ---- axis ----------------------------------------------------------------",
         f"AXIS_POINT = Vector{vec(c)}",
         f"AXIS_DIR = Vector{vec(d)}",
         f"RADIAL = Vector{vec(e1)}          # the profile's u direction", "", "",
         "def P(u, v):",
         '    """Profile point: u from the axis, v along it."""',
         "    return AXIS_POINT + RADIAL * u + AXIS_DIR * v", "", "",
         "# ---- profile --------------------------------------------------------------",
         "outline = Wire.assembleEdges(["]
    L += [f"    {e}," for e in edges(loops[0])]
    L.append("])")
    inner = []
    for n, Li in enumerate(loops[1:], 1):
        L.append(f"inner_{n} = Wire.assembleEdges([")
        L += [f"    {e}," for e in edges(Li)]
        L.append("])")
        inner.append(f"inner_{n}")
    L += ["", f"profile = Face.makeFromWires(outline, [{', '.join(inner)}])",
          "part = Solid.revolve(profile, 360.0, AXIS_POINT, AXIS_POINT + AXIS_DIR)", ""]
    L += _cut_code(rv, f, vec)
    L += ['if __name__ == "__main__":',
          f'    part.exportStep(sys.argv[1] if len(sys.argv) > 1 else "{name}.step")']
    return "\n".join(L) + "\n"


def _cut_code(rv, f, vec):
    """Straight cuts across a body of revolution: an outline in the plane through the axis
    spanned by ACROSS (w) and AXIS_DIR, extruded along the cut direction through the body (both
    sides) or from the axis outward; identical outlines are one loop over the directions."""
    import numpy as np
    cuts = rv.get("cuts") or []
    if not cuts:
        return []
    key = lambda cu: (cu["both"], tuple((e[0],) + tuple(tuple(np.round(np.asarray(q, float), 6)) for q in e[1:])
                                        for L_ in cu["loops"] for e in L_.edges))
    groups = []
    for cu in cuts:
        for g in groups:
            if key(g[0]) == key(cu):
                g.append(cu)
                break
        else:
            groups.append([cu])
    reach = 1.2 * max(abs(float(q[0])) for L_ in rv["loops"] for e in L_.edges for q in e[1:])
    L = ["# ---- straight cuts across the body: outline (w across the cut, v along the axis) "
         "extruded along the cut direction ----",
         f"REACH = {f(reach)}                 # beyond the body's largest radius", "", "",
         "def cut_prism(direction, across, outline, through):",
         '    """Outline in the plane through the axis spanned by across and AXIS_DIR, extruded along',
         '    direction: through the body (both sides) or from the axis outward."""',
         "    Q = lambda w, v: AXIS_POINT + Vector(*across) * w + AXIS_DIR * v",
         "    face = Face.makeFromWires(Wire.assembleEdges([",
         "        Edge.makeLine(Q(*a), Q(*b)) for a, b in outline]), [])",
         "    D = Vector(*direction)",
         "    if through:",
         "        face = face.translate(D * -REACH)",
         "    return Solid.extrudeLinear(face, D * (2 * REACH if through else REACH))", ""]
    for n, g in enumerate(groups, 1):
        cu = g[0]
        if any(e[0] != "line" for L_ in cu["loops"] for e in L_.edges):
            raise ValueError("a cut outline with arcs")
        L.append(f"# cut {n}" + (f" ({len(g)} directions)" if len(g) > 1 else ""))
        L.append("outline = [")
        L += [f"    (({f(e[1][0])}, {f(e[1][1])}), ({f(e[2][0])}, {f(e[2][1])}))," for L_ in cu["loops"] for e in L_.edges]
        L.append("]")
        L.append("for direction, across in [")
        L += [f"    ({vec(x['D'])}, {vec(x['w'])})," for x in g]
        L.append("]:")
        L.append(f"    part = part - cut_prism(direction, across, outline, {bool(cu['both'])})")
        L.append("")
    return L


def stacked_script(p25, name, nd=4):
    """Feature script for a 2.5D result: a sequence of extrusions that add material and cuts that
    remove it, each an outline (lines and arcs) between two levels; round cuts are listed as
    holes with named diameters."""
    k = p25["axis"]
    feats = p25["features"]
    axes = "xyz"
    i, j = [m for m in range(3) if m != k]
    f = lambda x: _f(x, nd)
    is_hole = lambda ft: ft["op"] == "cut" and len(ft["loop"].edges) == 1 and ft["loop"].edges[0][0] == "circle"

    def P(p):
        return f"P({f(p[0])}, {f(p[1])})"

    def edges(loop, which=None):
        out = []
        for e in (loop.edges if which is None else which):
            if e[0] == "line":
                out.append(f"cf.segment({P(e[1])}, {P(e[2])})")
            elif e[0] == "arc":
                out.append(f"Edge.makeThreePointArc({P(e[1])}, {P(e[2])}, {P(e[3])})")
            elif e[0] == "ellipse":
                deg = math.degrees(float(e[4]))
                xdir = ("P(1, 0)" if abs(deg) < 1e-9 else
                        f"P(math.cos(math.radians({f(deg)})), math.sin(math.radians({f(deg)})))")
                out.append(f"Edge.makeEllipse({f(e[2])}, {f(e[3])}, {P(e[1])}, DIR, {xdir})")
            else:
                out.append(f"Edge.makeCircle({f(e[2])}, {P(e[1])}, DIR)")
        return out
    diams = sorted({round(2 * ft["loop"].edges[0][2], nd) for ft in feats if is_hole(ft)})
    dname = {d: (f"hole_d{n + 1}" if len(diams) > 1 else "hole_d") for n, d in enumerate(diams)}
    zs = sorted({round(ft["z0"], nd) for ft in feats} | {round(ft["z1"], nd) for ft in feats})
    zname = {z: f"z{n}" for n, z in enumerate(zs)}
    Z = lambda z: zname[round(z, nd)]
    D = ["0", "0", "0"]
    D[k] = "1"
    pt = ["0", "0", "0"]
    pt[i], pt[j] = "a", "b"
    n_add = sum(ft["op"] == "add" for ft in feats)
    L = [f'"""{name}: parametric CadQuery model generated by stl2cad (exact mode, 2.5D).\n',
         f"The part is built along +{axes[k]} from {n_add} extrusion(s) and {len(feats) - n_add} cut(s), "
         f"applied in order;",
         f"each is an outline in the {axes[i]}{axes[j]} plane between two levels. Units: mm.",
         f"Run:  python {name}.py [out.step]\n\"\"\"\n",
         "import math\nimport sys\n", "from cadquery import Edge, Vector", "from cadquery import func as cf\n",
         "# ---- parameters ----------------------------------------------------------"]
    for z in zs:
        L.append(f"{zname[z]} = {f(z)}    # level along {axes[k]}")
    for d in diams:
        L.append(f"{dname[d]} = {f(d)}")
    L += ["", f"DIR = ({', '.join(D)})", "", "",
          "def P(a, b):",
          f'    """Point ({axes[i]}, {axes[j]}) in the sketch plane ({axes[k]} = 0)."""',
          f"    return ({', '.join(pt)})", "", "",
          "def prism(outline, lo, hi):",
          '    """Outline extruded from level lo to level hi."""',
          "    body = cf.extrude(cf.face(outline), [(hi - lo) * c for c in DIR])",
          "    return body.translate(Vector(*[lo * c for c in DIR]))", "", "",
          "def hole(a, b, d, lo, hi):",
          "    return prism(Edge.makeCircle(d / 2, P(a, b), DIR), lo, hi)", "", ""]
    if any(ft["loop"].pattern is not None for ft in feats):
        L += ["def turn(p, centre, deg):",
              '    """Point p turned by deg about the axis through centre along DIR."""',
              "    from cadquery import Vertex",
              "    q = Vertex.makeVertex(*p).rotate(Vector(*centre), Vector(*centre) + Vector(*DIR), deg)",
              "    return q.toTuple()", "", ""]
    m = 0
    n = 0
    while m < len(feats):
        ft = feats[m]
        if is_hole(ft):
            grp = [ft]
            while m + len(grp) < len(feats) and is_hole(feats[m + len(grp)]):
                grp.append(feats[m + len(grp)])
            L.append("# ---- round holes: (" + f"{axes[i]}, {axes[j]}, diameter, from, to) " + "-" * 24)
            L.append("for a, b, d, lo, hi in [")
            for g in grp:
                c, r = g["loop"].edges[0][1], g["loop"].edges[0][2]
                L.append(f"    ({f(c[0])}, {f(c[1])}, {dname[round(2 * r, nd)]}, {Z(g['z0'])}, {Z(g['z1'])}),")
            L.append("]:")
            L.append("    part = part - hole(a, b, d, lo, hi)")
            L.append("")
            m += len(grp)
            continue
        n += 1
        verb = "extrusion" if ft["op"] == "add" else "cut"
        L.append(f"# ---- {n}. {verb}: {axes[k]} {f(ft['z0'])} to {f(ft['z1'])} " + "-" * 30)
        pat = ft["loop"].pattern
        if pat is not None:
            pc, N, one = pat
            L.append(f"count_{n} = {N}                  # the outline repeats {N} times about the centre")
            L.append(f"centre_{n} = {P(pc)}")
            vs = [e[1] for e in one]
            L.append(f"v = [                            # vertices of one period")
            L += [f"    {P(v)}," for v in vs]
            L.append("]")
            L.append(f"v.append(turn(v[0], centre_{n}, 360.0 / count_{n}))   # the period ends where the next begins")
            L.append(f"period_{n} = [")
            for m_, e in enumerate(one):
                if e[0] == "line":
                    L.append(f"    cf.segment(v[{m_}], v[{m_ + 1}]),")
                else:
                    L.append(f"    Edge.makeThreePointArc(v[{m_}], {P(e[2])}, v[{m_ + 1}]),")
            L.append("]")
            L.append(f"outline_{n} = cf.wire(*[")
            L.append(f"    e.rotate(Vector(*centre_{n}), Vector(*centre_{n}) + Vector(*DIR), 360.0 * m / count_{n})")
            L.append(f"    for m in range(count_{n}) for e in period_{n}")
            L.append("])")
        elif getattr(ft["loop"], "linear", None) is not None:
            p_, N_, v_ = ft["loop"].linear
            E_ = ft["loop"].edges
            L.append(f"count_{n} = {N_}                  # the period repeats {N_} times, each moved by step_{n}")
            L.append(f"step_{n} = Vector(*{P(v_)})")
            L.append(f"period_{n} = [")
            L += [f"    {c}," for c in edges(ft["loop"], E_[:p_])]
            L.append("]")
            L.append(f"rest_{n} = [                     # the outline after the last copy, back to the start")
            L += [f"    {c}," for c in edges(ft["loop"], E_[p_ * N_:])]
            L.append("]")
            L.append(f"outline_{n} = cf.wire(*[e.translate(step_{n} * m) for m in range(count_{n}) "
                     f"for e in period_{n}], *rest_{n})")
        else:
            L.append(f"outline_{n} = cf.wire(")
            L += [f"    {c}," for c in edges(ft["loop"])]
            L.append(")")
        expr = f"prism(outline_{n}, {Z(ft['z0'])}, {Z(ft['z1'])})"
        if m == 0:
            L.append(f"part = {expr}")
        else:
            L.append(f"part = part {'+' if ft['op'] == 'add' else '-'} {expr}")
        L.append("")
        m += 1
    revs = p25.get("revolved") or []
    if revs:
        L += _revolved_code(revs, k, f, P)
        L[L.index("from cadquery import Edge, Vector")] = "from cadquery import Edge, Face, Solid, Vector, Wire"
    cyls = p25.get("cylinders") or []
    if cyls:
        import numpy as np
        L.append("# ---- off-axis bosses and holes: (start point, axis direction, diameter, length) " + "-" * 4)
        L.append("for p, u, d, h, op in [")
        for cy in cyls:
            c, u = np.asarray(cy["c"]), np.asarray(cy["d"])
            p0 = c + cy["t0"] * u
            L.append(f"    (({f(p0[0])}, {f(p0[1])}, {f(p0[2])}), ({_f(u[0], 10)}, {_f(u[1], 10)}, {_f(u[2], 10)}), "
                     f"{f(2 * cy['r'])}, {f(cy['t1'] - cy['t0'])}, \"{cy['op']}\"),")
        L.append("]:")
        L.append("    c = Solid.makeCylinder(d / 2, h, Vector(*p), Vector(*u))")
        L.append("    part = part + c if op == \"add\" else part - c")
        L.append("")
        if "from cadquery import Edge, Vector" in L:
            L[L.index("from cadquery import Edge, Vector")] = "from cadquery import Edge, Solid, Vector"
    L.append("part = part.clean()")
    L += ["", 'if __name__ == "__main__":',
          f'    part.exportStep(sys.argv[1] if len(sys.argv) > 1 else "{name}.step")']
    return "\n".join(L) + "\n"


def _revolved_code(revs, k, f, P):
    """Script lines for revolved features about axes along k: each profile (u from the axis, v the
    level) revolved 360 degrees and added or cut; identical profiles about several axes are one
    loop over the axis points."""
    import numpy as np
    i, j = [m for m in range(3) if m != k]
    axes = "xyz"
    rad = ["0", "0", "0"]
    rad[i] = "1"

    def key(ft):
        return (ft["op"], ft["mode"], tuple((e[0],) + tuple(tuple(np.round(np.asarray(q, float), 9)) if m else q
                                                            for m, q in enumerate(e[1:], 1))
                                            for L_ in ft["loops"] for e in L_.edges))
    groups = []
    for ft in revs:
        for g in groups:
            if key(g[0]) == key(ft):
                g.append(ft)
                break
        else:
            groups.append([ft])

    def Q(q):
        return f"Q({f(q[0])}, {f(q[1])})"

    def edges(Lp):
        out = []
        for e in Lp.edges:
            if e[0] == "line":
                out.append(f"Edge.makeLine({Q(e[1])}, {Q(e[2])})")
            elif e[0] == "arc":
                out.append(f"Edge.makeThreePointArc({Q(e[1])}, {Q(e[2])}, {Q(e[3])})")
            else:
                raise ValueError("circle in a revolved profile")
        return out
    L = ["# ---- revolved features about axes along " + axes[k] + ": a profile in the half plane through the",
         "# axis (u = distance from the axis, v = level), revolved 360 degrees " + "-" * 14,
         f"RADIAL = ({', '.join(rad)})", "", "",
         "def revolved(axis, rings):",
         '    """Profile rings (outline first, then holes) revolved about the axis through axis along DIR."""',
         "    face = Face.makeFromWires(Wire.assembleEdges(rings[0]), [Wire.assembleEdges(r) for r in rings[1:]])",
         "    return Solid.revolve(face, 360.0, Vector(*axis), Vector(*axis) + Vector(*DIR))", ""]
    for n, g in enumerate(groups, 1):
        ft = g[0]
        what = {"add": "boss", "cut": "cut"}[ft["op"]]
        pts = [ft_["a3"] for ft_ in g]
        L.append(f"# revolved {what} {n}" + (f" ({len(g)} axes)" if len(g) > 1 else ""))
        L.append("for a, b in [" + ", ".join(f"({f(p[i])}, {f(p[j])})" for p in pts) + "]:")
        L.append("    axis = P(a, b)")
        L.append("    Q = lambda u, v: Vector(*axis) + Vector(*RADIAL) * u + Vector(*DIR) * v")
        L.append("    rings = [")
        for Lp in ft["loops"]:
            L.append("        [")
            L += [f"            {e}," for e in edges(Lp)]
            L.append("        ],")
        L.append("    ]")
        if ft["op"] == "cut":
            L.append("    part = part - revolved(axis, rings)")
        else:
            L.append("    part = part + revolved(axis, rings)")
        L.append("")
    return L


def generate(res, name, tol, log=print):
    """Return (script, kind, message); script is None if no form could be verified."""
    if res.mode != "exact" or res.shape is None:
        return None, None, "scripts are generated for exact mode only"
    msgs = []
    code = None
    rv = res.extras.get("revolve")
    if rv is not None:
        for nd in (4, 10):
            code = revolve_script(rv, name, nd)
            ok, msg = verify(code, res.shape, tol)
            if ok:
                return code, "feature", msg
        msgs.append("revolve script " + msg)
        code = None
    p25 = res.extras.get("prism25")
    if p25 is not None:
        for nd in (4, 10):                       # readable values first, then full precision
            code = stacked_script(p25, name, nd)
            ok, msg = verify(code, res.shape, tol)
            if ok:
                return code, "feature", msg
        msgs.append("stacked-prism script " + msg)
        code = None
    try:
        code = prism_script(res.shape, name, tol)
    except Exception as e:                       # a detector failure falls through to the next form
        msgs.append(f"feature script: {e}")
    if code is not None:
        ok, msg = verify(code, res.shape, tol)
        if ok:
            return code, "feature", msg
        msgs.append("feature script " + msg)
    info = res.extras.get("brep")
    if info is not None:
        code = construction_script(info, name)
        ok, msg = verify(code, res.shape, tol)
        if ok:
            return code, "construction", msg
        msgs.append("construction script " + msg)
    return None, None, "; ".join(msgs) or "no script form applies"
