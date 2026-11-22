"""B-rep construction helpers shared by exact mode and the generated CadQuery scripts.

Every function takes plain numbers, lists and tuples, so the generated scripts can
carry a verbatim copy of this module and run without stl2cad installed.
"""
import os

import numpy as np
from cadquery import Compound, Shape, Solid, Vector
from cadquery import func as cf
from OCP.BOPAlgo import BOPAlgo_Splitter
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon,
                                BRepBuilderAPI_MakeWire)
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeCone, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism,
                             BRepPrimAPI_MakeRevol, BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeTorus)
from OCP.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCP.Geom import Geom_Plane
from OCP.GeomAPI import GeomAPI_Interpolate
from OCP.TColgp import TColgp_HArray1OfPnt
from OCP.TopTools import TopTools_ListOfShape
from OCP.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Dir, gp_Pnt, gp_Vec

TWO_PI = 2 * np.pi


def _pnt(p):
    return gp_Pnt(*map(float, p))


def _ax3(o, z, x=None):
    z = np.asarray(z, float) / np.linalg.norm(z)
    if x is None:
        t = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
        x = np.cross(z, t)
    x = np.asarray(x, float)
    x = x - (x @ z) * z
    x /= np.linalg.norm(x)
    return gp_Ax3(_pnt(o), gp_Dir(*map(float, z)), gp_Dir(*map(float, x)))


def plane_face(n, o, u, rect):
    """Rectangle rect = (u0, u1, v0, v1) on the plane through o with normal n and in-plane axis u."""
    return BRepBuilderAPI_MakeFace(Geom_Plane(_ax3(o, n, u)), *map(float, rect), 1e-7).Face()


def _ax2(o, z):
    z = np.asarray(z, float) / np.linalg.norm(z)
    return gp_Ax2(_pnt(o), gp_Dir(*map(float, z)))


def _lateral(solid_maker, kind):
    """The curved face of a primitive solid (it carries a proper seam edge)."""
    shape = Shape.cast(solid_maker.Shape())
    faces = [f for f in shape.Faces() if f.geomType() == kind]
    return faces[0].wrapped


def cylinder_face(c, a, r, h):
    """Full cylinder band of radius r about the axis (c, a), for axial coordinates h = (h0, h1)."""
    a = np.asarray(a, float) / np.linalg.norm(a)
    base = np.asarray(c, float) + h[0] * a
    return _lateral(BRepPrimAPI_MakeCylinder(_ax2(base, a), float(r), float(h[1] - h[0])), "CYLINDER")


def cone_face(v, a, t, h):
    """Cone with apex v, axis a, half-angle t (rad), between axial distances h = (h0, h1) from the apex."""
    a = np.asarray(a, float) / np.linalg.norm(a)
    base = np.asarray(v, float) + h[0] * a
    r0, r1 = h[0] * np.tan(t), h[1] * np.tan(t)
    return _lateral(BRepPrimAPI_MakeCone(_ax2(base, a), float(r0), float(r1), float(h[1] - h[0])), "CONE")


def sphere_face(c, r):
    return _lateral(BRepPrimAPI_MakeSphere(_pnt(c), float(r)), "SPHERE")


def torus_face(c, a, R, r):
    return _lateral(BRepPrimAPI_MakeTorus(_ax2(c, a), float(R), float(r)), "TORUS")


def extrusion_face(points, vec, closed=False):
    """B-spline curve through `points` swept along `vec`."""
    pts = [_pnt(p) for p in points]
    arr = TColgp_HArray1OfPnt(1, len(pts))
    for i, p in enumerate(pts):
        arr.SetValue(i + 1, p)
    itp = GeomAPI_Interpolate(arr, bool(closed), 1e-7)
    itp.Perform()
    curve = itp.Curve()
    edge = BRepBuilderAPI_MakeEdge(curve).Edge()
    return BRepPrimAPI_MakePrism(edge, gp_Vec(*map(float, vec))).Shape()


def revolution_face(points, c, a):
    """B-spline profile through `points` (in a plane containing the axis) revolved 360 deg about (c, a)."""
    arr = TColgp_HArray1OfPnt(1, len(points))
    for i, p in enumerate(points):
        arr.SetValue(i + 1, _pnt(p))
    itp = GeomAPI_Interpolate(arr, False, 1e-7)
    itp.Perform()
    edge = BRepBuilderAPI_MakeEdge(itp.Curve()).Edge()
    ax = gp_Ax1(_pnt(c), gp_Dir(*map(float, np.asarray(a, float) / np.linalg.norm(a))))
    return BRepPrimAPI_MakeRevol(edge, ax, TWO_PI).Shape()


def profile_face(edges, vec):
    """Line/arc profile swept along `vec`. edges: ("line", P, Q) or ("arc", P, M, Q) in 3D."""
    mw = BRepBuilderAPI_MakeWire()
    for e in edges:
        if e[0] == "line":
            seg = GC_MakeSegment(_pnt(e[1]), _pnt(e[2])).Value()
        else:
            seg = GC_MakeArcOfCircle(_pnt(e[1]), _pnt(e[2]), _pnt(e[3])).Value()
        mw.Add(BRepBuilderAPI_MakeEdge(seg).Edge())
    return BRepPrimAPI_MakePrism(mw.Wire(), gp_Vec(*map(float, vec))).Shape()


def facet_faces(tris):
    out = []
    for tri in tris:
        poly = BRepBuilderAPI_MakePolygon(*[_pnt(q) for q in tri], True)
        out.append(BRepBuilderAPI_MakeFace(poly.Wire(), True).Face())
    return out


def _split(box, tools, fuzzy, nondestructive):
    args, tl = TopTools_ListOfShape(), TopTools_ListOfShape()
    args.Append(box.wrapped)
    for t in tools:
        for f in (t if isinstance(t, list) else [t]):
            tl.Append(f)
    sp = BOPAlgo_Splitter()
    sp.SetArguments(args)
    sp.SetTools(tl)
    if fuzzy:
        sp.SetFuzzyValue(fuzzy)
    sp.SetRunParallel(True)
    sp.SetNonDestructive(nondestructive)
    sp.Perform()
    if sp.HasErrors():
        return []
    return Shape.cast(sp.Shape()).Solids()


def _run_worker(args, timeout):
    """Run `python -m stl2cad._worker ...`; independent of the caller's main module."""
    import subprocess
    import sys
    try:
        r = subprocess.run([sys.executable, "-m", "stl2cad._worker", *map(str, args)], timeout=timeout,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return r.returncode == 0
    except subprocess.TimeoutExpired:
        return False


def _split_isolated(box, tools, fuzzy, nondestructive, timeout=900):
    """Run the OCC splitter in a child process: it can crash (segfault) on some inputs."""
    import os
    import tempfile
    from cadquery import Compound
    d = tempfile.mkdtemp()
    src, dst = os.path.join(d, "in.brep"), os.path.join(d, "out.brep")
    flat = [box.wrapped]
    for t in tools:
        flat += t if isinstance(t, list) else [t]
    try:
        Compound.makeCompound([Shape.cast(x) for x in flat]).exportBrep(src)
        if not _run_worker(["split", src, dst, fuzzy, int(nondestructive)], timeout) or not os.path.exists(dst):
            return []
        return Shape.importBrep(dst).Solids()
    finally:
        for f in (src, dst):
            if os.path.exists(f):
                os.remove(f)
        os.rmdir(d)


def split_box(lo, hi, tools, fuzzy=0.0, isolated=True, return_ok=False):
    """Split the box [lo, hi] by all tool faces; returns the cells (solids).
    The OCC splitter occasionally fails (or crashes) on exact configurations; it runs in a child
    process, a sequence of tiny fuzzy values is tried, and a result is accepted only if its cells
    add up to the box volume."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    size = float(np.linalg.norm(hi - lo))
    vbox = float(np.prod(hi - lo))
    cells, ok = [], False
    best = (np.inf, [])
    for fz, nd in [(fuzzy, True), (max(fuzzy, 1e-8 * size), True), (max(fuzzy, 1e-6 * size), True),
                   (max(fuzzy, 1e-5 * size), True)]:
        box = Solid.makeBox(*(hi - lo), Vector(*lo))
        cells = _split_isolated(box, tools, fz, nd) if isolated else _split(box, tools, fz, nd)
        # a failed split may return the uncut box (or a lost cell); valid cells must also sum to the box
        uncut = bool(tools) and len(cells) == 1 and len(cells[0].Faces()) == 6
        err = abs(sum(c.Volume() for c in cells) - vbox) / vbox if cells else np.inf
        if not uncut and err < best[0]:
            best = (err, cells)
        if cells and err < 1e-4 and not uncut:
            ok = True
            break
    if not ok:
        # a small volume defect usually comes from a lost sliver cell and is harmless; an uncut box
        # or a large defect means the split failed
        err, cells = best
        ok = err < 1e-3
        os.environ.get("STL2CAD_DEBUG") and print(f"split volume defect {err:.2e}")
    out = []
    for c in cells:
        if not c.isValid():                      # repair occasional bad faces from the splitter
            try:
                from OCP.ShapeFix import ShapeFix_Shape
                fx = ShapeFix_Shape(c.wrapped)
                fx.Perform()
                f = Shape.cast(fx.Shape())
                if f.isValid() and f.Solids() and f.Volume() > 0:
                    c = f.Solids()[0] if len(f.Solids()) == 1 else f
            except Exception:
                pass
        out.append(c)
    return (out, ok) if return_ok else out


def sew_cells(cells, tol=1e-6):
    """Assemble cells without booleans: faces shared by two cells are internal; the remaining
    faces are sewn into shells and closed into solids."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing
    from OCP.ShapeFix import ShapeFix_Solid
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    count = {}
    for c in cells:
        for f in c.Faces():
            h = f.hashCode()
            count.setdefault(h, [f, 0])
            count[h][1] += 1
    sew = BRepBuilderAPI_Sewing(tol)
    for f, n in count.values():
        if n == 1:
            sew.Add(f.wrapped)
    sew.Perform()
    sewn = sew.SewedShape()
    solids = []
    ex = TopExp_Explorer(sewn, TopAbs_SHELL)
    while ex.More():
        sf = ShapeFix_Solid()
        solids.append(Shape.cast(sf.SolidFromShell(TopoDS.Shell_s(ex.Current()))))
        ex.Next()
    if not solids:
        return None
    return solids[0] if len(solids) == 1 else cf.compound(solids)


def safe_clean(shape, mem_gb=2.5, timeout=120):
    """Unify coplanar / co-surface faces in a child process with a memory limit and a timeout:
    OCC's face unification occasionally exhausts memory on valid inputs."""
    import os
    import tempfile
    d = tempfile.mkdtemp()
    src, dst = os.path.join(d, "in.brep"), os.path.join(d, "out.brep")
    try:
        shape.exportBrep(src)
        if not _run_worker(["clean", src, dst, int(mem_gb * 1e9)], timeout) or not os.path.exists(dst):
            return shape
        out = Shape.importBrep(dst)
        if out.Solids() and abs(out.Volume() - shape.Volume()) < 1e-4 * abs(shape.Volume()):
            return out.Solids()[0] if len(out.Solids()) == 1 else out
        return shape
    except Exception:
        return shape
    finally:
        for f in (src, dst):
            if os.path.exists(f):
                os.remove(f)
        os.rmdir(d)


def _drop_dangling(faces):
    """Drop faces with a free (non-seam) edge, repeatedly: sheets hanging into a solid."""
    from OCP.BRep import BRep_Tool
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
    from OCP.TopoDS import TopoDS
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    if True:
        comp = TopoDS_Compound()
        bld = BRep_Builder()
        bld.MakeCompound(comp)
        for f in faces:
            bld.Add(comp, f.wrapped)
        emap = TopTools_IndexedDataMapOfShapeListOfShape()
        TopExp.MapShapesAndAncestors_s(comp, TopAbs_EDGE, TopAbs_FACE, emap)
        bad = []
        for k in range(1, emap.Extent() + 1):
            e = TopoDS.Edge_s(emap.FindKey(k))
            if BRep_Tool.Degenerated_s(e):
                continue
            uniq = []
            for g in emap.FindFromIndex(k):
                if not any(g.IsSame(x) for x in uniq):
                    uniq.append(g)
            if len(uniq) == 1 and not BRep_Tool.IsClosed_s(e, TopoDS.Face_s(uniq[0])):
                bad.append(uniq[0])
        # single pass: removing a face exposes its neighbours' edges, so iterating would peel a
        # closed shell away
        return [f for f in faces if not any(f.wrapped.IsSame(b) for b in bad)]


def clean_mesh(shape):
    """Drop stored triangulations (left by meshing for classification); stale polygons on edges
    make BRepCheck report otherwise valid shapes as invalid."""
    from OCP.BRepTools import BRepTools
    if shape is None:
        return shape
    try:
        BRepTools.Clean_s(shape.wrapped, True)
    except TypeError:
        BRepTools.Clean_s(shape.wrapped)
    return shape


def strip_internal(shape):
    """Remove INTERNAL faces (dangling pieces of tool faces that end inside a cell) from a solid or
    compound of solids; they make edges non-manifold and the solid invalid."""
    clean_mesh(shape)
    from OCP.TopAbs import TopAbs_INTERNAL, TopAbs_EXTERNAL
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
    from OCP.ShapeFix import ShapeFix_Solid
    from OCP.TopoDS import TopoDS
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer
    out = []
    changed = False
    for so in shape.Solids():
        faces = so.Faces()
        keep = [f for f in faces if f.wrapped.Orientation() not in (TopAbs_INTERNAL, TopAbs_EXTERNAL)]
        keep = _drop_dangling(keep)
        if len(keep) == len(faces):
            out.append(so)
            continue
        sew = BRepBuilderAPI_Sewing(1e-7)
        for f in keep:
            sew.Add(f.wrapped)
        sew.Perform()
        ex = TopExp_Explorer(sew.SewedShape(), TopAbs_SHELL)
        mk = BRepBuilderAPI_MakeSolid()
        n = 0
        while ex.More():
            mk.Add(TopoDS.Shell_s(ex.Current()))
            n += 1
            ex.Next()
        if n == 0:
            out.append(so)
            continue
        fx = ShapeFix_Solid(mk.Solid())
        fx.Perform()
        new = Shape.cast(fx.Solid())
        if new.isValid() and abs(new.Volume() - so.Volume()) < 1e-4 * abs(so.Volume()) + 1e-12:
            out.append(new)
            changed = True
        else:
            out.append(so)
    if not changed:
        return shape
    return out[0] if len(out) == 1 else Compound.makeCompound(out)


def fuse_cells(cells):
    """Glue cells that share faces into one solid and unify coplanar / co-surface faces."""
    cells = [strip_internal(clean_mesh(c)) for c in cells]
    target = sum(c.Volume() for c in cells)
    ok = lambda s: s is not None and s.Solids() and abs(s.Volume() - target) < 1e-3 * target

    def attempt(fn):
        try:
            return fn()
        except Exception:
            return None
    solid = cells[0] if len(cells) == 1 else attempt(lambda: cf.fuse(*cells, glue="full"))
    clean_mesh(solid)
    if not ok(solid) or not solid.isValid():
        sewn = attempt(lambda: sew_cells(cells))
        if ok(sewn):
            solid = sewn
    if not ok(solid):
        solid = attempt(lambda: cf.fuse(*cells))
    if not ok(solid):
        solid = cells[0]
        for c in cells[1:]:
            solid = solid + c
    if not ok(solid) or len(solid.Solids()) != 1 or len(solid.Faces()) > 5000:
        return solid                       # face unification can exhaust memory on broken results
    cleaned = safe_clean(solid)
    solid = cleaned if ok(cleaned) else solid
    solid = strip_internal(solid)
    return fix_shape(solid)


def fix_shape(shape):
    """ShapeFix pass for an invalid result (typically a shell with the wrong orientation)."""
    clean_mesh(shape)
    if shape.isValid():
        return shape
    from OCP.ShapeFix import ShapeFix_Shape
    try:
        sf = ShapeFix_Shape(shape.wrapped)
        sf.Perform()
        new = Shape.cast(sf.Shape())
    except Exception:
        return shape
    if new.isValid() and abs(abs(new.Volume()) - abs(shape.Volume())) < 1e-4 * abs(shape.Volume()) + 1e-12:
        return new
    return shape


def keep_cells(cells, points):
    """The cells that contain at least one of the given interior points, fused."""
    keep = [c for c in cells if any(c.isInside(Vector(*map(float, p)), 1e-9) for p in points)]
    if not keep:
        raise RuntimeError("no cell contains the interior points")
    return fuse_cells(keep)
