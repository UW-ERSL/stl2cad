"""SolidWorks feature-tree script for 2.5D (stacked-prism) results.

The generated Python script drives SolidWorks through its COM API (pywin32, Windows): one sketch
and one Extruded Boss or Extruded Cut per feature, in the order stl2cad applies them, with round
holes of one height range grouped into one cut. The result is a native feature tree, not an
imported body.

Sketch points use SolidWorks' base-plane orientations: Front (x, y), Top (x, -z), Right (-z, y).
The script does not assume SolidWorks' default extrusion directions: after every feature the
part's volume and centre of mass are compared with the values stl2cad expects at that step; if
they differ, the feature is deleted and rebuilt with the other direction (and the other
start-offset side), so a wrong guess cannot pass silently. Each sketch is first drawn with
AddToDB; if SolidWorks makes no feature from it, it is drawn again in the normal mode.
"""
import math

import numpy as np

AX = "xyz"
PLANES = {2: "Front Plane", 1: "Top Plane", 0: "Right Plane"}       # normal z, y, x
PLANE_ORDER = {2: 0, 1: 1, 0: 2}                                    # position among the first planes


# STL frame -> build frame: the extrusion axis goes to z by a rotation (det +1), and the rotation
# (rad about x, y, z) that turns the built body back
BUILD_FRAME = {
    2: ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (0.0, 0.0, 0.0)),
    1: ([[1, 0, 0], [0, 0, -1], [0, 1, 0]], (-math.pi / 2, 0.0, 0.0)),     # +90 deg about x
    0: ([[0, 0, -1], [0, 1, 0], [1, 0, 0]], (0.0, math.pi / 2, 0.0)),      # -90 deg about y
}


def _f(x, nd=6):
    v = round(float(x), nd)
    return repr(0.0 if v == 0 else v)


def _is_hole(ft):
    return ft["op"] == "cut" and len(ft["loop"].edges) == 1 and ft["loop"].edges[0][0] == "circle"


def _steps(p25):
    """Features in order; consecutive round holes with the same range become one step."""
    feats = p25["features"]
    out, m = [], 0
    while m < len(feats):
        ft = feats[m]
        if _is_hole(ft):
            grp = [ft]
            while (m + len(grp) < len(feats) and _is_hole(feats[m + len(grp)])
                   and abs(feats[m + len(grp)]["z0"] - ft["z0"]) < 1e-9
                   and abs(feats[m + len(grp)]["z1"] - ft["z1"]) < 1e-9):
                grp.append(feats[m + len(grp)])
            out.append(dict(op="cut", z0=ft["z0"], z1=ft["z1"], loops=[g["loop"] for g in grp], holes=True))
            m += len(grp)
        else:
            out.append(dict(op=ft["op"], z0=ft["z0"], z1=ft["z1"], loops=[ft["loop"]], holes=False))
            m += 1
    return out


def expected_volumes(p25, steps):
    """Volume (mm^3) and centre of mass (mm) after each step, from the same construction in
    CadQuery. Volume alone does not fix a feature's position along the axis (a cut of the same
    section one slab lower removes the same volume); the centre of mass does."""
    import cadquery as cq
    from .prism25 import _wire
    k = p25["axis"]
    solid, vols = None, []
    for st in steps:
        for loop in st["loops"]:
            f = cq.Face.makeFromWires(_wire(loop, k, st["z0"]), [])
            D = [0.0, 0.0, 0.0]
            D[k] = st["z1"] - st["z0"]
            body = cq.Solid.extrudeLinear(f, cq.Vector(*D))
            if solid is None:
                solid = body
            elif st["op"] == "add":
                solid = solid.fuse(body)
            else:
                solid = solid.cut(body)
        c = solid.Center()                      # centre of mass
        vols.append((float(solid.Volume()), (c.x, c.y, c.z)))
    return vols


def cylinder_steps(p25, steps):
    """Off-axis bosses and holes, applied after the prisms (in the STL frame): for each, the
    point to pick its start face, the circle centre on that face, diameter, length, and the
    volume and centre of mass expected after it. None if a feature cannot be sketched on a face
    perpendicular to its axis (no face point was found)."""
    import cadquery as cq
    from .prism25 import _wire, cylinder_solid
    k = p25["axis"]
    solid = None
    for st in steps:
        for loop in st["loops"]:
            f = cq.Face.makeFromWires(_wire(loop, k, st["z0"]), [])
            D = [0.0, 0.0, 0.0]
            D[k] = st["z1"] - st["z0"]
            body = cq.Solid.extrudeLinear(f, cq.Vector(*D))
            solid = body if solid is None else (solid.fuse(body) if st["op"] == "add" else solid.cut(body))
    out = []
    for n, cy in enumerate(p25.get("cylinders") or [], 1):
        if cy.get("face_pt") is None:
            return None
        # pick point and circle centre on the model's planar face (the mesh lies within the fit
        # tolerance of it, the API picks on the model)
        cen, d = np.asarray(cy["c"]), np.asarray(cy["d"])
        pick, ctr = np.asarray(cy["face_pt"], float), cen + cy["tf"] * d
        best = None
        for fc in solid.Faces():
            if fc.geomType() != "PLANE":
                continue
            nrm = np.array(fc.normalAt(cq.Vector(*pick)).toTuple())
            if abs(nrm @ d) < np.cos(np.radians(1.0)):
                continue
            dist = cq.Vertex.makeVertex(*pick).distance(fc)
            if best is None or dist < best[0]:
                best = (dist, fc, nrm)
        if best is None or best[0] > 0.5:
            return None
        o = np.array(best[1].Center().toTuple())
        nrm = best[2]
        pick = pick - ((pick - o) @ nrm) * nrm
        ctr = ctr - ((ctr - o) @ nrm) * nrm
        body = cylinder_solid(cy)
        solid = solid.fuse(body).clean() if cy["op"] == "add" else solid.cut(body).clean()
        c = solid.Center()
        out.append(dict(op=cy["op"], pick=tuple(map(float, pick)), normal=tuple(map(float, nrm)),
                        centre=tuple(map(float, ctr)), dia=2 * cy["r"],
                        depth=cy["t1"] - cy["tf"], depth2=(cy["tf"] - cy["t0"]) if cy.get("spans") else 0.0,
                        volume=float(solid.Volume()),
                        centroid=(c.x, c.y, c.z), through=cy.get("through", False),
                        name=f"{'Boss' if cy['op'] == 'add' else 'Hole'} {n} dia {2 * cy['r']:g} "
                             f"along ({d[0]:.3g}, {d[1]:.3g}, {d[2]:.3g})"))
    return out


def _loop_code(loop, nd):
    """Python literal of a loop: list of edges in sketch-plane coordinates (a, b) in mm, or a
    pattern of one period repeated about a centre."""
    def P(p):
        return f"({_f(p[0], nd)}, {_f(p[1], nd)})"

    def edge(e):
        if e[0] == "line":
            return f"('line', {P(e[1])}, {P(e[2])})"
        if e[0] == "arc":
            return f"('arc', {P(e[1])}, {P(e[2])}, {P(e[3])})"
        if e[0] == "ellipse":
            return f"('ellipse', {P(e[1])}, {_f(e[2], nd)}, {_f(e[3], nd)}, {_f(e[4], 12)})"
        return f"('circle', {P(e[1])}, {_f(e[2], nd)})"
    if getattr(loop, "linear", None) is not None:
        p, N, v = loop.linear
        return (f"dict(count={N}, step={P(v)}, period=[\n"
                + "".join(f"            {edge(e)},\n" for e in loop.edges[:p])
                + "        ], rest=[\n"
                + "".join(f"            {edge(e)},\n" for e in loop.edges[p * N:]) + "        ])")
    if loop.pattern is not None:
        c, N, one = loop.pattern
        return ("dict(centre=" + P(c) + f", count={N}, period=[\n"
                + "".join(f"            {edge(e)},\n" for e in one) + "        ])")
    return "[\n" + "".join(f"            {edge(e)},\n" for e in loop.edges) + "        ]"


_RUNTIME = r'''
# ---- SolidWorks helpers (pywin32) -------------------------------------------------------------
import math
import sys
import time

import pythoncom
import win32com.client
from win32com.client import VARIANT

MM = 0.001                        # SolidWorks API lengths are in metres
NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)


def call0(obj, name):
    """Call a COM method that takes no arguments. With late binding, pywin32 may already run such
    a method when the attribute is read (it then returns the result, not a method), so the result
    is used as it comes."""
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


def settle(seconds=0.2):
    """Process pending Windows messages and wait briefly."""
    end = time.time() + seconds
    while time.time() < end:
        pythoncom.PumpWaitingMessages()
        time.sleep(0.02)


def make(create, *args):
    """Run a sketch-entity creation call; stop with the entity named if SolidWorks returns nothing."""
    r = create(*args)
    if r is None:
        raise RuntimeError(f"SolidWorks did not create the sketch entity {tuple(round(a, 6) for a in args)}")
    return r


HELPER_R = 0.5                    # mm, radius of the helper circle


SW_INPUT_DIM_VAL_ON_CREATE = 10          # swUserPreferenceToggle_e.swInputDimValOnCreate


def first_point(loops):
    loop = loops[0]
    e = expand(loop)[0]
    return e[1]


def model_point(a, b):
    """Point (a, b) of the sketch plane in model coordinates (mm), level 0 along the axis."""
    p = [0.0, 0.0, 0.0]
    p[AXIS_I], p[AXIS_J] = a, b
    return p


def expand(loop):
    """Edges of a loop; a pattern is one period turned COUNT times about its centre."""
    if isinstance(loop, list):
        return loop
    if "step" in loop:                     # linear pattern: COUNT copies of the period, then the rest
        (dx, dy), n = loop["step"], loop["count"]
        out = [(e[0],) + tuple((p[0] + m * dx, p[1] + m * dy) for p in e[1:])
               for m in range(n) for e in loop["period"]] + [tuple(e) for e in loop["rest"]]
        for m in range(1, len(out)):       # identical shared vertices, so the outline is closed
            out[m] = (out[m][0], out[m - 1][-1]) + tuple(out[m][2:])
        out[-1] = tuple(out[-1][:-1]) + (out[0][1],)
        return out
    (cx, cy), n = loop["centre"], loop["count"]
    out = []
    for m in range(n):
        t = 2 * math.pi * m / n
        c, s = math.cos(t), math.sin(t)
        rot = lambda p: (cx + c * (p[0] - cx) - s * (p[1] - cy), cy + s * (p[0] - cx) + c * (p[1] - cy))
        for e in loop["period"]:
            out.append((e[0],) + tuple(rot(p) for p in e[1:]))
    for m in range(1, len(out)):           # identical shared vertices, so the outline is closed
        out[m] = (out[m][0], out[m - 1][-1]) + tuple(out[m][2:])
    out[-1] = tuple(out[-1][:-1]) + (out[0][1],)
    return out


class SW:
    def __init__(self):
        self.app = win32com.client.Dispatch("SldWorks.Application")
        self.app.Visible = True
        self.dims_ok, self.dims_failed = 0, []
        self.rels_ok, self.rels_failed, self.status = 0, [], []
        try:                                     # no Modify dialog when a dimension is created
            self.input_dim_pref = self.app.GetUserPreferenceToggle(SW_INPUT_DIM_VAL_ON_CREATE)
            self.app.SetUserPreferenceToggle(SW_INPUT_DIM_VAL_ON_CREATE, False)
        except Exception:
            self.input_dim_pref = None
        self.turned_back = False
        self.ents, self._picked = [], []
        self.linkable = []                       # (R or D, value mm, dimension name) for linking
        self.use_patterns = True
        self.model = call0(self.app, "NewPart")
        if self.model is None:                   # no default part template set: use it explicitly
            template = self.app.GetUserPreferenceStringValue(8)      # swDefaultTemplatePart
            self.model = self.app.NewDocument(template, 0, 0, 0)

    def selected(self):
        return self.model.SelectionManager.GetSelectedObjectCount2(-1)

    def base_plane(self):
        """Select the base plane, checking that the selection took: by name, then the feature of
        that name, then the n-th reference plane of the template."""
        tried = []
        self.model.ClearSelection2(True)
        r = self.model.Extension.SelectByID2(PLANE_NAME, "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
        tried.append(f"SelectByID2('{PLANE_NAME}') returned {r}")
        if self.selected() > 0:
            return
        f = self.model.FeatureByName(PLANE_NAME)
        tried.append(f"FeatureByName('{PLANE_NAME}') returned {'a feature' if f is not None else None}")
        if f is not None:
            f.Select2(False, 0)
            if self.selected() > 0:
                return
        f, n, names = call0(self.model, "FirstFeature"), 0, []
        while f is not None:
            if call0(f, "GetTypeName2") == "RefPlane":
                names.append(f.Name)
                if n == PLANE_INDEX:
                    f.Select2(False, 0)
                    if self.selected() > 0:
                        return
                n += 1
            f = call0(f, "GetNextFeature")
        tried.append(f"reference planes found: {names}")
        raise RuntimeError("could not select the base plane; " + "; ".join(tried))

    def to_sketch(self, a, b):
        """Model point (mm) -> Front Plane sketch coordinates (m). Over pywin32, SolidWorks 2025
        put API sketches on the Front Plane whichever plane was selected, so every part is built
        on the Front Plane with its extrusion axis along z: ROT turns the STL frame into the build
        frame (identity for parts extruded along z)."""
        x = [v * MM for v in model_point(a, b)]
        X = [sum(ROT[r][c] * x[c] for c in range(3)) for r in range(3)]
        return X[0], X[1]

    def sketch(self, loops, add_to_db=True, plans=None):
        self.base_plane()
        sm = self.model.SketchManager
        sm.InsertSketch(True)
        settle()                             # let SolidWorks finish entering the sketch
        sm.AddToDB = add_to_db               # True: no snapping or automatic relations while drawing
        if not add_to_db:                    # normal mode: no inferencing (snapped a 5 mm slot shut)
            try:
                sm.InferenceMode = False
            except Exception as ex:
                print(f"  could not turn sketch inferencing off: {type(ex).__name__}")
        sm.DisplayWhenAdded = not add_to_db
        # Over pywin32, SolidWorks 2025 creates no line in a sketch until another entity exists;
        # a small helper circle is drawn first and removed (or made construction) at the end.
        p0 = self.to_sketch(*HELPER_AB)          # away from the part, so selections do not hit it
        helper = sm.CreateCircleByRadius(p0[0], p0[1], 0.0, HELPER_R * MM)    # may return None
        self.ents = []                           # (kind, sketch geometry (m), entity) as drawn
        patterns = []                            # (entities of the period, count, step) to repeat
        circles = []                             # (recorded period entities, loop) to turn
        for loop in loops:
            linear = isinstance(loop, dict) and "step" in loop
            circular = isinstance(loop, dict) and "centre" in loop and self.use_patterns
            first = len(self.ents)
            todo = expand(loop)
            if linear:                           # the period and the rest; the copies by pattern
                todo = todo[:len(loop["period"])] + todo[len(loop["period"]) * loop["count"]:]
            elif circular:                       # the period; the copies by a circular pattern
                todo = todo[:len(loop["period"])]
            for e in todo:
                if e[0] == "line":
                    (x1, y1), (x2, y2) = self.to_sketch(*e[1]), self.to_sketch(*e[2])
                    g = make(sm.CreateLine, x1, y1, 0.0, x2, y2, 0.0)
                    self.ents.append(("line", ((x1, y1), (x2, y2)), g))
                elif e[0] == "circle":
                    (xc, yc) = self.to_sketch(*e[1])
                    g = make(sm.CreateCircleByRadius, xc, yc, 0.0, e[2] * MM)
                    self.ents.append(("circle", ((xc, yc), e[2] * MM), g))
                elif e[0] == "ellipse":
                    (cx, cy), a, b, ang = e[1], e[2], e[3], e[4]
                    u = (math.cos(ang), math.sin(ang))
                    pc = self.to_sketch(cx, cy)
                    pa = self.to_sketch(cx + a * u[0], cy + a * u[1])
                    pb = self.to_sketch(cx - b * u[1], cy + b * u[0])
                    g = make(sm.CreateEllipse, pc[0], pc[1], 0.0, pa[0], pa[1], 0.0, pb[0], pb[1], 0.0)
                    self.ents.append(("ellipse", (pc, pa, pb), g))
                else:
                    s, m, t = (self.to_sketch(*p) for p in e[1:])
                    xc, yc = circle_centre(s, m, t)
                    turn = (m[0] - s[0]) * (t[1] - s[1]) - (m[1] - s[1]) * (t[0] - s[0])
                    g = make(sm.CreateArc, xc, yc, 0.0, s[0], s[1], 0.0, t[0], t[1], 0.0, 1 if turn > 0 else -1)
                    self.ents.append(("arc", (s, t, (xc, yc), math.dist(s, (xc, yc)), m), g))
            if linear:
                patterns.append(([g for _, _, g in self.ents[first:first + len(loop["period"])]],
                                 loop["count"], loop["step"]))
            elif circular:
                circles.append((self.ents[first:first + len(loop["period"])], loop))
        for seed, count, step in patterns:
            self.sketch_pattern(seed, count, step)
        for ents, loop in circles:
            self.circular_pattern(ents, loop)
        sm.AddToDB = False
        sm.DisplayWhenAdded = True
        if helper is not None:
            self.remove_helper(helper)
        self.dimension_holes(loops)              # every full circle: diameter and position
        self.dimension_ellipses()                # centre, axis direction, semi-axes
        for plan in plans or []:
            self.apply_plan(plan)
        self.region_report(add_to_db)
        sm.InsertSketch(True)                # close the sketch
        if not add_to_db:
            try:
                sm.InferenceMode = True
            except Exception:
                pass
        self.last_loops = loops
        self.report_status()
        return self.model.FeatureByPositionReverse(0)

    def region_report(self, add_to_db):
        """Diagnostic: closed regions SolidWorks finds in the open sketch (0 means it sees the
        outline as open, for example endpoints that are not merged)."""
        try:
            sk = call0(self.model, "GetActiveSketch2")
            n = call0(sk, "GetSketchRegionCount") if sk is not None else None
            print(f"  sketch ({'AddToDB' if add_to_db else 'normal mode'}): {n} closed region(s)")
        except Exception as ex:
            print(f"  sketch region count not available: {type(ex).__name__}")

    def dump_sketch(self, sk, everything=False):
        """Diagnostic after a failed check: sketch segments whose endpoints are not at the
        intended vertices (sketch coordinates, mm)."""
        try:
            want = []
            for loop in self.last_loops:
                for e in expand(loop):
                    want += [self.to_sketch(*q) for q in e[1:] if isinstance(q, (tuple, list)) and len(q) == 2]
            W = [(x / MM, y / MM) for x, y in want]
            segs = call0(call0(sk, "GetSpecificFeature2"), "GetSketchSegments") or []
            bad = 0
            for g in segs:
                if call0(g, "ConstructionGeometry"):
                    continue
                if call0(g, "GetType") != 0:           # lines only
                    if everything:
                        print(f"  segment of type {call0(g, 'GetType')}")
                    continue
                pts = [call0(g, "GetStartPoint2"), call0(g, "GetEndPoint2")]
                P = [(p.X / MM, p.Y / MM) for p in pts]
                off = [min(math.dist(q, w) for w in W) for q in P]
                if everything:
                    print(f"  line ({P[0][0]:.5f}, {P[0][1]:.5f}) -> ({P[1][0]:.5f}, {P[1][1]:.5f})"
                          + ("   MOVED" if max(off) > 1e-4 else ""))
                if max(off) > 1e-4:
                    bad += 1
                    print(f"  moved line: ({P[0][0]:.4f}, {P[0][1]:.4f}) -> ({P[1][0]:.4f}, {P[1][1]:.4f}), "
                          f"endpoints {off[0]:.4f} and {off[1]:.4f} mm from the nearest intended vertex")
            print(f"  sketch read back: {len(segs)} segment(s), {bad} line(s) with moved endpoints")
        except Exception as ex:
            print(f"  sketch read-back failed: {type(ex).__name__}: {ex}")

    def restore_preferences(self):
        if self.input_dim_pref is not None:
            try:
                self.app.SetUserPreferenceToggle(SW_INPUT_DIM_VAL_ON_CREATE, self.input_dim_pref)
            except Exception:
                pass

    def dimension_holes(self, loops):
        """Diameter of each hole, and its centre's horizontal and vertical distance from the
        origin. Entities are selected by coordinates (empty name), as the API help advises for
        dimensioning. A failure is reported, not fatal: the geometry is already complete."""
        ext = self.model.Extension
        first = []                               # (radius, entity) of the first circle of each size
        for loop in loops:
            e = expand(loop)[0]
            if e[0] != "circle":
                continue
            cx, cy = self.to_sketch(*e[1])
            r = e[2] * MM
            tag = f"hole at ({e[1][0]:g}, {e[1][1]:g})"
            g = next((t[2] for t in self.ents if t[0] == "circle" and math.dist(t[1][0], (cx, cy)) < 1e-9
                      and abs(t[1][1] - r) < 1e-12), None)
            same = next((h for rr, h in first if abs(rr - r) < 1e-9 * max(1.0, r / MM)), None)
            try:
                if same is not None and g is not None:   # same size as an earlier circle: Equal
                    self.model.ClearSelection2(True)
                    if same.Select4(False, NOTHING) and g.Select4(True, NOTHING):
                        self.model.SketchAddConstraints("sgSAMELENGTH")
                        self.rels_ok += 1
                    else:
                        self.rels_failed.append(tag + " equal diameter (selection)")
                else:
                    if g is not None:
                        first.append((r, g))
                    self.model.ClearSelection2(True)
                    ok = ext.SelectByID2("", "SKETCHSEGMENT", cx + r, cy, 0.0, False, 0, NOTHING, 0)
                    d = self.model.AddDimension2(cx + 1.6 * r, cy + 1.6 * r, 0.0) if ok else None
                    self._count(d, tag + " diameter")
                for kind, tx, ty in (("H", cx / 2, cy - 3 * r), ("V", cx - 3 * r, cy / 2)):
                    if abs(cx if kind == "H" else cy) < 1e-9:
                        continue                 # centre on the origin's line: no distance needed
                    self.model.ClearSelection2(True)
                    ok1 = ext.SelectByID2("", "SKETCHPOINT", cx, cy, 0.0, False, 0, NOTHING, 0)
                    ok2 = ext.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0, True, 0, NOTHING, 0)
                    d = None
                    if ok1 and ok2:
                        d = (self.model.AddHorizontalDimension2(tx, ty, 0.0) if kind == "H"
                             else self.model.AddVerticalDimension2(tx, ty, 0.0))
                    self._count(d, tag + (" x" if kind == "H" else " y") + " position")
            except Exception as ex:
                self.dims_failed.append(f"{tag}: {type(ex).__name__}: {ex}")
        self.model.ClearSelection2(True)

    def dimension_ellipses(self):
        """Each ellipse fully defined: centre on the origin (Coincident) or dimensioned from it,
        the major axis Horizontal or Vertical (relation between the centre and the major-axis
        point), and the semi-axes as distances from the centre to the axis points. The points
        are the ellipse's own (GetCenterPoint2, GetMajorPoint2, GetMinorPoint2), selected
        directly. A failure is reported, not fatal."""
        ext = self.model.Extension
        for kind, (pc, pa, pb), g in [t for t in self.ents if t[0] == "ellipse"]:
            tag = f"ellipse at ({pc[0] / MM:g}, {pc[1] / MM:g})"
            try:
                c, ma, mi = (call0(g, n) for n in ("GetCenterPoint2", "GetMajorPoint2", "GetMinorPoint2"))

                def pick(*pts):
                    self.model.ClearSelection2(True)
                    ok = True
                    for n, p in enumerate(pts):
                        if p == "origin":
                            ok = ok and ext.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0,
                                                        n > 0, 0, NOTHING, 0)
                        else:
                            ok = ok and p.Select4(n > 0, NOTHING)
                    return ok

                def rel(name, what, *pts):
                    if pick(*pts):
                        self.model.SketchAddConstraints(name)
                        self.rels_ok += 1
                    else:
                        self.rels_failed.append(tag + " " + what)
                if math.hypot(*pc) < 1e-9:
                    rel("sgCOINCIDENT", "centre on the origin", c, "origin")
                else:
                    for k_, (tx, ty) in (("H", (pc[0] / 2, pc[1] - 0.005)), ("V", (pc[0] - 0.005, pc[1] / 2))):
                        if abs(pc[0] if k_ == "H" else pc[1]) < 1e-9:
                            continue
                        d = None
                        if pick(c, "origin"):
                            d = (self.model.AddHorizontalDimension2(tx, ty, 0.0) if k_ == "H"
                                 else self.model.AddVerticalDimension2(tx, ty, 0.0))
                        self._count(d, tag + (" x" if k_ == "H" else " y") + " of the centre")
                dx, dy = pa[0] - pc[0], pa[1] - pc[1]
                if abs(dy) < 1e-9 * max(1.0, abs(dx)):
                    rel("sgHORIZONTALPOINTS2D", "major axis horizontal", c, ma)
                elif abs(dx) < 1e-9 * max(1.0, abs(dy)):
                    rel("sgVERTICALPOINTS2D", "major axis vertical", c, ma)
                for p, q, what in ((ma, pa, "major"), (mi, pb, "minor")):
                    d = self.model.AddDimension2((pc[0] + q[0]) / 2 + 0.003, (pc[1] + q[1]) / 2 + 0.003, 0.0) \
                        if pick(c, p) else None
                    self._count(d, f"{tag} {what} semi-axis")
            except Exception as ex:
                self.dims_failed.append(f"{tag}: {type(ex).__name__}: {ex}")
        self.model.ClearSelection2(True)

    # ---- outline relations and dimensions -------------------------------------------------------
    def _select(self, typ, a, b, append=False):
        """Select a sketch segment or point at model point (a, b) (mm) of the build frame. The
        entity objects recorded while drawing are selected directly (Select4): SelectByID2 picks
        whatever lies within its pick radius, which depends on the zoom, and inside a 5 mm slot
        it picked the neighbouring line, so relations and dimensions went to the wrong entity and
        the solver made edges cross. SelectByID2 remains the fall-back."""
        x, y = self.to_sketch(a, b)
        if not append:
            self._picked = []
        try:
            obj = self._entity_at(typ, x, y)
            if obj is not None and obj.Select4(append, NOTHING):
                self._picked.append(obj)
                return True
        except Exception as ex:
            print(f"  direct selection failed ({type(ex).__name__}); selecting by coordinates")
        return self.model.Extension.SelectByID2("", typ, x, y, 0.0, append, 0, NOTHING, 0)

    def _entity_at(self, typ, x, y, tol=1e-6):
        """Recorded segment through (x, y), or an endpoint object at (x, y) not yet picked in this
        selection (m)."""
        P = (x, y)
        if typ == "SKETCHSEGMENT":
            best = None
            for kind, geo, g in self.ents:
                if kind in ("line", "cline"):
                    (ax, ay), (bx, by) = geo
                    dx, dy = bx - ax, by - ay
                    t = max(0.0, min(1.0, ((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy)))
                    d = math.dist(P, (ax + t * dx, ay + t * dy))
                    if kind == "cline":          # a centerline only where no outline edge is
                        d += 0.5 * tol
                elif kind == "arc":
                    d = abs(math.dist(P, geo[2]) - geo[3])
                elif kind == "circle":
                    d = abs(math.dist(P, geo[0]) - geo[1])
                else:
                    continue
                if d < tol and (best is None or d < best[0]):
                    best = (d, g)
            return None if best is None else best[1]
        if typ == "SKETCHPOINT":
            for kind, geo, g in self.ents:
                if kind not in ("line", "arc"):
                    continue
                if min(math.dist(P, geo[0]), math.dist(P, geo[1])) > tol:
                    continue
                for name in ("GetStartPoint2", "GetEndPoint2"):   # an arc may store its ends swapped
                    pt = call0(g, name)
                    if (pt is not None and math.dist(P, (pt.X, pt.Y)) < tol
                            and not any(pt == q for q in self._picked)):
                        return pt
        return None

    def _relation(self, rel, picks, what):
        try:
            self.model.ClearSelection2(True)
            ok = True
            for k, (t, p) in enumerate(picks):
                if t in ("EXTSKETCHPOINT_ORIGIN", "ORIGIN"):
                    ok = ok and self.model.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0,
                                                                  True, 0, NOTHING, 0)
                else:
                    ok = ok and self._select(t, *p, append=k > 0)
            if ok:
                self.model.SketchAddConstraints(rel)
                self.rels_ok += 1
            else:
                self.rels_failed.append(what + " (selection)")
        except Exception as ex:
            self.rels_failed.append(f"{what}: {type(ex).__name__}: {ex}")

    def _dimension(self, picks, text, kind, what):
        try:
            self.model.ClearSelection2(True)
            ok = True
            for k, (t, p) in enumerate(picks):
                if t == "ORIGIN":
                    ok = ok and self.model.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0,
                                                                  True, 0, NOTHING, 0)
                else:
                    ok = ok and self._select(t, *p, append=k > 0)
            d = None
            if ok:
                tx, ty = self.to_sketch(*text)
                d = {"any": self.model.AddDimension2, "H": self.model.AddHorizontalDimension2,
                     "V": self.model.AddVerticalDimension2}[kind](tx, ty, 0.0)
            self._count(d, what)
        except Exception as ex:
            self.dims_failed.append(f"{what}: {type(ex).__name__}: {ex}")

    def centerline(self, ends, what):
        """Construction centerline through the sketch origin along a symmetry axis (ends in model
        mm), fixed by a Horizontal or Vertical relation, the origin at its middle, and its length."""
        (x1, y1), (x2, y2) = self.to_sketch(*ends[0]), self.to_sketch(*ends[1])
        sm = self.model.SketchManager
        g = sm.CreateCenterLine(x1, y1, 0.0, x2, y2, 0.0)
        if g is None:
            self.rels_failed.append(what + " (not created)")
            return
        self.ents.append(("cline", ((x1, y1), (x2, y2)), g))
        rel = "sgVERTICAL2D" if abs(x2 - x1) < abs(y2 - y1) else "sgHORIZONTAL2D"
        try:
            self.model.ClearSelection2(True)
            if g.Select4(False, NOTHING):
                self.model.SketchAddConstraints(rel)
                self.rels_ok += 1
            self.model.ClearSelection2(True)
            ok = self.model.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0, False, 0,
                                                  NOTHING, 0) and g.Select4(True, NOTHING)
            if ok:
                self.model.SketchAddConstraints("sgATMIDDLE")
                self.rels_ok += 1
            else:
                self.rels_failed.append(what + " middle on the origin (selection)")
            self.model.ClearSelection2(True)
            d = self.model.AddDimension2(x2 + 0.002, y2 + 0.002, 0.0) if g.Select4(False, NOTHING) else None
            self._count(d, what + " length")
        except Exception as ex:
            self.rels_failed.append(f"{what}: {type(ex).__name__}: {ex}")
        self.model.ClearSelection2(True)

    def circular_pattern(self, ents, loop):
        """Turn the period's entities COUNT times about the loop centre with a sketch circular
        pattern (SketchManager.CreateCircularSketchStepAndRepeat). SolidWorks 2025 (sw_probe10-12)
        puts the pattern centre at the centre of the seed's bounding box plus ArcRadius at
        ArcAngle (radians, 0..2 pi); the copies then share end points. If the call is refused, the
        copies are drawn one by one."""
        xs, ys = [], []
        for kind, geo, g in ents:
            if kind == "line":
                xs += [geo[0][0], geo[1][0]]
                ys += [geo[0][1], geo[1][1]]
            else:                                # arc: ends plus the axis extremes it sweeps
                s, t, c, r, m = geo
                xs += [s[0], t[0]]
                ys += [s[1], t[1]]
                a0, a1 = math.atan2(s[1] - c[1], s[0] - c[0]), math.atan2(t[1] - c[1], t[0] - c[0])
                am = math.atan2(m[1] - c[1], m[0] - c[0])
                sw_ = (a1 - a0) % (2 * math.pi)
                ccw = (am - a0) % (2 * math.pi) < sw_
                for k in range(4):
                    ang = k * math.pi / 2
                    d = (ang - a0) % (2 * math.pi) if ccw else (a0 - ang) % (2 * math.pi)
                    if d < (sw_ if ccw else 2 * math.pi - sw_):
                        xs.append(c[0] + r * math.cos(ang))
                        ys.append(c[1] + r * math.sin(ang))
        bc = ((max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2)
        cx, cy = self.to_sketch(*loop["centre"])
        R, A = math.hypot(cx - bc[0], cy - bc[1]), math.atan2(cy - bc[1], cx - bc[0]) % (2 * math.pi)
        n = loop["count"]
        p0, p1 = self.to_sketch(0.0, 0.0), self.to_sketch(1.0, 0.0)
        q1 = self.to_sketch(0.0, 1.0)
        turn = 1.0 if (p1[0] - p0[0]) * (q1[1] - p0[1]) - (p1[1] - p0[1]) * (q1[0] - p0[0]) > 0 else -1.0
        self.model.ClearSelection2(True)
        ok = all(g.Select4(True, NOTHING) for _, _, g in ents)
        r = self.model.SketchManager.CreateCircularSketchStepAndRepeat(
            R, A, n, turn * 2 * math.pi / n, True, "", False, False, False) if ok else False
        self.model.ClearSelection2(True)
        if r:
            print(f"  circular pattern: {len(ents)} entities x {n}")
            return
        print(f"  circular pattern refused; drawing the {n - 1} copies")
        sm = self.model.SketchManager
        for e in expand(loop)[len(loop["period"]):]:
            if e[0] == "line":
                (x1, y1), (x2, y2) = self.to_sketch(*e[1]), self.to_sketch(*e[2])
                make(sm.CreateLine, x1, y1, 0.0, x2, y2, 0.0)
            else:
                s_, m_, t_ = (self.to_sketch(*q) for q in e[1:])
                xc, yc = circle_centre(s_, m_, t_)
                tn = (m_[0] - s_[0]) * (t_[1] - s_[1]) - (m_[1] - s_[1]) * (t_[0] - s_[0])
                make(sm.CreateArc, xc, yc, 0.0, s_[0], s_[1], 0.0, t_[0], t_[1], 0.0, 1 if tn > 0 else -1)

    def sketch_pattern(self, seed, count, step):
        """Repeat the period's entities count times by step (model mm) with a sketch linear
        pattern (SketchManager.CreateLinearSketchStepAndRepeat). In SolidWorks 2025 the copies
        share their end points with each other and with the rest of the outline (sw_probe8)."""
        p0, p1 = self.to_sketch(0.0, 0.0), self.to_sketch(*step)
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        self.model.ClearSelection2(True)
        ok = all(g.Select4(True, NOTHING) for g in seed)
        r = self.model.SketchManager.CreateLinearSketchStepAndRepeat(
            count, 1, math.hypot(dx, dy), 0.0, math.atan2(dy, dx), math.pi / 2, "",
            False, False, False, False, False) if ok else False
        self.model.ClearSelection2(True)
        if not r:
            raise RuntimeError(f"sketch linear pattern of {len(seed)} entities x {count} was not created")
        print(f"  sketch pattern: {len(seed)} entities x {count}, step {math.hypot(dx, dy) / MM:g} mm")

    def apply_plan(self, plan):
        """Relations and dimensions planned by stl2cad for one outline, in order."""
        for op in plan:
            if op[0] == "centerline":
                self.centerline(op[2], op[3])
                continue
            if op[0] == "rel":
                _, rel, picks, what = op
                self._relation(rel, picks, what)
            else:
                _, kind, picks, text, what = op
                self._dimension(picks, text, kind, what)
        self.model.ClearSelection2(True)

    def report_status(self):
        """Constrained status of the sketch just closed, if SolidWorks returns it."""
        try:
            f = self.model.FeatureByPositionReverse(0)
            sk = call0(f, "GetSpecificFeature2")
            st = call0(sk, "GetConstrainedStatus")
            self.status.append({1: "unknown", 2: "under defined", 3: "fully defined", 4: "over defined",
                                5: "no solution", 6: "invalid solution"}.get(st, str(st)))
        except Exception as ex:
            self.status.append(f"not available ({type(ex).__name__})")

    def _count(self, d, what):
        if d is None:
            self.dims_failed.append(what)
            return
        self.dims_ok += 1
        if "radius" in what or "diameter" in what:   # candidates for linking across sketches
            try:
                dd = d.GetDimension2(0)
                name = "@".join(str(dd.FullName).split("@")[:2])         # D1@Sketch2
                self.linkable.append(("R" if "radius" in what else "D", float(dd.SystemValue) / MM, name))
            except Exception:
                pass

    def link_dimensions(self):
        """Equal radii and equal diameters in different sketches driven by one global variable each
        (Equation Manager, Add2; sw_probe11: the dimensions follow after a rebuild). A failure is
        reported, not fatal."""
        groups = {}
        for kind, val, name in self.linkable:
            groups.setdefault((kind, round(val, 6)), []).append(name)
        groups = {k: sorted(set(v)) for k, v in groups.items() if len(set(v)) > 1}
        if not groups:
            return "linked dimensions: none (no equal radius or diameter in two sketches)"
        try:
            em = call0(self.model, "GetEquationMgr")
            n = 0
            for (kind, val), names in sorted(groups.items()):
                var = f"{'radius' if kind == 'R' else 'dia'}_{val:g}".replace(".", "_").replace("-", "m")
                em.Add2(-1, f'"{var}" = {val:g}', True)
                for nm in names:
                    em.Add2(-1, f'"{nm}" = "{var}"', True)
                    n += 1
            call0(self.model, "EditRebuild3")
            return f"linked dimensions: {n} in {len(groups)} global variable(s)"
        except Exception as ex:
            return f"linked dimensions: failed ({type(ex).__name__}: {ex})"

    def remove_helper(self, seg):
        """Delete the helper circle; if that fails, make it construction geometry, which
        extrusions ignore."""
        try:
            self.model.ClearSelection2(True)
            if seg.Select4(False, NOTHING) and self.selected() > 0:
                self.model.EditDelete()
                if self.selected() == 0:
                    return
        except Exception:
            pass
        seg.ConstructionGeometry = True

    def rotate_back(self):
        """Turn the body from the build frame back to the STL's orientation with a Move/Copy
        Body feature (body selected with Mark 1). Which angle argument of InsertMoveCopyBody2
        turns about which axis is not assumed (in SolidWorks 2025 the first angle turned about z):
        a quarter turn of each sign in each argument is tried, and the one that brings the centre
        of mass to its place in the STL frame is kept. Otherwise the part stays in the build frame."""
        self.turned_back = True
        if all(abs(ROT[r][c] - (r == c)) < 1e-12 for r in range(3) for c in range(3)):
            return "part is in the STL's orientation"
        for slot in range(3):
            for sign in (1, -1):
                ang = [0.0, 0.0, 0.0]
                ang[slot] = sign * math.pi / 2
                f = None
                try:
                    bodies = self.model.GetBodies2(0, True)          # swSolidBody
                    self.model.ClearSelection2(True)
                    sd = call0(self.model.SelectionManager, "CreateSelectData")
                    sd.Mark = 1
                    bodies[0].Select2(False, sd)
                    f = self.model.FeatureManager.InsertMoveCopyBody2(
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ang[0], ang[1], ang[2], False, 1)
                except Exception as e:
                    print(f"  rotation attempt failed: {type(e).__name__}: {e}")
                if f is None:
                    continue
                v, c = self.mass()
                if math.dist(c, FINAL_CENTROID) <= 1e-4 * SIZE + 1e-4:
                    f.Name = "Rotate to STL orientation"
                    return f"rotated back to the STL's orientation (Move/Copy Body, angle argument {slot + 1}, {'+' if sign > 0 else '-'}90 deg)"
                self.model.ClearSelection2(True)
                f.Select2(False, 0)
                self.model.Extension.DeleteSelection2(0)
        self.turned_back = False
        return ("NOTE: part left in the build frame (extrusion axis along z); no quarter turn brought it "
                "to the STL's orientation")

    def in_model_frame(self, p):
        """STL-frame point (mm) in the part's frame: the STL frame once the body is turned back,
        else the build frame."""
        if self.turned_back:
            return tuple(p)
        return tuple(sum(ROT[r][c] * p[c] for c in range(3)) for r in range(3))

    def select_face(self, cy):
        """Select the planar face the feature starts on. SelectByID2 with a point picks the first
        face along the line of sight through the point (seen in SolidWorks 2025: a hole meant for
        a side face was sketched on the top face), so the face is picked with a ray fired at the
        point from just outside the face along its inward normal (SelectByRay), and the selected
        face's normal is checked; SelectByID2 is the fall-back."""
        ext = self.model.Extension
        p = [v * MM for v in self.in_model_frame(cy["pick"])]
        n = list(self.in_model_frame(cy["normal"]))
        tried = []
        for how in ("ray", "id"):
            self.model.ClearSelection2(True)
            try:
                if how == "ray":
                    q = [p[i] + 0.001 * n[i] for i in range(3)]          # 1 mm outside the face
                    ok = ext.SelectByRay(q[0], q[1], q[2], -n[0], -n[1], -n[2], 0.0001, 2, False, 0, 0)
                else:
                    ok = ext.SelectByID2("", "FACE", p[0], p[1], p[2], False, 0, NOTHING, 0)
            except Exception as ex:
                tried.append(f"{how}: {type(ex).__name__}")
                continue
            if not ok:
                tried.append(f"{how}: nothing selected")
                continue
            try:
                fc = self.model.SelectionManager.GetSelectedObject6(1, -1)
                fn = list(call0(fc, "Normal"))
                if sum(a * b for a, b in zip(fn, n)) > math.cos(math.radians(1.0)):
                    return
                tried.append(f"{how}: face normal {tuple(round(v, 3) for v in fn)}")
            except Exception as ex:                 # normal not readable: trust the ray pick
                if how == "ray":
                    return
                tried.append(f"{how}: {type(ex).__name__}")
        raise RuntimeError(f"{cy['name']}: could not select the start face at {cy['pick']} ({'; '.join(tried)})")

    def face_feature(self, cy):
        """Boss or hole along the normal of a planar face: the face is picked at a point next to
        the circle, the circle is placed through the face sketch's ModelToSketchTransform (in
        SolidWorks 2025 a face sketch takes its own axes: on the x = const face, sketch x is -z
        and sketch y is +y), dimensioned, and extruded blind; both directions are tried against
        the expected volume and centre of mass."""
        ext, sm, fm = self.model.Extension, self.model.SketchManager, self.model.FeatureManager
        self.select_face(cy)
        sm.InsertSketch(True)
        settle()
        sk = call0(self.model, "GetActiveSketch2")
        if sk is None:
            sk = sm.ActiveSketch
        if sk is None:
            raise RuntimeError(f"{cy['name']}: no active sketch on the face")
        A = list(call0(sk, "ModelToSketchTransform").ArrayData)
        q = [v * MM for v in self.in_model_frame(cy["centre"])]
        cx, cyy = (sum(q[r] * A[3 * r + i] for r in range(3)) + A[9 + i] for i in range(2))
        r = cy["dia"] / 2 * MM
        sm.CreateCircleByRadius(cx, cyy, 0.0, r)
        tag = cy["name"]
        try:
            self.model.ClearSelection2(True)
            ok = ext.SelectByID2("", "SKETCHSEGMENT", cx + r, cyy, 0.0, False, 0, NOTHING, 0)
            self._count(self.model.AddDimension2(cx + 1.6 * r, cyy + 1.6 * r, 0.0) if ok else None,
                        tag + " diameter")
            for kind, tx, ty in (("H", cx / 2, cyy - 3 * r), ("V", cx - 3 * r, cyy / 2)):
                if abs(cx if kind == "H" else cyy) < 1e-9:
                    continue
                self.model.ClearSelection2(True)
                ok1 = ext.SelectByID2("", "SKETCHPOINT", cx, cyy, 0.0, False, 0, NOTHING, 0)
                ok2 = ext.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0, True, 0, NOTHING, 0)
                d = None
                if ok1 and ok2:
                    d = (self.model.AddHorizontalDimension2(tx, ty, 0.0) if kind == "H"
                         else self.model.AddVerticalDimension2(tx, ty, 0.0))
                self._count(d, tag + (" x" if kind == "H" else " y") + " position on the face")
        except Exception as ex:
            self.dims_failed.append(f"{tag}: {type(ex).__name__}: {ex}")
        sm.InsertSketch(True)
        self.report_status()
        sk = self.model.FeatureByPositionReverse(0)
        expected, centroid = cy["volume"], self.in_model_frame(cy["centroid"])
        h = cy["depth"] * MM
        for rev in (False, True):
            self.model.ClearSelection2(True)
            sk.Select2(False, 0)
            if cy["op"] == "add":
                h2 = cy.get("depth2", 0.0) * MM          # a boss through the body: both directions
                f = fm.FeatureExtrusion3(h2 <= 0, False, rev, 0, 0, h, h2, False, False, False, False,
                                         0.0, 0.0, False, False, False, False, True, True, True,
                                         0, 0.0, False)
            else:
                f = fm.FeatureCut4(True, False, rev, 0, 0, h, 0.0, False, False, False, False,
                                   0.0, 0.0, False, False, False, False, False, True, True,
                                   False, False, False, 0, 0.0, False, False)
            if f is None:
                print(f"  try reverse={rev}: SolidWorks returned no feature")
                continue
            v, c = self.mass()
            print(f"  try reverse={rev}: volume {v:.3f} mm^3, centre of mass ({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}) mm")
            if abs(v - expected) <= 1e-4 * expected + 1e-3 and math.dist(c, centroid) <= 1e-5 * SIZE + 1e-4:
                f.Name = cy["name"]
                return v
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(0)
        raise RuntimeError(f"{cy['name']}: neither direction gives the expected volume {expected:.3f} mm^3")

    def mass(self):
        """Volume (mm^3) and centre of mass (mm) of the part."""
        self.model.ForceRebuild3(False)
        mp = call0(self.model.Extension, "CreateMassProperty")
        if mp is None:
            return 0.0, (0.0, 0.0, 0.0)
        return mp.Volume / MM ** 3, tuple(c / MM for c in mp.CenterOfMass)

    def feature(self, sk, op, z0, z1, expected, centroid, name):
        """Extruded boss or cut of sketch sk from level z0 to z1. Each combination of direction
        and start-offset side is tried until the part volume and centre of mass equal the
        expected ones."""
        fm = self.model.FeatureManager
        h, off = (z1 - z0) * MM, abs(z0) * MM
        t0 = 3 if abs(z0) > 1e-9 else 0      # swStartOffset, or swStartSketchPlane
        # SolidWorks 2025 (observed): a boss goes along the plane normal, a cut against it, so a
        # cut tries the reversed direction first; the other combinations remain as fall-backs
        first = (True, False) if op == "cut" else (False, True)
        tries = [(d, fl) for fl in (z0 < 0, z0 >= 0) for d in first] if t0 else [(first[0], False), (first[1], False)]
        nothing = 0
        for rev, flip in tries:
            self.model.ClearSelection2(True)
            sk.Select2(False, 0)
            if op == "add":
                f = fm.FeatureExtrusion3(True, False, rev, 0, 0, h, 0.0, False, False, False, False,
                                         0.0, 0.0, False, False, False, False, True, True, True,
                                         t0, off, flip)
            else:
                f = fm.FeatureCut4(True, False, rev, 0, 0, h, 0.0, False, False, False, False,
                                   0.0, 0.0, False, False, False, False, False, True, True,
                                   False, False, False, t0, off, flip, False)
            if f is None:
                print(f"  try reverse={rev} flip_start={flip}: SolidWorks returned no feature")
                nothing += 1
                continue
            v, c = self.mass()
            dc = math.dist(c, centroid)
            print(f"  try reverse={rev} flip_start={flip}: volume {v:.3f} mm^3, centre of mass "
                  f"({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}) mm")
            if abs(v - expected) <= 1e-4 * expected + 1e-3 and dc <= 1e-5 * SIZE + 1e-4:
                f.Name = name
                return v
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(0)       # keeps the sketch for the next try
        if nothing == len(tries):
            raise NoFeature(name)
        self.dump_sketch(sk)
        raise RuntimeError(f"{name}: no extrusion direction gives the expected volume {expected:.3f} mm^3 "
                           f"and centre of mass {centroid}")


class NoFeature(Exception):
    """SolidWorks created no feature from the sketch (typically an open or self-touching outline)."""


def _tangent(p, q):
    """Line and arc meet tangentially at their common corner (within 0.5 deg)."""
    ln, ar = (p, q) if p[0] == "line" else (q, p)
    c = circle_centre(*ar[1:])
    corner = q[1]
    d = (ln[2][0] - ln[1][0], ln[2][1] - ln[1][1])
    rad = (corner[0] - c[0], corner[1] - c[1])
    cosang = abs(d[0] * rad[0] + d[1] * rad[1]) / (math.hypot(*d) * math.hypot(*rad) + 1e-30)
    return cosang < math.sin(math.radians(0.5))


def circle_centre(p, q, r):
    ax, ay = p
    bx, by = q
    cx, cy = r
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    return ux, uy


def build():
    sw = SW()
    print(f"connected to SolidWorks; new part created; building {len(FEATURES)} feature(s)")
    for n, st in enumerate(FEATURES, 1):
        # entities are first added directly (AddToDB, no snapping); if SolidWorks then makes no
        # feature, the sketch is drawn again in the normal mode, which joins coincident endpoints
        # attempts: AddToDB with the planned relations and dimensions, AddToDB with none (the
        # solver cannot move anything), normal mode with none
        modes = [(True, True), (True, False), (False, False)]
        for m, (add_to_db, use_plans) in enumerate(modes):
            sk = sw.sketch(st["loops"], add_to_db, plans=st.get("plans") if use_plans else None)
            try:
                v = sw.feature(sk, st["op"], st["lo"], st["hi"], st["volume"], st["centroid"], st["name"])
                break
            except NoFeature:
                sw.dump_sketch(sk, everything=True)
                if m == len(modes) - 1:
                    raise RuntimeError(f"{st['name']}: SolidWorks makes no feature from the sketch")
                nxt = modes[m + 1]
                print(f"  no feature; drawing the sketch again ({'AddToDB' if nxt[0] else 'normal mode'}, "
                      f"{'with' if nxt[1] else 'without'} relations and dimensions)")
                sw.model.ClearSelection2(True)
                sk.Select2(False, 0)
                sw.model.Extension.DeleteSelection2(0)
        print(f"{st['name']}: volume {v:.3f} mm^3 (expected {st['volume']:.3f})")
    print(sw.rotate_back())
    for cy in CYLINDERS:
        v = sw.face_feature(cy)
        print(f"{cy['name']}: volume {v:.3f} mm^3 (expected {cy['volume']:.3f})")
    v0 = sw.mass()[0]
    print(sw.link_dimensions())
    v1 = sw.mass()[0]
    if abs(v1 - v0) > 1e-4 * v0 + 1e-3:
        raise RuntimeError(f"linking dimensions changed the volume from {v0:.3f} to {v1:.3f} mm^3")
    print(f"dimensions: {sw.dims_ok} added" + (f", {len(sw.dims_failed)} failed: {sw.dims_failed}" if sw.dims_failed else ""))
    print(f"relations: {sw.rels_ok} added" + (f", {len(sw.rels_failed)} failed: {sw.rels_failed}" if sw.rels_failed else ""))
    print("sketch status: " + ", ".join(sw.status))
    print("planned: " + ", ".join(st.get("expect", "?") for st in FEATURES))
    sw.restore_preferences()
    return sw


if __name__ == "__main__":
    sw = build()
    if len(sys.argv) > 1:
        sw.model.SaveAs3(sys.argv[1], 0, 0)
'''


def _round_plans(plans, nd=6):
    def r(v):
        if isinstance(v, (float, np.floating)):
            return round(float(v), nd)
        if isinstance(v, (list, tuple)):
            return type(v)(r(u) for u in v)
        return v
    return r(plans)


def solidworks_script(p25, name, nd=6):
    from .sketchplan import plan_outline
    k = p25["axis"]
    i, j = [m for m in range(3) if m != k]
    steps = _steps(p25)
    vols = expected_volumes(p25, steps)
    cyls = cylinder_steps(p25, steps)
    if cyls is None:
        raise ValueError("an off-axis feature does not start on a face perpendicular to its axis")
    R, back = BUILD_FRAME[k]
    allp = np.array([np.asarray(q, float) for ft in p25["features"]
                     for e in (ft["loop"].pattern[2] if ft["loop"].pattern is not None else ft["loop"].edges)
                     for q in e[1:] if np.ndim(q) == 1])
    helper = allp.min(axis=0) - 10.0
    size = float(np.linalg.norm(np.ptp(np.array([v[1] for v in vols]), axis=0))) + 1.0
    feats_pts = [np.asarray(q) for ft in p25["features"] for e in ft["loop"].edges for q in e[1:] if np.ndim(q) == 1]
    if feats_pts:
        P = np.array(feats_pts)
        size = float(np.linalg.norm(np.r_[np.ptp(P, axis=0), max(ft["z1"] for ft in p25["features"]) - min(ft["z0"] for ft in p25["features"])]))
    n_add = sum(s["op"] == "add" for s in steps)
    L = [f'"""{name}: SolidWorks feature tree generated by stl2cad (exact mode, 2.5D).',
         "",
         f"Builds the part in SolidWorks through its API: {n_add} Extruded Boss and "
         f"{len(steps) - n_add} Extruded Cut feature(s) along +{AX[k]},",
         ("each from a sketch on the Front Plane. Units: mm; the API is fed metres." if k == 2 else
          f"each from a sketch on the Front Plane, built with {AX[k]} along z and turned back by a "
          "Move/Copy Body feature. Units: mm."),
         "After every feature the part volume and centre of mass are checked against stl2cad's values.",
         "",
         "Needs: Windows, SolidWorks, Python with pywin32 (pip install pywin32).",
         f"Run:   python {name}.py [C:\\\\path\\\\part.SLDPRT]",
         '"""',
         "",
         f"AXIS_I, AXIS_J = {i}, {j}           # model axes of the sketch coordinates (a, b)",
         'PLANE_NAME = "Front Plane"',
         "PLANE_INDEX = 0                 # its position among the template planes",
         f"ROT = {[list(map(float, r)) for r in R]}   # STL frame -> build frame (extrusion axis {AX[k]} -> z)",
         f"ROT_BACK = {tuple(back)}   # rotation (rad) about x, y, z back to the STL frame",
         f"FINAL_CENTROID = ({vols[-1][1][0]:.5f}, {vols[-1][1][1]:.5f}, {vols[-1][1][2]:.5f})   # STL frame, mm",
         f"SIZE = {_f(size, 3)}                  # part diagonal (mm), scale of the centre-of-mass check",
         f"HELPER_AB = ({_f(helper[0], 3)}, {_f(helper[1], 3)})   # helper circle (a, b), outside the part",
         "",
         "# ---- features, in order: sketch loops (a, b in mm), levels, expected volume (mm^3) ----",
         "FEATURES = ["]
    cnt = {"add": 0, "cut": 0, "hole": 0}
    for st, v in zip(steps, vols):
        if st["holes"]:
            cnt["hole"] += 1
            ds = sorted({round(2 * l.edges[0][2], 4) for l in st["loops"]})
            label = f"Holes {cnt['hole']} (dia {', '.join(f'{d:g}' for d in ds)})"
        elif st["op"] == "add":
            cnt["add"] += 1
            label = f"Boss {cnt['add']}"
        else:
            cnt["cut"] += 1
            label = f"Cut {cnt['cut']}"
        label += f" {AX[k]} {st['z0']:.4g} to {st['z1']:.4g}"
        cb = np.asarray(R, float) @ np.asarray(v[1])          # expected centroid in the build frame
        L.append("    dict(")
        L.append(f'        name="{label}", op="{st["op"]}", lo={_f(st["z0"], nd)}, hi={_f(st["z1"], nd)}, '
                 f"volume={v[0]:.4f}, centroid=({cb[0]:.5f}, {cb[1]:.5f}, {cb[2]:.5f}), holes={st['holes']},")
        plans, full = [], True
        for loop in st["loops"]:
            if not (len(loop.edges) == 1 and loop.edges[0][0] in ("circle", "ellipse")):
                lin = None if getattr(loop, "linear", None) is None else loop.linear[:2]
                rot = None if loop.pattern is None else (len(loop.pattern[2]), loop.pattern[1],
                                                         tuple(map(float, loop.pattern[0])))
                ops, ok = plan_outline([tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in loop.edges],
                                       linear=lin, rotational=rot)
                plans.append(ops)
                full = full and ok
        L.append(f"        expect=\"{'fully defined' if full else 'under defined'}\",")
        L.append("        plans=" + repr(_round_plans(plans)) + ",")
        L.append("        loops=[")
        for loop in st["loops"]:
            L.append("        " + _loop_code(loop, nd) + ",")
        L.append("        ],")
        L.append("    ),")
    L.append("]")
    L.append("")
    L.append("# ---- off-axis bosses and holes, after the prisms, each sketched on its start face (STL frame, mm) ----")
    L.append("CYLINDERS = [")
    for cy in cyls:
        L.append(f'    dict(name="{cy["name"]}", op="{cy["op"]}", pick={tuple(round(v, nd) for v in cy["pick"])}, '
                 f'normal={tuple(round(v, 9) for v in cy["normal"])}, '
                 f'centre={tuple(round(v, nd) for v in cy["centre"])},')
        L.append(f'         dia={_f(cy["dia"], nd)}, depth={_f(cy["depth"], nd)}, depth2={_f(cy["depth2"], nd)}, volume={cy["volume"]:.4f}, '
                 f'centroid=({cy["centroid"][0]:.5f}, {cy["centroid"][1]:.5f}, {cy["centroid"][2]:.5f})),')
    L.append("]")
    return "\n".join(L) + "\n" + _RUNTIME


def verify_with_standin(script, ref_shape, tol, cut_reversed=True, boss_reversed=False, timeout=900):
    """Run the script against the CadQuery stand-in for the SolidWorks API (swmock) and compare
    the saved body with ref_shape. This checks the script's geometry and its direction handling,
    not SolidWorks itself. Returns (ok, message)."""
    import os
    import subprocess
    import sys
    import tempfile
    import cadquery as cq
    from .codegen import _max_distance
    d = tempfile.mkdtemp(prefix="stl2cad_sw_")
    py, out = os.path.join(d, "part.py"), os.path.join(d, "part.step")
    with open(py, "w") as fh:
        fh.write(script)
    env = dict(os.environ, PYTHONPATH=os.path.join(os.path.dirname(__file__), "swmock"),
               SWMOCK_CUT_REVERSED="1" if cut_reversed else "0",
               SWMOCK_BOSS_REVERSED="1" if boss_reversed else "0")
    r = subprocess.run([sys.executable, py, out], capture_output=True, text=True, timeout=timeout, env=env)
    if r.returncode != 0 or not os.path.exists(out):
        tail = (r.stderr or r.stdout).strip().splitlines()[-1:] or ["no output"]
        return False, "stand-in run failed: " + tail[0][:200]
    import json
    import re
    s = cq.importers.importStep(out).val()
    dv = abs(s.Volume() - ref_shape.Volume()) / ref_shape.Volume()
    m = re.search(r"dimensions: (\d+) added(, (\d+) failed)?", r.stdout)
    m2 = re.search(r"relations: (\d+) added(, (\d+) failed)?", r.stdout)
    if m is None or m.group(3) or m2 is None or m2.group(3):
        return False, ("dimensions or relations failed on the API stand-in: "
                       + (m.group(0) if m else "no report") + "; " + (m2.group(0) if m2 else "no report"))
    dims = json.load(open(out + ".dims.json"))["dims"]
    dist = max(_max_distance(ref_shape, s, tol), _max_distance(s, ref_shape, tol))
    ok = dv < 1e-4 and dist < 0.5 * tol
    return ok, f"API stand-in: volume {100 * dv:.4f} %, max distance {dist:.4f} mm, {len(dims)} dimension(s), {m2.group(1)} relation(s)"
