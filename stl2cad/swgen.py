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

from .core import mass_props

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
        vols.append(mass_props(solid))           # volume and centre of mass
    return vols


def cylinder_steps(p25, steps):
    """Off-axis bosses and holes, applied after the prisms (in the STL frame): for each, the
    point to pick its start face, the circle centre on that face, diameter, length, and the
    volume and centre of mass expected after it. None if a feature cannot be sketched on a face
    perpendicular to its axis (no face point was found)."""
    import cadquery as cq
    from .prism25 import cylinder_solid
    solid = _prism_solid(p25, steps)
    for ft in p25.get("revolved") or []:
        solid = _apply_revolved(solid, ft)
    out = []
    for n, cy in enumerate(p25.get("cylinders") or [], 1):
        cen, d = np.asarray(cy["c"]), np.asarray(cy["d"])
        offset = 0.0

        def face_at(pick):
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
            return None if best is None or best[0] > 0.5 else best
        # pick point and circle centre on the model's planar face (the mesh lies within the fit
        # tolerance of it, the API picks on the model)
        best = None if cy.get("face_pt") is None else face_at(np.asarray(cy["face_pt"], float))
        if best is not None:
            pick = np.asarray(cy["face_pt"], float)
        else:
            # no flat face where the feature starts, in the model at this step (a cross bore
            # entering a curved wall, a step face not cut yet): sketch on another planar face
            # perpendicular to the axis and start the feature at an offset from it
            alt = _sketch_face(solid, cen, d, cy["t0"])
            if alt is None:
                return None
            pick, s_ = alt
            offset = cy["t0"] - s_
            cy = dict(cy, tf=s_)
            best = face_at(pick)
            if best is None:
                return None
        ctr = cen + cy["tf"] * d
        o = np.array(best[1].Center().toTuple())
        nrm = best[2]
        pick = pick - ((pick - o) @ nrm) * nrm
        ctr = ctr - ((ctr - o) @ nrm) * nrm
        tip = cy.get("tip")
        body = cylinder_solid(dict(cy, tip=None))
        solid = solid.fuse(body).clean() if cy["op"] == "add" else solid.cut(body).clean()
        vol_, c = mass_props(solid)
        depth = cy["t1"] - cy["tf"] - offset
        out.append(dict(op=cy["op"], pick=tuple(map(float, pick)), normal=tuple(map(float, nrm)),
                        centre=tuple(map(float, ctr)), dia=2 * cy["r"], offset=offset,
                        depth=depth, depth2=(cy["tf"] - cy["t0"]) if cy.get("spans") else 0.0,
                        volume=vol_,
                        centroid=tuple(c), through=cy.get("through", False),
                        name=f"{'Boss' if cy['op'] == 'add' else 'Hole'} {n} dia {2 * cy['r']:g} "
                             f"along ({d[0] + 0:.3g}, {d[1] + 0:.3g}, {d[2] + 0:.3g})"))
        if tip:
            # drill point: a cut from the hole's end, tapered to a
            # point (the wall at the half angle to the axis; a hair short of the point: an
            # extrusion cannot taper to nothing)
            import cadquery as cq
            p1 = cen + cy["t1"] * d
            solid = solid.cut(cq.Solid.makeCone(cy["r"], 0.0, tip["h"], cq.Vector(*p1), cq.Vector(*d))).clean()
            vol_, c = mass_props(solid)
            out.append(dict(op="cut", pick=tuple(map(float, pick)), normal=tuple(map(float, nrm)),
                            centre=tuple(map(float, ctr)), dia=2 * cy["r"], offset=offset + depth,
                            depth=tip["h"] * (1 - 1e-6), depth2=0.0, volume=vol_, centroid=tuple(c),
                            through=False, draft=tip["angle"] / 2, dims=False,
                            name=f"Drill point {tip['angle']:g} deg of hole {n}"))
    return out


def _sketch_face(solid, c, d, t0):
    """A planar face of solid perpendicular to the axis (c, d), to sketch a feature on that starts
    at c + t0 d: the one whose plane is nearest to t0 (behind it preferred). Returns (a point on
    the face, the plane's position s along the axis) or None."""
    best = None
    for fc in solid.Faces():
        if fc.geomType() != "PLANE":
            continue
        o = np.array(fc.Center().toTuple())
        nrm = np.array(fc.normalAt(fc.Center()).toTuple())
        if abs(nrm @ d) < np.cos(np.radians(0.5)):
            continue
        s_ = float((o - c) @ d)
        key = (s_ > t0 + 1e-9, abs(t0 - s_))
        if best is None or key < best[0]:
            verts, tris = fc.tessellate(1e-3)
            V = np.array([v.toTuple() for v in verts])
            T = max(tris, key=lambda t: np.linalg.norm(np.cross(V[t[1]] - V[t[0]], V[t[2]] - V[t[0]])))
            best = (key, V[list(T)].mean(axis=0), s_)
    return None if best is None else (best[1], best[2])


def _prism_solid(p25, steps):
    """The extrusions and cuts of the steps, in the STL frame (CadQuery)."""
    import cadquery as cq
    from .prism25 import _wire
    k = p25["axis"]
    solid = None
    for st in steps:
        for loop in st["loops"]:
            f = cq.Face.makeFromWires(_wire(loop, k, st["z0"]), [])
            D = [0.0, 0.0, 0.0]
            D[k] = st["z1"] - st["z0"]
            body = cq.Solid.extrudeLinear(f, cq.Vector(*D))
            solid = body if solid is None else (solid.fuse(body) if st["op"] == "add" else solid.cut(body))
    return solid


def _apply_revolved(solid, ft):
    from .revfeat import region_solid
    body = region_solid(ft["loops"], np.asarray(ft["a3"], float), np.asarray(ft["e1"], float), ft["k"])
    if ft["op"] == "cut":
        return solid.cut(body)
    return solid.fuse(body)


def revolve_steps(p25, steps, nd=6):
    """Revolved features about axes along the prism axis, for the SolidWorks script. API sketches
    land on the Front Plane, so the body is turned with the prism axis along y (REVOLVE_FRAME) and
    moved, for each axis, so that the axis is the y axis; the profile is then sketched on the Front
    Plane with u along +x (a 360-degree revolve does not depend on the half plane) and a vertical
    centerline at x = 0, and the body is moved back. Consecutive features on one axis share the
    move; consecutive ones of one kind (boss, cut) share a sketch.
    Returns (turn, groups): turn = expected centre of mass after the turn into the revolve frame
    (None when the prism axis is y: the STL frame is that frame), groups as dicts."""
    from .prism25 import Loop
    from .sketchplan import plan_outline
    revs = p25.get("revolved") or []
    if not revs:
        return None, []
    k = p25["axis"]
    Rk = np.asarray(REVOLVE_FRAME[k][0], float)
    solid = _prism_solid(p25, steps)
    com = lambda s_: np.array(mass_props(s_)[1])
    turn = None if k == 1 else tuple(map(float, Rk @ com(solid)))
    groups = []
    for ft in revs:
        a3 = np.asarray(ft["a3"], float)
        q = Rk @ a3
        shift = np.array([-q[0], 0.0, -q[2]])
        if not groups or np.linalg.norm(groups[-1]["shift_v"] - shift) > 1e-9:
            groups.append(dict(shift_v=shift, steps=[], before=Rk @ com(solid) + shift))
        g = groups[-1]
        solid = _apply_revolved(solid, ft)
        stp = g["steps"][-1] if g["steps"] and g["steps"][-1]["op"] == ft["op"] else None
        if stp is None:
            stp = dict(op=ft["op"], loops=[])
            g["steps"].append(stp)
        stp["loops"] += [Loop([(e[0],) + tuple(np.asarray(p, float) for p in e[1:]) for e in L.edges], L.area)
                         for L in ft["loops"]]
        stp["volume"] = mass_props(solid)[0]
        stp["centroid"] = tuple(map(float, Rk @ com(solid) + shift))
        g["after"] = Rk @ com(solid)                  # centre of mass once moved back
        g["axis"] = tuple(map(float, a3))
    n = 0
    for g in groups:
        for stp in g["steps"]:
            n += 1
            plans, full = [], True
            for m_, L in enumerate(stp["loops"]):
                E_ = [tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in L.edges]
                ops, ok = plan_outline(E_, revolve_axis=(m_ == 0))
                plans.append(ops)
                full = full and ok
            P = np.array([q for L in stp["loops"] for e in L.edges for q in e[1:]])
            stp["helper"] = tuple(map(float, P.min(axis=0) - [0.0, 10.0]))
            stp["plans"], stp["full"] = plans, full
            stp["name"] = (f"Revolved {'Boss' if stp['op'] == 'add' else 'Cut'} {n} about "
                           f"({', '.join(f'{v:g}' for v in g['axis'])}) along {AX[k]}")
    return turn, groups


def seed_copy(E, p, N):
    """Index of the pattern copy used as the seed: the one at the left (else bottom) end."""
    a, b = tuple(E[0][1]), tuple(E[(N - 1) * p][1])
    return N - 1 if (b[0], b[1]) < (a[0], a[1]) else 0


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


def seed_copy(E, p, N):
    """Index of the pattern copy used as the seed: the one at the left (else bottom) end."""
    a, b = tuple(E[0][1]), tuple(E[(N - 1) * p][1])
    return N - 1 if (b[0], b[1]) < (a[0], a[1]) else 0


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
        self.direct, self.helper_ab = False, None    # revolve profiles: sketch (x, y) = (a, b)
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
        if self.direct:                          # revolve profiles: Front Plane (x, y) as given
            return a * MM, b * MM
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
        p0 = self.to_sketch(*(self.helper_ab or HELPER_AB))   # away from the part, so selections do not hit it
        helper = sm.CreateCircleByRadius(p0[0], p0[1], 0.0, HELPER_R * MM)    # may return None
        self.ents = []                           # (kind, sketch geometry (m), entity) as drawn
        patterns = []                            # (entities of the period, count, step) to repeat
        circles = []                             # (recorded period entities, loop) to turn
        for loop in loops:
            linear = isinstance(loop, dict) and "step" in loop
            circular = isinstance(loop, dict) and "centre" in loop and self.use_patterns
            first = len(self.ents)
            todo = expand(loop)
            if linear:                           # the seed copy and the rest; the others by pattern
                p_, n_ = len(loop["period"]), loop["count"]
                s_ = seed_copy(todo, p_, n_)
                todo = todo[s_ * p_:(s_ + 1) * p_] + todo[p_ * n_:]
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
                step = loop["step"] if s_ == 0 else tuple(-v for v in loop["step"])
                patterns.append(([g for _, _, g in self.ents[first:first + p_]], n_, step))
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
                    ok = (g is not None and g.Select4(False, NOTHING)) or \
                        ext.SelectByID2("", "SKETCHSEGMENT", cx + r, cy, 0.0, False, 0, NOTHING, 0)
                    d = self.model.AddDimension2(cx + 1.6 * r, cy + 1.6 * r, 0.0) if ok else None
                    self._count(d, tag + " diameter")
                def pick_centre(append=False):
                    try:                         # the circle's own centre point, else by coordinates
                        c_ = call0(g, "GetCenterPoint2") if g is not None else None
                        if c_ is not None and c_.Select4(append, NOTHING):
                            return True
                    except Exception:
                        pass
                    return ext.SelectByID2("", "SKETCHPOINT", cx, cy, 0.0, append, 0, NOTHING, 0)
                self.centre_to_origin(tag, cx, cy, pick_centre)
                for kind, tx, ty in (("H", cx / 2, cy - 3 * r), ("V", cx - 3 * r, cy / 2)):
                    if abs(cx if kind == "H" else cy) < 1e-9:
                        continue                 # centre on the origin's line: a relation above
                    self.model.ClearSelection2(True)
                    ok1 = pick_centre()
                    ok2 = ext.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0, True, 0, NOTHING, 0)
                    d = None
                    if ok1 and ok2:
                        d = (self.model.AddHorizontalDimension2(tx, ty, 0.0) if kind == "H"
                             else self.model.AddVerticalDimension2(tx, ty, 0.0))
                    self._count(d, tag + (" x" if kind == "H" else " y") + " position")
            except Exception as ex:
                self.dims_failed.append(f"{tag}: {type(ex).__name__}: {ex}")
        self.model.ClearSelection2(True)

    def centre_to_origin(self, tag, cx, cy, pick_centre):
        """A circle centre on the origin, or on its vertical or horizontal line, needs a relation in
        place of the distance that is zero: Coincident, Vertical or Horizontal with the origin."""
        zx, zy = abs(cx) < 1e-9, abs(cy) < 1e-9
        if not (zx or zy):
            return
        rel, what = (("sgCOINCIDENT", "on the origin") if zx and zy else
                     ("sgVERTICALPOINTS2D", "vertical with the origin") if zx else
                     ("sgHORIZONTALPOINTS2D", "horizontal with the origin"))
        self.model.ClearSelection2(True)
        ok = pick_centre(False) and self.model.Extension.SelectByID2(
            "Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0, True, 0, NOTHING, 0)
        if ok:
            self.model.SketchAddConstraints(rel)
            self.rels_ok += 1
        else:
            self.rels_failed.append(f"{tag} centre {what} (selection)")
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
                    if not _on_arc(P, geo[0], geo[4], geo[1], geo[2]):
                        d += 1.0                 # another arc on the same circle (snapring)
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

    def _swapped(self):
        """Whether the plan's first axis is the sketch's vertical (parts extruded along x: a plan
        point (a, b) is the sketch point (-b, a)); then the plan's horizontal is vertical."""
        x0, y0 = self.to_sketch(0.0, 0.0)
        x1, y1 = self.to_sketch(1.0, 0.0)
        return abs(y1 - y0) > abs(x1 - x0)

    def _relation(self, rel, picks, what):
        if self._swapped():
            rel = {"sgHORIZONTAL2D": "sgVERTICAL2D", "sgVERTICAL2D": "sgHORIZONTAL2D",
                   "sgHORIZONTALPOINTS2D": "sgVERTICALPOINTS2D",
                   "sgVERTICALPOINTS2D": "sgHORIZONTALPOINTS2D"}.get(rel, rel)
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
                if self._swapped():                  # the plan's horizontal is the sketch's vertical
                    kind = {"H": "V", "V": "H"}.get(kind, kind)
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
        pattern (SketchManager.CreateLinearSketchStepAndRepeat) with its spacing and count
        dimensions. In SolidWorks 2025 the copies share their end points with each other and with
        the rest of the outline (sw_probe8), and follow changes of the seed only when the pattern
        has its dimensions (sw_probe13)."""
        p0, p1 = self.to_sketch(0.0, 0.0), self.to_sketch(*step)
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        self.model.ClearSelection2(True)
        ok = all(g.Select4(True, NOTHING) for g in seed)
        r = self.model.SketchManager.CreateLinearSketchStepAndRepeat(
            count, 1, math.hypot(dx, dy), 0.0, math.atan2(dy, dx), math.pi / 2, "",
            True, False, False, True, False) if ok else False
        self.model.ClearSelection2(True)
        if not r:
            raise RuntimeError(f"sketch linear pattern of {len(seed)} entities x {count} was not created")
        print(f"  sketch pattern: {len(seed)} entities x {count}, step {math.hypot(dx, dy) / MM:g} mm")

    def apply_plan(self, plan):
        """Relations and dimensions planned by stl2cad for one outline, in order. A step that acts
        on the wrong entity shows as a feature of the wrong shape; the sketch is then drawn again
        without its plan (build_step). (Checking each step and undoing it with EditUndo2 did more
        harm than good in SolidWorks 2025: an undo could take back more than the step.)"""
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
                if not self.direct:              # revolve profiles: re-solving one by an equation
                    # turned topcover's dish into a huge body in SolidWorks 2025; equal radii there
                    # are coincidences, not design intent
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
            self.n_equations = call0(em, "GetCount")     # equations before linking
            n = 0
            for (kind, val), names in sorted(groups.items()):
                var = f"{'radius' if kind == 'R' else 'dia'}_{val:g}".replace(".", "_").replace("-", "m")
                em.Add2(-1, f'"{var}" = {val:g}mm', True)   # units explicit: the template may be in metres
                for nm in names:
                    em.Add2(-1, f'"{nm}" = "{var}"', True)
                    n += 1
            self.linked = [(nm, val) for (kind, val), names in groups.items() for nm in names]
            call0(self.model, "EditRebuild3")
            return f"linked dimensions: {n} in {len(groups)} global variable(s)"
        except Exception as ex:
            return f"linked dimensions: failed ({type(ex).__name__}: {ex})"

    def unlink_dimensions(self):
        """Delete the equations link_dimensions added (from the last), restore the linked
        dimensions' values, and rebuild."""
        em = call0(self.model, "GetEquationMgr")
        for idx in range(call0(em, "GetCount") - 1, getattr(self, "n_equations", 0) - 1, -1):
            em.Delete(idx)
        for nm, val in getattr(self, "linked", []):  # an equation leaves its value behind
            try:
                self.model.Parameter(nm).SystemValue = val * MM
            except Exception as ex:
                print(f"  could not restore {nm}: {type(ex).__name__}: {ex}")
        call0(self.model, "EditRebuild3")

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
        r = self.quarter_turn(FINAL_CENTROID, "Rotate to STL orientation")
        if r is not None:
            return "rotated back to the STL's orientation " + r
        self.turned_back = False
        return ("NOTE: part left in the build frame (extrusion axis along z); no quarter turn brought it "
                "to the STL's orientation")

    def turn_y(self, angle, target, name):
        """Move/Copy Body turning the body by the signed angle (rad, right-handed) about the y
        axis: angle argument 2, which turned about y in the stand-in's convention (argument 1
        about z and argument 3 about x were seen in SolidWorks 2025). The sign is not searched:
        a body symmetric about y has the same centre of mass either way, so a search could turn
        it the wrong way. The centre of mass is checked; other arguments are the fall-back."""
        f = None
        try:
            bodies = self.model.GetBodies2(0, True)
            self.model.ClearSelection2(True)
            sd = call0(self.model.SelectionManager, "CreateSelectData")
            sd.Mark = 1
            bodies[0].Select2(False, sd)
            f = self.model.FeatureManager.InsertMoveCopyBody2(
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, angle, 0.0, False, 1)
        except Exception as e:
            print(f"  turn about y failed: {type(e).__name__}: {e}")
        if f is not None:
            v, c = self.mass()
            if math.dist(c, target) <= 1e-4 * SIZE + 1e-4:
                f.Name = name
                return f"(Move/Copy Body, angle argument 2, {math.degrees(angle):+g} deg)"
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(0)
        return self.quarter_turn(target, name, angle=abs(angle))

    def turn_about(self, m, angle, target, name):
        """Move/Copy Body turning the body by the signed angle (rad, right-handed) about model
        axis m: angle argument 1 turned about z, 2 about y and 3 about x in SolidWorks 2025. The
        centre of mass is checked; a quarter-turn search is the fall-back."""
        slot = {2: 0, 1: 1, 0: 2}[m]
        ang = [0.0, 0.0, 0.0]
        ang[slot] = angle
        f = None
        try:
            bodies = self.model.GetBodies2(0, True)
            self.model.ClearSelection2(True)
            sd = call0(self.model.SelectionManager, "CreateSelectData")
            sd.Mark = 1
            bodies[0].Select2(False, sd)
            f = self.model.FeatureManager.InsertMoveCopyBody2(
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ang[0], ang[1], ang[2], False, 1)
        except Exception as e:
            print(f"  turn about {'xyz'[m]} failed: {type(e).__name__}: {e}")
        if f is not None:
            v, c = self.mass()
            if math.dist(c, target) <= 1e-4 * SIZE + 1e-4:
                f.Name = name
                return f"(Move/Copy Body, angle argument {slot + 1}, {math.degrees(angle):+g} deg)"
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(0)
        return self.quarter_turn(target, name, angle=abs(angle))

    def quarter_turn(self, target, name, angle=math.pi / 2):
        """Move/Copy Body turning the body by angle (a quarter turn by default) about an axis
        through the origin such that its centre of mass comes to target (mm); each sign in each
        angle argument is tried. Returns a description, or None if no turn does it."""
        for slot in range(3):
            for sign in (1, -1):
                ang = [0.0, 0.0, 0.0]
                ang[slot] = sign * angle
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
                if math.dist(c, target) <= 1e-4 * SIZE + 1e-4:
                    f.Name = name
                    return f"(Move/Copy Body, angle argument {slot + 1}, {sign * math.degrees(angle):+g} deg)"
                self.model.ClearSelection2(True)
                f.Select2(False, 0)
                self.model.Extension.DeleteSelection2(0)
        return None

    def translate(self, shift, target, name):
        """Move/Copy Body moving the body by shift (mm). The meaning of the translation arguments
        of InsertMoveCopyBody2 is not assumed: the vector itself (metres), then a unit direction
        with the distance, are tried; the one that brings the centre of mass to target is kept."""
        d = math.sqrt(sum(v * v for v in shift))
        tries = [("vector", [v * MM for v in shift], 0.0), ("direction and distance", [v / d for v in shift], d * MM)]
        for how, t, dist in tries:
            f = None
            try:
                bodies = self.model.GetBodies2(0, True)
                self.model.ClearSelection2(True)
                sd = call0(self.model.SelectionManager, "CreateSelectData")
                sd.Mark = 1
                bodies[0].Select2(False, sd)
                f = self.model.FeatureManager.InsertMoveCopyBody2(
                    t[0], t[1], t[2], dist, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 1)
            except Exception as e:
                print(f"  move attempt ({how}) failed: {type(e).__name__}: {e}")
            if f is None:
                continue
            v, c = self.mass()
            if math.dist(c, target) <= 1e-4 * SIZE + 1e-4:
                f.Name = name
                return f"{name}: moved by ({', '.join(f'{x:g}' for x in shift)}) mm (Move/Copy Body, {how})"
            print(f"  move ({how}): centre of mass ({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}), expected {tuple(round(x, 3) for x in target)}")
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(0)
        raise RuntimeError(f"{name}: no Move/Copy Body form moves the body by {tuple(shift)} mm")

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
        circ = sm.CreateCircleByRadius(cx, cyy, 0.0, r)
        tag = cy["name"]

        def to_model(sx, sy):
            """Sketch point (m) -> model point (m): the inverse of the transform above (selection
            and dimension placement take model coordinates; on a face sketch they differ)."""
            v = (sx - A[9], sy - A[10], -A[11])
            return tuple(sum(v[i] * A[3 * r_ + i] for i in range(3)) for r_ in range(3))

        def pick_circle(append):
            if circ is not None and circ.Select4(append, NOTHING):
                return True
            return ext.SelectByID2("", "SKETCHSEGMENT", *to_model(cx + r, cyy), append, 0, NOTHING, 0)

        def pick_centre(append):
            try:
                c = call0(circ, "GetCenterPoint2")
                if c is not None and c.Select4(append, NOTHING):
                    return True
            except Exception:
                pass
            return ext.SelectByID2("", "SKETCHPOINT", *to_model(cx, cyy), append, 0, NOTHING, 0)
        try:
            if not cy.get("dims", True):             # a drill point: its circle is the hole's
                raise StopIteration
            self.model.ClearSelection2(True)
            ok = pick_circle(False)
            self._count(self.model.AddDimension2(*to_model(cx + 1.6 * r, cyy + 1.6 * r)) if ok else None,
                        tag + " diameter")
            self.centre_to_origin(tag, cx, cyy, pick_centre)
            for kind, tx, ty in (("H", cx / 2, cyy - 3 * r), ("V", cx - 3 * r, cyy / 2)):
                if abs(cx if kind == "H" else cyy) < 1e-9:
                    continue
                self.model.ClearSelection2(True)
                ok1 = pick_centre(False)
                ok2 = ext.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0, True, 0, NOTHING, 0)
                d = None
                if ok1 and ok2:
                    d = (self.model.AddHorizontalDimension2(*to_model(tx, ty)) if kind == "H"
                         else self.model.AddVerticalDimension2(*to_model(tx, ty)))
                self._count(d, tag + (" x" if kind == "H" else " y") + " position on the face")
        except StopIteration:
            pass
        except Exception as ex:
            self.dims_failed.append(f"{tag}: {type(ex).__name__}: {ex}")
        sm.InsertSketch(True)
        if cy.get("dims", True):
            self.report_status()
        sk = self.model.FeatureByPositionReverse(0)
        expected, centroid = cy["volume"], self.in_model_frame(cy["centroid"])
        h = cy["depth"] * MM
        off = abs(cy.get("offset", 0.0)) * MM       # start offset from the sketch face (swStartOffset)
        t0 = 3 if off > 0 else 0
        tries = [(rev, flip) for rev in (False, True) for flip in ((False, True) if off > 0 else (False,))]
        dang = math.radians(cy.get("draft", 0.0))
        if dang:                                     # tapered cut (drill point): inward first
            tries = [(rev, flip, out) for out in (False, True) for rev, flip in tries]
        else:
            tries = [(rev, flip, False) for rev, flip in tries]
        for rev, flip, outward in tries:
            self.model.ClearSelection2(True)
            sk.Select2(False, 0)
            if dang:
                f = fm.FeatureCut4(True, False, rev, 0, 0, h, 0.0, True, False, outward, False,
                                   dang, 0.0, False, False, False, False, False, True, True,
                                   False, False, False, t0, off, flip, False)
            elif cy["op"] == "add":
                h2 = cy.get("depth2", 0.0) * MM          # a boss through the body: both directions
                f = fm.FeatureExtrusion3(h2 <= 0, False, rev, 0, 0, h, h2, False, False, False, False,
                                         0.0, 0.0, False, False, False, False, True, True, True,
                                         t0, off, flip)
            else:
                f = fm.FeatureCut4(True, False, rev, 0, 0, h, 0.0, False, False, False, False,
                                   0.0, 0.0, False, False, False, False, False, True, True,
                                   False, False, False, t0, off, flip, False)
            if f is None:
                print(f"  try reverse={rev} flip_start={flip}: SolidWorks returned no feature")
                continue
            v, c = self.mass()
            print(f"  try reverse={rev} flip_start={flip}: volume {v:.3f} mm^3, centre of mass ({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}) mm")
            if abs(v - expected) <= 1e-4 * expected + 1e-3 and math.dist(c, centroid) <= 1e-5 * SIZE + 1e-4:
                f.Name = cy["name"]
                return v
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(0)
        raise RuntimeError(f"{cy['name']}: neither direction gives the expected volume {expected:.3f} mm^3")

    def blend(self, b):
        """Fillet or chamfer on the edges through b's points (model coordinates): the edges are
        picked with SelectByID2 (mark 1, as a recorded fillet macro does; then mark 0), and the
        feature is checked against the expected volume and centre of mass. A fillet is tried with
        FeatureFillet3 (array arguments empty, then Nothing) and with the older FeatureFillet; a
        chamfer with InsertFeatureChamfer, equal distances."""
        fm, ext = self.model.FeatureManager, self.model.Extension
        E = VARIANT(pythoncom.VT_EMPTY, None)
        r = b["size"] * MM
        if b["kind"] == "fillet":
            makes = [("FeatureFillet3", lambda: fm.FeatureFillet3(195, r, r, 0, 0, 0, 0, E, E, E, E, E, E, E)),
                     ("FeatureFillet3 (Nothing)", lambda: fm.FeatureFillet3(195, r, r, 0, 0, 0, 0, NOTHING, NOTHING,
                                                                            NOTHING, NOTHING, NOTHING, NOTHING, NOTHING)),
                     ("FeatureFillet", lambda: fm.FeatureFillet(195, r, 0, 0, E, E, E))]
        else:
            makes = [("InsertFeatureChamfer", lambda: fm.InsertFeatureChamfer(0, 2, r, math.radians(45.0), r, 0, 0, 0)),
                     ("InsertFeatureChamfer (tangent propagation)",
                      lambda: fm.InsertFeatureChamfer(1, 2, r, math.radians(45.0), r, 0, 0, 0))]
        for mark in (1, 0):
            for how, make in makes:
                self.model.ClearSelection2(True)
                ok = all(ext.SelectByID2("", "EDGE", *(v * MM for v in self.in_model_frame(p)), True, mark, NOTHING, 0)
                         for p in b["points"])
                if not ok:
                    print(f"  {b['name']}: an edge was not picked (mark {mark})")
                    continue
                if mark == 1 and how == makes[0][0]:
                    self.report_selection(b)
                f = None
                try:
                    f = make()
                except Exception as ex:
                    print(f"  {b['name']}: {how} failed: {type(ex).__name__}: {ex}")
                if f is None:
                    print(f"  {b['name']}: {how} (mark {mark}) made no feature")
                    continue
                ok_ = lambda v, c: (abs(v - b["volume"]) <= 1e-4 * b["volume"] + 1e-3
                                    and math.dist(c, self.in_model_frame(b["centroid"])) <= 1e-5 * SIZE + 1e-4)
                v, c = self.mass()
                if ok_(v, c):
                    f.Name = b["name"]
                    print(f"  {b['name']}: {how}, mark {mark}")
                    return v
                print(f"  {b['name']}: {how} (mark {mark}): volume {v:.3f} mm^3, expected {b['volume']:.3f}")
                # what SolidWorks made, read back from the feature definition; for a fillet, the
                # radius is set there and the feature rebuilt
                if self.fix_blend(f, b, r):
                    v, c = self.mass()
                    if ok_(v, c):
                        f.Name = b["name"]
                        print(f"  {b['name']}: {how}, radius set through the feature definition")
                        return v
                    print(f"  {b['name']}: after setting the radius: volume {v:.3f} mm^3")
                self.model.ClearSelection2(True)
                f.Select2(False, 0)
                self.model.Extension.DeleteSelection2(0)
        raise RuntimeError(f"{b['name']}: no fillet or chamfer call gives the expected volume")

    def report_selection(self, b):
        """Print what the edge picks selected: count, entity type, and for an edge its end points
        (mm), to see whether SelectByID2 took the edge at the point or another one."""
        sm = self.model.SelectionManager
        try:
            n = sm.GetSelectedObjectCount2(-1)
        except Exception as ex:
            print(f"    selection count failed: {type(ex).__name__}: {ex}")
            return
        print(f"    {n} selected (picked {len(b['points'])} point(s))")
        for i in range(1, n + 1):
            try:
                typ = sm.GetSelectedObjectType3(i, -1)
                o = sm.GetSelectedObject6(i, -1)
                ends = ""
                try:
                    sv, ev = call0(o, "GetStartVertex"), call0(o, "GetEndVertex")
                    p0 = [v / MM for v in call0(sv, "GetPoint")] if sv is not None else None
                    p1 = [v / MM for v in call0(ev, "GetPoint")] if ev is not None else None
                    ends = f", ends {None if p0 is None else tuple(round(v, 3) for v in p0)} -> " \
                           f"{None if p1 is None else tuple(round(v, 3) for v in p1)}"
                except Exception as ex:
                    ends = f", ends unavailable ({type(ex).__name__})"
                print(f"    #{i}: type {typ}{ends}")
            except Exception as ex:
                print(f"    #{i}: {type(ex).__name__}: {ex}")

    def fix_blend(self, f, b, r):
        """Print the fillet or chamfer definition SolidWorks made (radius or distance, edge
        count); for a fillet whose radius differs, set it and modify the feature. True if the
        feature was modified."""
        try:
            d = call0(f, "GetDefinition")
        except Exception as ex:
            print(f"    no definition: {type(ex).__name__}: {ex}")
            return False
        info = []
        for nm in ("DefaultRadius", "DefaultConicRhoOrRadius", "Distance", "EdgeChamferDistance", "Type", "FilletType"):
            try:
                info.append(f"{nm}={getattr(d, nm)}")
            except Exception:
                pass
        for nm in ("GetEdgeCount", "GetEntitiesToChamferCount", "GetFaceCount"):
            try:
                info.append(f"{nm}={getattr(d, nm)() if nm != 'GetFaceCount' else ''}")
            except Exception:
                pass
        print("    definition: " + ", ".join(info))
        if b["kind"] != "fillet":
            return False
        try:
            if abs(float(d.DefaultRadius) - r) < 1e-9:
                return False
            d.DefaultRadius = r
            return bool(f.ModifyDefinition(d, self.model, NOTHING))
        except Exception as ex:
            print(f"    setting the radius failed: {type(ex).__name__}: {ex}")
            return False

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
        raise WrongShape(f"{name}: no extrusion direction gives the expected volume {expected:.3f} mm^3 "
                           f"and centre of mass {centroid}")


    def revolve(self, sk, expected, centroid, name, cut=False):
        """Revolve feature of 360 degrees about the sketch's construction centerline. The axis is
        taken from the sketch's single centerline; if SolidWorks makes no feature that way, the
        centerline is selected as well (Mark 4, then 16)."""
        fm = self.model.FeatureManager
        cl = next((g for kind, _, g in self.ents if kind == "cline"), None)
        for mark in (None, 4, 16):
            self.model.ClearSelection2(True)
            sk.Select2(False, 0)
            if mark is not None:
                if cl is None:
                    break
                sd = call0(self.model.SelectionManager, "CreateSelectData")
                sd.Mark = mark
                cl.Select4(True, sd)
            f = fm.FeatureRevolve2(True, True, False, cut, False, False, 0, 0, 2 * math.pi, 0.0,
                                   False, False, 0.0, 0.0, 0, 0.0, 0.0, True, True, True)
            if f is None:
                print(f"  revolve{'' if mark is None else f' (centerline Mark {mark})'}: SolidWorks returned no feature")
                continue
            v, c = self.mass()
            print(f"  revolve: volume {v:.3f} mm^3, centre of mass ({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}) mm")
            if abs(v - expected) <= 1e-4 * expected + 1e-3 and math.dist(c, centroid) <= 1e-5 * SIZE + 1e-4:
                f.Name = name
                return v
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(0)
        raise WrongShape(f"{name}: no Revolve feature gives the expected volume {expected:.3f} mm^3")


def _on_arc(P, s, m, t, c):
    """Whether P (on the circle about c) lies on the arc from s through m to t."""
    ang = lambda q: math.atan2(q[1] - c[1], q[0] - c[0])
    a0, am, a1, ap = ang(s), ang(m), ang(t), ang(P)
    ccw = lambda a, b: (b - a) % (2 * math.pi)
    sweep = ccw(a0, a1)
    if ccw(a0, am) > sweep:                      # the arc runs clockwise from s to t
        return ccw(ap, a0) <= ccw(a1, a0) + 1e-9
    return ccw(a0, ap) <= sweep + 1e-9


class WrongShape(RuntimeError):
    """The feature was made, but no direction gives the expected volume and centre of mass."""


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


def build_step(sw, st):
    """Sketch and feature of one step. Entities are first added directly (AddToDB, no snapping);
    if SolidWorks then makes no feature, the sketch is drawn again in the normal mode, which joins
    coincident endpoints. Attempts: AddToDB with the planned relations and dimensions, AddToDB
    with none (the solver cannot move anything), normal mode with none."""
    modes = [(True, True), (True, False), (False, False)]
    for m, (add_to_db, use_plans) in enumerate(modes):
        sk = sw.sketch(st["loops"], add_to_db, plans=st.get("plans") if use_plans else None)
        try:
            if st["op"] == "revolve":
                v = sw.revolve(sk, st["volume"], st["centroid"], st["name"])
            elif sw.direct and "angle" not in st:    # revolved boss or cut about the y axis
                v = sw.revolve(sk, st["volume"], st["centroid"], st["name"], cut=st["op"] == "cut")
            else:
                v = sw.feature(sk, st["op"], st["lo"], st["hi"], st["volume"], st["centroid"], st["name"])
            break
        except (NoFeature, WrongShape) as ex:
            # a planned relation or dimension can pick the wrong entity (picks by coordinates
            # depend on the view scale): the sketch is drawn again without them
            sw.dump_sketch(sk, everything=True)
            if m == len(modes) - 1:
                if isinstance(ex, WrongShape):
                    raise
                raise RuntimeError(f"{st['name']}: SolidWorks makes no feature from the sketch")
            nxt = modes[m + 1]
            print(f"  {'wrong shape' if isinstance(ex, WrongShape) else 'no feature'}; drawing the sketch again ({'AddToDB' if nxt[0] else 'normal mode'}, "
                  f"{'with' if nxt[1] else 'without'} relations and dimensions)")
            sw.model.ClearSelection2(True)
            sk.Select2(False, 0)
            sw.model.Extension.DeleteSelection2(0)
    print(f"{st['name']}: volume {v:.3f} mm^3 (expected {st['volume']:.3f})")


def build():
    sw = SW()
    print(f"connected to SolidWorks; new part created; building {len(FEATURES)} feature(s)")
    for n, st in enumerate(FEATURES, 1):
        build_step(sw, st)
    for st in CUTS:                              # straight cuts across a revolved body
        if abs(st["angle"]) > 1e-9:
            r = sw.turn_y(st["angle"], st["turned"], "Turn cut direction to z")
            if r is None:
                raise RuntimeError(f"{st['name']}: no turn about the axis brings the body to the expected place")
            print(f"turned {math.degrees(st['angle']):+g} deg about the axis " + r)
        sw.direct, sw.helper_ab = True, st["helper"]     # outline in Front Plane coordinates
        build_step(sw, st)
        sw.direct, sw.helper_ab = False, None
        if abs(st["angle"]) > 1e-9:
            r = sw.turn_y(-st["angle"], st["back"], "Turn back")
            if r is None:
                raise RuntimeError(f"{st['name']}: no turn brings the body back")
            print("turned back " + r)
    print(sw.rotate_back())
    if REVOLVES:
        if not sw.turned_back:
            raise RuntimeError("the revolved features need the part in the STL's orientation")
        if REV_TURN is not None:
            r = sw.quarter_turn(REV_TURN, "Turn prism axis to y")
            if r is None:
                raise RuntimeError("no quarter turn brings the prism axis to y")
            print("turned the prism axis to y " + r)
        for g in REVOLVES:
            moved = any(abs(v) > 1e-12 for v in g["shift"])
            if moved:
                print(sw.translate(g["shift"], g["before"], "Move revolve axis to origin"))
            sw.direct = True
            for st in g["steps"]:
                sw.helper_ab = st["helper"]
                build_step(sw, st)
            sw.direct, sw.helper_ab = False, None
            if moved:
                print(sw.translate([-v for v in g["shift"]], g["after"], "Move back"))
        if REV_TURN is not None:
            r = sw.quarter_turn(REV_DONE, "Turn back to STL orientation")
            if r is None:
                raise RuntimeError("no quarter turn brings the part back to the STL's orientation")
            print("turned back to the STL's orientation " + r)
    for cy in CYLINDERS:
        v = sw.face_feature(cy)
        print(f"{cy['name']}: volume {v:.3f} mm^3 (expected {cy['volume']:.3f})")
    for st in XCUTS:                             # cuts along a second axis
        if not sw.turned_back:
            raise RuntimeError("the cuts along a second axis need the part in the STL's orientation")
        if st["axis"] is not None:
            r = sw.turn_about(st["axis"], st["angle"], st["turned"], "Turn cut axis to z")
            if r is None:
                raise RuntimeError(f"{st['name']}: no turn brings the cut axis to z")
            print(f"turned the cut axis to z " + r)
        sw.direct, sw.helper_ab = True, st["helper"]
        build_step(sw, st)
        sw.direct, sw.helper_ab = False, None
        if st["axis"] is not None:
            r = sw.turn_about(st["axis"], -st["angle"], st["back"], "Turn back")
            if r is None:
                raise RuntimeError(f"{st['name']}: no turn brings the body back")
            print("turned back " + r)
    for b in BLENDS:                             # fillets and chamfers, last
        if not sw.turned_back:
            raise RuntimeError("the fillets and chamfers need the part in the STL's orientation")
        v = sw.blend(b)
        print(f"{b['name']}: volume {v:.3f} mm^3 (expected {b['volume']:.3f})")
    v0 = sw.mass()[0]
    print(sw.link_dimensions())
    v1 = sw.mass()[0]
    if abs(v1 - v0) > 1e-4 * v0 + 1e-3:
        print(f"NOTE: linking {[nm for nm, _ in getattr(sw, 'linked', [])]} changed the volume from {v0:.3f} to {v1:.3f} mm^3; "
              "the equations are removed")
        sw.unlink_dimensions()
        v2 = sw.mass()[0]
        if abs(v2 - v0) > 1e-4 * v0 + 1e-3:
            raise RuntimeError(f"after removing the equations the volume is {v2:.3f} mm^3, not {v0:.3f}")
    print(f"dimensions: {sw.dims_ok} added" + (f", {len(sw.dims_failed)} failed: {sw.dims_failed}" if sw.dims_failed else ""))
    print(f"relations: {sw.rels_ok} added" + (f", {len(sw.rels_failed)} failed: {sw.rels_failed}" if sw.rels_failed else ""))
    print("sketch status: " + ", ".join(sw.status))
    print("planned: " + ", ".join(st.get("expect", "?") for st in FEATURES + CUTS + [x for g in REVOLVES for x in g["steps"]] + XCUTS))
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
    rev_turn, rev_groups = revolve_steps(p25, steps, nd)
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
         f"{len(steps) - n_add} Extruded Cut feature(s) along +{AX[k]},"
         + (f" then {sum(len(g['steps']) for g in rev_groups)} Revolved Boss/Cut feature(s) about axes "
            f"along {AX[k]}," if rev_groups else ""),
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
                rot = None if loop.pattern is None else (len(loop.pattern[2]), loop.pattern[1],
                                                         tuple(map(float, loop.pattern[0])))
                E_ = [tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in loop.edges]
                cps = None
                if getattr(loop, "linear", None) is not None:
                    # a sketch linear pattern with its spacing and count dimensions (sw_probe13: the
                    # copies then follow the seed); the seed is the copy at the left (else bottom)
                    # end of the run, the plan starts with it
                    p_, N_ = loop.linear[0], loop.linear[1]
                    s_ = seed_copy(E_, p_, N_)
                    n_ = len(E_)
                    E_ = E_[s_ * p_:] + E_[:s_ * p_]
                    cps = (p_, [((m * p_ - s_ * p_) % n_, m - s_) for m in range(N_) if m != s_])
                ops, ok = plan_outline(E_, rotational=rot, copies=cps)
                plans.append(ops)
                full = full and ok
        L.append(f"        expect=\"{'fully defined' if full else 'under defined'}\",")
        L.append("        plans=" + repr(_round_plans(plans)) + ",")
        L.append("        loops=[")
        for loop in st["loops"]:
            L.append("        " + _loop_code(loop, max(nd, 9)) + ",")
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
        extra = (f'draft={_f(cy["draft"], nd)}, dims=False, ' if cy.get("draft") else "")
        L.append(f'         dia={_f(cy["dia"], nd)}, offset={_f(cy["offset"], nd)}, depth={_f(cy["depth"], 12)}, depth2={_f(cy["depth2"], nd)}, {extra}volume={cy["volume"]:.4f}, '
                 f'centroid=({cy["centroid"][0]:.5f}, {cy["centroid"][1]:.5f}, {cy["centroid"][2]:.5f})),')
    L.append("]")
    L += _revolves_code(p25, rev_turn, rev_groups, nd)
    L.append("CUTS = []")
    L += _xcuts_code(p25, steps, nd)
    L += _blends_code(p25, steps, nd)
    return "\n".join(L) + "\n" + _RUNTIME


def _blends_code(p25, steps, nd):
    """Fillets and chamfers, last (STL frame): the edges picked at their midpoints, with the
    volume and centre of mass expected after each."""
    from .core import mass_props
    from .crosscut import cut_solid
    from .blends import apply_blends
    from .prism25 import cylinder_solid
    bl = p25.get("blends") or []
    L = ["", "# ---- fillets and chamfers: edges picked at their midpoints (STL frame, mm) ----", "BLENDS = ["]
    if not bl:
        return L[:1] + ["BLENDS = []"]
    solid = _prism_solid(p25, steps)
    for ft in p25.get("revolved") or []:
        solid = _apply_revolved(solid, ft)
    for cy in p25.get("cylinders") or []:
        solid = solid.fuse(cylinder_solid(cy)) if cy["op"] == "add" else solid.cut(cylinder_solid(cy))
    for cu in p25.get("crosscuts") or []:
        solid = solid.cut(cut_solid(cu))
    solid = solid.clean()
    t3 = lambda v: "(" + ", ".join(f"{float(x):.6f}" for x in v) + ")"
    for n, b in enumerate(bl, 1):
        solid = apply_blends(solid, [b]).clean()
        vol, com = mass_props(solid)
        L.append(f'    dict(name="{b["kind"].capitalize()} {n} {b["size"]:g} mm", kind="{b["kind"]}", '
                 f'size={_f(b["size"], 9)}, volume={vol:.4f}, centroid={t3(com)},')
        L.append("         points=[" + ", ".join(t3(p) for p in b["points"]) + "]),")
    L.append("]")
    return L


XCUT_TURN = {   # cut axis -> (model axis turned about, signed angle) bringing it to +z
    2: (None, 0.0),
    0: (1, -math.pi / 2),        # about y: x -> z
    1: (0, math.pi / 2),         # about x: y -> z
}


def _rot(m, ang):
    c, s_ = math.cos(ang), math.sin(ang)
    if m is None:
        return np.eye(3)
    if m == 0:
        return np.array([[1, 0, 0], [0, c, -s_], [0, s_, c]])
    if m == 1:
        return np.array([[c, 0, s_], [0, 1, 0], [-s_, 0, c]])
    return np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])


def _xcuts_code(p25, steps, nd):
    """Cuts along a second model axis, after the off-axis cylinders (STL frame): the body is
    turned so the cut axis is +z, the outline sketched on the Front Plane and cut from lo to hi
    along z, and the body turned back; expected volumes and centres of mass from CadQuery."""
    from .core import mass_props
    from .crosscut import plane_frame, cut_solid
    from .prism25 import Loop, cylinder_solid
    from .sketchplan import plan_outline
    xc = p25.get("crosscuts") or []
    L = ["", "# ---- cuts along a second axis: the body is turned so the cut axis is z ----", "XCUTS = ["]
    if not xc:
        return L[:1] + ["XCUTS = []"]
    solid = _prism_solid(p25, steps)
    for ft in p25.get("revolved") or []:
        solid = _apply_revolved(solid, ft)
    for cy in p25.get("cylinders") or []:
        solid = solid.fuse(cylinder_solid(cy)) if cy["op"] == "add" else solid.cut(cylinder_solid(cy))
    solid = solid.clean()
    t3 = lambda v: "(" + ", ".join(f"{float(x):.5f}" for x in v) + ")"
    for n, cu in enumerate(xc, 1):
        m, ang = XCUT_TURN[cu["j"]]
        R = _rot(m, ang)
        w, h = plane_frame(cu["j"])
        com0 = np.array(mass_props(solid)[1])
        solid = solid.cut(cut_solid(cu)).clean()
        vol, com = mass_props(solid)
        com = np.array(com)
        to_sk = lambda q: np.array((R @ (float(q[0]) * w + float(q[1]) * h))[:2]) + 0.0
        loops = []
        for Lc in cu["loops"]:
            es = []
            for e in Lc.edges:
                if e[0] == "circle":
                    es.append(("circle", to_sk(e[1]), e[2]))
                else:
                    es.append((e[0],) + tuple(to_sk(q) for q in e[1:]))
            loops.append(Loop(es, Lc.area))
        plans, full = [], True
        for Lc in loops:
            if len(Lc.edges) == 1:
                plans.append([])
                continue
            E_ = [tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in Lc.edges]
            ops, ok = plan_outline(E_)
            plans.append(ops)
            full = full and ok
        P = np.array([q for Lc in loops for e in Lc.edges for q in e[1:] if np.ndim(q) == 1])
        L.append("    dict(")
        L.append(f'        name="Cut {n} along {AX[cu["j"]]} {cu["a"]:.4g} to {cu["b"]:.4g}", op="cut", '
                 f"axis={m}, angle={ang:.12f}, lo={_f(cu['lo'], nd)}, hi={_f(cu['hi'], nd)}, volume={vol:.4f},")
        L.append(f"        turned={t3(R @ com0)}, centroid={t3(R @ com)}, back={t3(com)}, holes=False,")
        L.append(f"        helper=({_f(P[:, 0].min() - 10.0, 3)}, {_f(P[:, 1].min() - 10.0, 3)}),")
        L.append(f"        expect=\"{'fully defined' if full else 'under defined'}\",")
        L.append("        plans=" + repr(_round_plans(plans)) + ",")
        L.append("        loops=[")
        for Lp in loops:
            L.append("        " + _loop_code(Lp, max(nd, 9)) + ",")
        L.append("        ],")
        L.append("    ),")
    L.append("]")
    return L


def _revolves_code(p25, turn, groups, nd):
    """Script data for the revolved features (see revolve_steps)."""
    k = p25["axis"]
    t3 = lambda v: "(" + ", ".join(f"{float(x):.5f}" for x in v) + ")"
    L = ["",
         f"# ---- revolved features about axes along {AX[k]}: the body is turned with {AX[k]} along y "
         "(if it is not),",
         "# moved so that each axis is the y axis, the profile (a = distance from the axis, b = level) "
         "is sketched",
         "# on the Front Plane and revolved about a vertical centerline at x = 0, and the body is moved "
         "back ----"]
    if not groups:
        return L + ["REV_TURN = None", "REV_DONE = None", "REVOLVES = []"]
    Rk = np.asarray(REVOLVE_FRAME[k][0], float)
    done = Rk.T @ groups[-1]["after"]
    L.append(f"REV_TURN = {None if turn is None else t3(turn)}   # centre of mass once {AX[k]} is along y (mm)")
    L.append(f"REV_DONE = {t3(done)}   # centre of mass after the revolved features, STL frame (mm)")
    L.append("REVOLVES = [")
    for g in groups:
        L.append(f"    dict(shift={t3(g['shift_v'])}, before={t3(g['before'])}, after={t3(g['after'])}, steps=[")
        for st in g["steps"]:
            L.append("        dict(")
            L.append(f'            name="{st["name"]}", op="{st["op"]}", volume={st["volume"]:.4f}, '
                     f"centroid={t3(st['centroid'])}, helper=({_f(st['helper'][0], 3)}, {_f(st['helper'][1], 3)}),")
            L.append(f"            expect=\"{'fully defined' if st['full'] else 'under defined'}\",")
            L.append("            plans=" + repr(_round_plans(st["plans"])) + ",")
            L.append("            loops=[")
            for Lp in st["loops"]:
                L.append("            " + _loop_code(Lp, max(nd, 9)) + ",")
            L.append("            ],")
            L.append("        ),")
        L.append("    ]),")
    L.append("]")
    return L


def _cuts_code(rv, shape, R, nd):
    """Straight cuts across a revolved body (build frame: axis along y through the origin). For
    each cut the body is turned about y by the angle that brings the cut direction to +z, the
    outline is sketched on the Front Plane (x = -w, y = v) and cut from level lo to hi along z,
    and the body is turned back; expected volumes and centres of mass from CadQuery."""
    from .core import mass_props
    from .prism25 import Loop
    from .sketchplan import plan_outline
    cuts = rv.get("cuts") or []
    L = ["", "# ---- straight cuts across the revolved body: turned about the axis so the cut runs along z ----",
         "CUTS = ["]
    if not cuts:
        return L[:1] + ["CUTS = []"]
    import cadquery as cq
    from .revolve import _solid, _cut_wire
    c, d, e1 = (np.asarray(rv[k_], float) for k_ in ("c", "d", "e1"))
    Rb = np.asarray(R, float)
    solid = _solid(rv["loops"], c, d, e1)
    reach = 1.2 * max(abs(float(q[0])) for L_ in rv["loops"] for e in L_.edges for q in e[1:])
    Ry = lambda a: np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    for n, cu in enumerate(cuts, 1):
        Db = Rb @ np.asarray(cu["D"], float)
        psi = float(np.arctan2(-Db[0], Db[2]))       # Ry(psi) Db = +z
        com0 = Rb @ np.array(mass_props(solid)[1])
        D, w = np.asarray(cu["D"], float), np.asarray(cu["w"], float)
        body = None
        for Lc in cu["loops"]:
            f = cq.Face.makeFromWires(_cut_wire(Lc, c, w, d), [])
            if cu["both"]:
                piece = cq.Solid.extrudeLinear(f.translate(cq.Vector(*(-reach * D))), cq.Vector(*(2 * reach * D)))
            else:
                piece = cq.Solid.extrudeLinear(f, cq.Vector(*(reach * D)))
            body = piece if body is None else body.fuse(piece)
        solid = solid.cut(body).clean()
        vol, com = mass_props(solid)
        comb = Rb @ np.array(com)
        # the outline in the turned frame: w turns to -x, the axis is y
        wt = Ry(psi) @ (Rb @ w)
        sx = float(np.sign(wt[0]))
        loops = [Loop([(e[0],) + tuple(np.array([sx * q[0], q[1]], float) for q in e[1:]) for e in Lc.edges], Lc.area)
                 for Lc in cu["loops"]]
        plans, full = [], True
        for Lc in loops:
            E_ = [tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in Lc.edges]
            ops, ok = plan_outline(E_)
            plans.append(ops)
            full = full and ok
        lo, hi = (-reach, reach) if cu["both"] else (0.0, reach)
        P = np.array([q for Lc in loops for e in Lc.edges for q in e[1:]])
        t3 = lambda v: "(" + ", ".join(f"{float(x):.5f}" for x in v) + ")"
        L.append("    dict(")
        L.append(f'        name="Cut {n} across the axis ({", ".join(f"{x + 0:g}" for x in np.round(D, 4))})", op="cut", '
                 f"angle={psi:.12f}, lo={_f(lo, nd)}, hi={_f(hi, nd)}, volume={vol:.4f},")
        L.append(f"        turned={t3(Ry(psi) @ com0)}, centroid={t3(Ry(psi) @ comb)}, back={t3(comb)}, holes=False,")
        L.append(f"        helper=({_f(P[:, 0].min() - 10.0, 3)}, {_f(P[:, 1].min() - 10.0, 3)}),")
        L.append(f"        expect=\"{'fully defined' if full else 'under defined'}\",")
        L.append("        plans=" + repr(_round_plans(plans)) + ",")
        L.append("        loops=[")
        for Lp in loops:
            L.append("        " + _loop_code(Lp, max(nd, 9)) + ",")
        L.append("        ],")
        L.append("    ),")
    L.append("]")
    return L


REVOLVE_FRAME = {      # STL frame -> build frame turning the revolve axis onto +y, a quarter turn
    1: ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (0.0, 0.0, 0.0)),
    0: ([[0, -1, 0], [1, 0, 0], [0, 0, 1]], (0.0, 0.0, -math.pi / 2)),   # +90 deg about z
    2: ([[1, 0, 0], [0, 0, 1], [0, -1, 0]], (math.pi / 2, 0.0, 0.0)),    # -90 deg about x
}


def solidworks_revolve_script(rv, shape, name, nd=6):
    """SolidWorks script for a body of revolution: the profile sketched on the Front Plane with the
    axis as a vertical construction centerline through the origin, one Revolve feature of 360
    degrees, then the body turned to the STL's orientation. Needs the axis along a model axis and
    through the origin."""
    from .sketchplan import plan_outline
    k = rv["axis"]
    c, d = np.asarray(rv["c"], float), np.asarray(rv["d"], float)
    size = float(np.linalg.norm(np.ptp(np.array([np.asarray(q, float) for L in rv["loops"]
                                                    for e in L.edges for q in e[1:]]), axis=0))) + 1.0
    if np.linalg.norm(c) > 1e-6 * size or abs(abs(d[k]) - 1.0) > 1e-12:
        raise ValueError("the revolve axis is not a model axis through the origin")
    R, back = REVOLVE_FRAME[k]
    Rt = np.asarray(R, float).T
    ai = int(np.argmax(np.abs(Rt @ np.array([1.0, 0.0, 0.0]))))     # model axis drawn as sketch x
    sgn = float(np.sign((Rt @ np.array([1.0, 0.0, 0.0]))[ai]))
    # profile (u, v) -> sketch (a, b) = model coordinates (axis ai, axis k); u is drawn on the
    # side where the build frame's x is positive
    loops = []
    from .prism25 import Loop
    for L in rv["loops"]:
        es = [(e[0],) + tuple(np.array([sgn * q[0], q[1]], float) for q in e[1:]) for e in L.edges]
        loops.append(Loop(es, L.area))
    if rv.get("cuts"):                           # the revolve step: the body before the cuts
        from .revolve import _solid
        vol, cen = mass_props(_solid(rv["loops"], np.asarray(rv["c"], float), np.asarray(rv["d"], float),
                                     np.asarray(rv["e1"], float)))
    else:
        vol, cen = mass_props(shape)
    cb = np.asarray(R, float) @ np.asarray(cen)
    allp = np.array([np.asarray(q, float) for L in loops for e in L.edges for q in e[1:]])
    helper = allp.min(axis=0) - 10.0
    plans, full = [], True
    for n_, L in enumerate(sorted(loops, key=lambda L: -abs(L.area))):
        E_ = [tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in L.edges]
        ops, ok = plan_outline(E_, revolve_axis=(n_ == 0))
        plans.append(ops)
        full = full and ok
    loops = sorted(loops, key=lambda L: -abs(L.area))
    L = [f'"""{name}: SolidWorks feature tree generated by stl2cad (exact mode, revolve).',
         "",
         f"Builds the part in SolidWorks through its API: one Revolve feature of 360 degrees about the "
         f"{AX[k]} axis,",
         "from a profile sketched on the Front Plane with the axis as a construction centerline; "
         + ("the part is built with the axis along y and turned back by a Move/Copy Body feature. "
            if k != 1 else "") + "Units: mm.",
         "The part volume and centre of mass are checked against stl2cad's values.",
         "",
         "Needs: Windows, SolidWorks, Python with pywin32 (pip install pywin32).",
         f"Run:   python {name}.py [C:\\\\path\\\\part.SLDPRT]",
         '"""',
         "",
         f"AXIS_I, AXIS_J = {ai}, {k}           # model axes of the sketch coordinates (a, b)",
         'PLANE_NAME = "Front Plane"',
         "PLANE_INDEX = 0                 # its position among the template planes",
         f"ROT = {[list(map(float, r)) for r in R]}   # STL frame -> build frame (revolve axis {AX[k]} -> y)",
         f"ROT_BACK = {tuple(back)}   # rotation (rad) about x, y, z back to the STL frame",
         f"FINAL_CENTROID = {tuple(round(float(x), 5) for x in shape.Center().toTuple())}   # STL frame, mm",
         f"SIZE = {_f(size, 3)}                  # part diagonal (mm), scale of the centre-of-mass check",
         f"HELPER_AB = ({_f(helper[0], 3)}, {_f(helper[1], 3)})   # helper circle (a, b), outside the part",
         "",
         "# ---- features: the revolve profile (a, b in mm), expected volume (mm^3) ----",
         "FEATURES = [",
         "    dict(",
         f'        name="Revolve 1 about {AX[k]}", op="revolve", lo=0.0, hi=0.0, '
         f"volume={vol:.4f}, centroid=({cb[0]:.5f}, {cb[1]:.5f}, {cb[2]:.5f}), holes=False,",
         f"        expect=\"{'fully defined' if full else 'under defined'}\",",
         "        plans=" + repr(_round_plans(plans)) + ",",
         "        loops=["]
    for Lp in loops:
        L.append("        " + _loop_code(Lp, max(nd, 9)) + ",")
    L += ["        ],", "    ),", "]", "", "CYLINDERS = []", "REV_TURN = None", "REV_DONE = None", "REVOLVES = []", "XCUTS = []",
          "BLENDS = []"]
    L += _cuts_code(rv, shape, R, nd)
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
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.path.join(os.path.dirname(__file__), "swmock"),
                                                    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))]),
               SWMOCK_CUT_REVERSED="1" if cut_reversed else "0",
               SWMOCK_BOSS_REVERSED="1" if boss_reversed else "0")
    r = subprocess.run([sys.executable, py, out], capture_output=True, text=True, timeout=timeout, env=env)
    if r.returncode != 0 or not os.path.exists(out):
        tail = (r.stderr or r.stdout).strip().splitlines()[-1:] or ["no output"]
        return False, "stand-in run failed: " + tail[0][:200]
    import json
    import re
    s = cq.importers.importStep(out).val()
    dv = abs(mass_props(s)[0] - mass_props(ref_shape)[0]) / mass_props(ref_shape)[0]
    m = re.search(r"dimensions: (\d+) added(, (\d+) failed)?", r.stdout)
    m2 = re.search(r"relations: (\d+) added(, (\d+) failed)?", r.stdout)
    if m is None or m.group(3) or m2 is None or m2.group(3):
        return False, ("dimensions or relations failed on the API stand-in: "
                       + (m.group(0) if m else "no report") + "; " + (m2.group(0) if m2 else "no report"))
    dims = json.load(open(out + ".dims.json"))["dims"]
    dist = max(_max_distance(ref_shape, s, tol), _max_distance(s, ref_shape, tol))
    ok = dv < 1e-4 and dist < 0.5 * tol
    return ok, f"API stand-in: volume {100 * dv:.4f} %, max distance {dist:.4f} mm, {len(dims)} dimension(s), {m2.group(1)} relation(s)"
