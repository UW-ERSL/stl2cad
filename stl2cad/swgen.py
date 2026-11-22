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
        return f"('circle', {P(e[1])}, {_f(e[2], nd)})"
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
        sm.DisplayWhenAdded = not add_to_db
        # Over pywin32, SolidWorks 2025 creates no line in a sketch until another entity exists;
        # a small helper circle is drawn first and removed (or made construction) at the end.
        p0 = self.to_sketch(*HELPER_AB)          # away from the part, so selections do not hit it
        helper = sm.CreateCircleByRadius(p0[0], p0[1], 0.0, HELPER_R * MM)    # may return None
        for loop in loops:
            for e in expand(loop):
                if e[0] == "line":
                    (x1, y1), (x2, y2) = self.to_sketch(*e[1]), self.to_sketch(*e[2])
                    make(sm.CreateLine, x1, y1, 0.0, x2, y2, 0.0)
                elif e[0] == "circle":
                    (xc, yc) = self.to_sketch(*e[1])
                    make(sm.CreateCircleByRadius, xc, yc, 0.0, e[2] * MM)
                else:
                    s, m, t = (self.to_sketch(*p) for p in e[1:])
                    xc, yc = circle_centre(s, m, t)
                    turn = (m[0] - s[0]) * (t[1] - s[1]) - (m[1] - s[1]) * (t[0] - s[0])
                    make(sm.CreateArc, xc, yc, 0.0, s[0], s[1], 0.0, t[0], t[1], 0.0, 1 if turn > 0 else -1)
        sm.AddToDB = False
        sm.DisplayWhenAdded = True
        if helper is not None:
            self.remove_helper(helper)
        self.dimension_holes(loops)              # every full circle: diameter and position
        for plan in plans or []:
            self.apply_plan(plan)
        sm.InsertSketch(True)                # close the sketch
        self.report_status()
        return self.model.FeatureByPositionReverse(0)

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
        for loop in loops:
            e = expand(loop)[0]
            if e[0] != "circle":
                continue
            cx, cy = self.to_sketch(*e[1])
            r = e[2] * MM
            tag = f"hole at ({e[1][0]:g}, {e[1][1]:g})"
            try:
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

    # ---- outline relations and dimensions -------------------------------------------------------
    def _select(self, typ, a, b, append=False):
        """Select a sketch segment or point at model point (a, b) (mm) of the build frame."""
        x, y = self.to_sketch(a, b)
        return self.model.Extension.SelectByID2("", typ, x, y, 0.0, append, 0, NOTHING, 0)

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

    def apply_plan(self, plan):
        """Relations and dimensions planned by stl2cad for one outline, in order."""
        for op in plan:
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
        else:
            self.dims_ok += 1

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
        return ("NOTE: part left in the build frame (extrusion axis along z); no quarter turn brought it "
                "to the STL's orientation")

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
            if abs(v - expected) <= 1e-4 * expected + 1e-3 and dc <= 1e-4 * SIZE + 1e-4:
                f.Name = name
                return v
            self.model.ClearSelection2(True)
            f.Select2(False, 0)
            self.model.Extension.DeleteSelection2(0)       # keeps the sketch for the next try
        if nothing == len(tries):
            raise NoFeature(name)
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
        for add_to_db in (True, False):
            sk = sw.sketch(st["loops"], add_to_db, plans=st.get("plans"))
            try:
                v = sw.feature(sk, st["op"], st["lo"], st["hi"], st["volume"], st["centroid"], st["name"])
                break
            except NoFeature:
                if not add_to_db:
                    raise RuntimeError(f"{st['name']}: SolidWorks makes no feature from the sketch")
                print("  no feature from the directly added sketch; drawing it again in normal mode")
                sw.model.ClearSelection2(True)
                sk.Select2(False, 0)
                sw.model.Extension.DeleteSelection2(0)
        print(f"{st['name']}: volume {v:.3f} mm^3 (expected {st['volume']:.3f})")
    print(sw.rotate_back())
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
            if loop.pattern is None and not (len(loop.edges) == 1 and loop.edges[0][0] == "circle"):
                ops, ok = plan_outline([tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in loop.edges])
                plans.append(ops)
                full = full and ok
            elif loop.pattern is not None:
                full = False
        L.append(f"        expect=\"{'fully defined' if full else 'under defined'}\",")
        L.append("        plans=" + repr(_round_plans(plans)) + ",")
        L.append("        loops=[")
        for loop in st["loops"]:
            L.append("        " + _loop_code(loop, nd) + ",")
        L.append("        ],")
        L.append("    ),")
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
