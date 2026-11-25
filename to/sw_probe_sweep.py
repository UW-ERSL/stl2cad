"""Probe: swept members (path of lines and tangent arcs, section normal to the path).
Run on the Windows machine with SolidWorks:  python sw_probe_sweep.py  > probe_sweep.log
and send the log. Each step prints what SolidWorks returned; nothing is assumed.

Path (Front Plane, mm): line (0,0)-(20,0), arc centre (20,10) to (30,10), line (30,10)-(30,30). Length 50.71.
Section at the start: rectangle 4 (in-plane, across the path, i.e. along y) x 2 (along z), centred on the path.
Expected swept volume 4 x 2 x 50.71 = 405.66 mm^3; centre of mass of the solid (z = 0).
Steps:
  1. path sketch on the Front Plane;
  2. a reference plane normal to the path at its start (InsertRefPlane: perpendicular to the first line,
     coincident with its start point), with several selection variants;
  3. a sketch on that plane: where do its entities land? (its ModelToSketchTransform, and the plane's
     sketch axes found by drawing a point and reading it back);
  4. the rectangle drawn there (several axis guesses), then a Swept Boss (InsertProtrusionSwept4; profile
     Mark 1, path Mark 4), volume and centre of mass compared;
  5. a second section (8 x 2) on a plane normal to the path at its end, and a Loft from the start section
     to the end section with the path as centerline (InsertProtrusionBlend2; profiles Mark 1, centerline
     Mark 4). Expected volume for a linear taper 4 -> 8: 6 x 2 x 50.71 = 608.5 mm^3;
  6. an ellipse section (CreateEllipse) instead of the rectangle, swept.
"""
import math
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
PERP, COINC = 2, 4                 # swRefPlaneReferenceConstraints_e: Perpendicular, Coincident


def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


def mass(model):
    model.ForceRebuild3(False)
    mp = call0(model.Extension, "CreateMassProperty")
    if mp is None:
        return 0.0, (0, 0, 0)
    return round(mp.Volume / MM ** 3, 3), tuple(round(c / MM, 3) for c in mp.CenterOfMass)


def step(label, fn):
    try:
        r = fn()
        print(f"{label}: {r}")
        return r
    except Exception as ex:
        print(f"{label}: FAILED {type(ex).__name__}: {ex}")
        return None


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True
try:
    app.SetUserPreferenceToggle(10, False)          # no Modify dialog on dimensions
except Exception:
    pass


def new_part():
    return call0(app, "NewPart")


def select_plane(model, name):
    model.ClearSelection2(True)
    return model.Extension.SelectByID2(name, "PLANE", 0, 0, 0, False, 0, NOTHING, 0)


def path_sketch(model):
    """Front Plane: line, tangent arc, line (helper circle first: lines need another entity)."""
    select_plane(model, "Front Plane")
    sm = model.SketchManager
    sm.InsertSketch(True)
    sm.AddToDB = True
    h = sm.CreateCircleByRadius(-0.02, -0.02, 0.0, 0.0005)
    a = sm.CreateLine(0, 0, 0, 0.020, 0, 0)
    b = sm.CreateArc(0.020, 0.010, 0, 0.020, 0, 0, 0.030, 0.010, 0, 1)
    c = sm.CreateLine(0.030, 0.010, 0, 0.030, 0.030, 0)
    sm.AddToDB = False
    model.ClearSelection2(True)
    if h is not None and h.Select4(False, NOTHING):
        model.EditDelete()
    sm.InsertSketch(True)
    f = model.FeatureByPositionReverse(0)
    print(f"  path sketch: line {a is not None}, arc {b is not None}, line {c is not None}; feature {f.Name}")
    return f


def normal_plane(model, path, at_end=False):
    """Reference plane normal to the path at its start (or end). Variants of selecting the references."""
    seg_pt = (0.010, 0.0, 0.0) if not at_end else (0.030, 0.020, 0.0)      # a point on the first / last line
    end_pt = (0.0, 0.0, 0.0) if not at_end else (0.030, 0.030, 0.0)
    ext, fm = model.Extension, model.FeatureManager
    variants = [
        ("EXTSKETCHSEGMENT mark 0 + EXTSKETCHPOINT mark 1",
         lambda: ext.SelectByID2("", "EXTSKETCHSEGMENT", *seg_pt, False, 0, NOTHING, 0)
         and ext.SelectByID2("", "EXTSKETCHPOINT", *end_pt, True, 1, NOTHING, 0)),
        ("EXTSKETCHSEGMENT mark 0 + EXTSKETCHPOINT mark 0",
         lambda: ext.SelectByID2("", "EXTSKETCHSEGMENT", *seg_pt, False, 0, NOTHING, 0)
         and ext.SelectByID2("", "EXTSKETCHPOINT", *end_pt, True, 0, NOTHING, 0)),
        ("EDGE-style SKETCHSEGMENT names: Line1@<sketch> mark 0 + Point1@<sketch> mark 1",
         lambda: ext.SelectByID2(f"Line{3 if at_end else 1}@{path.Name}", "EXTSKETCHSEGMENT", 0, 0, 0, False, 0, NOTHING, 0)
         and ext.SelectByID2("", "EXTSKETCHPOINT", *end_pt, True, 1, NOTHING, 0)),
    ]
    for label, sel in variants:
        model.ClearSelection2(True)
        try:
            ok = sel()
            n = model.SelectionManager.GetSelectedObjectCount2(-1)
            p = fm.InsertRefPlane(PERP, 0, COINC, 0, 0, 0)
            print(f"  plane ({label}): selected {ok} ({n}); plane {'yes' if p is not None else 'NO'}")
            if p is not None:
                f = model.FeatureByPositionReverse(0)
                return f
        except Exception as ex:
            print(f"  plane ({label}): FAILED {type(ex).__name__}: {ex}")
    return None


def sketch_on(model, plane_feature, draw):
    """Open a sketch on the plane feature, draw, close; report where the sketch is."""
    model.ClearSelection2(True)
    plane_feature.Select2(False, 0)
    sm = model.SketchManager
    sm.InsertSketch(True)
    sk = call0(sm, "ActiveSketch")
    print(f"  ActiveSketch: {'object' if sk is not None else None}")
    T = None
    if sk is not None:
        try:
            T = sk.ModelToSketchTransform
            print(f"  ModelToSketchTransform.ArrayData: {[round(v, 6) for v in T.ArrayData]}")
        except Exception as ex:
            print(f"  ModelToSketchTransform: FAILED {type(ex).__name__}: {ex}")
    sm.AddToDB = True
    h = sm.CreateCircleByRadius(-0.02, -0.02, 0.0, 0.0005)
    ents = draw(sm)
    sm.AddToDB = False
    model.ClearSelection2(True)
    if h is not None and h.Select4(False, NOTHING):
        model.EditDelete()
    # where did the first entity land? read its start point back in model coordinates if possible
    try:
        e = ents[0]
        p = e.GetStartPoint2()
        print(f"  first entity start point (sketch coords, mm): ({p.X / MM:.3f}, {p.Y / MM:.3f}, {p.Z / MM:.3f})")
    except Exception as ex:
        print(f"  first entity start point: {type(ex).__name__}: {ex}")
    sm.InsertSketch(True)
    f = model.FeatureByPositionReverse(0)
    print(f"  section sketch: {f.Name}, entities {[x is not None for x in ents]}")
    return f


def rect(w, d, axes):
    """Rectangle w (along the path's in-plane normal) x d (along model z), drawn in sketch coordinates
    under an axis guess: axes = (sx, sy) where sketch x is +-'n' or +-'z' etc."""
    def draw(sm):
        def P(a, b):           # a: along n (in-plane normal), b: along z
            m = {"n": a, "-n": -a, "z": b, "-z": -b}
            return m[axes[0]] * MM, m[axes[1]] * MM
        pts = [P(w / 2, d / 2), P(-w / 2, d / 2), P(-w / 2, -d / 2), P(w / 2, -d / 2)]
        return [sm.CreateLine(*pts[k], 0, *pts[(k + 1) % 4], 0) for k in range(4)]
    return draw


def sweep(model, profile, path):
    fm = model.FeatureManager
    model.ClearSelection2(True)
    sd = call0(model.SelectionManager, "CreateSelectData")
    sd.Mark = 1
    ok1 = profile.Select2(False, 1)
    ok2 = path.Select2(True, 4)
    print(f"  selected profile {ok1}, path {ok2}")
    tries = [
        ("InsertProtrusionSwept4 (17 args)", lambda: fm.InsertProtrusionSwept4(
            False, False, 0, False, False, 0, 0, False, 0.0, 0.0, 0, 0, True, True, True, 0.0, True)),
        ("InsertProtrusionSwept3 (16 args)", lambda: fm.InsertProtrusionSwept3(
            False, False, 0, False, False, 0, 0, False, 0.0, 0.0, 0, 0, True, True, True, 0.0)),
    ]
    for label, mk in tries:
        try:
            f = mk()
            print(f"  {label}: {'feature' if f is not None else 'None'}; mass {mass(model)}")
            if f is not None:
                return f
        except Exception as ex:
            print(f"  {label}: FAILED {type(ex).__name__}: {ex}")
        model.ClearSelection2(True)
        profile.Select2(False, 1)
        path.Select2(True, 4)
    return None


def loft(model, s0, s1, path):
    fm = model.FeatureManager
    model.ClearSelection2(True)
    s0.Select2(False, 1); s1.Select2(True, 1); path.Select2(True, 4)
    tries = [
        ("InsertProtrusionBlend2 (18 args)", lambda: fm.InsertProtrusionBlend2(
            False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False, 0.0, 0.0, 0, True, True, True, 0)),
        ("InsertProtrusionBlend (17 args)", lambda: fm.InsertProtrusionBlend(
            False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False, 0.0, 0.0, 0, True, True, True)),
    ]
    for label, mk in tries:
        try:
            f = mk()
            print(f"  {label}: {'feature' if f is not None else 'None'}; mass {mass(model)}")
            if f is not None:
                return f
        except Exception as ex:
            print(f"  {label}: FAILED {type(ex).__name__}: {ex}")
        model.ClearSelection2(True)
        s0.Select2(False, 1); s1.Select2(True, 1); path.Select2(True, 4)
    return None


print("== A. rectangle sections under four sketch-axis guesses, Swept Boss")
for axes in (("n", "z"), ("z", "n"), ("-n", "z"), ("z", "-n")):
    print(f"-- guess: sketch x = {axes[0]}, sketch y = {axes[1]}")
    try:
        m = new_part()
        path = path_sketch(m)
        pl = normal_plane(m, path)
        if pl is None:
            print("  no plane; next guess"); continue
        prof = sketch_on(m, pl, rect(4, 2, axes))
        f = sweep(m, prof, path)
        if f is not None:
            try:
                box = m.GetBodies2(0, True)[0].GetBodyBox()
                print(f"  body box (mm): {[round(v / MM, 3) for v in box]}   (expected x -2..32 or so, z -1..1)")
            except Exception as ex:
                print(f"  body box: {type(ex).__name__}: {ex}")
    except Exception:
        traceback.print_exc()

print("== B. loft 4x2 -> 8x2 along the path as centerline")
try:
    m = new_part()
    path = path_sketch(m)
    p0, p1 = normal_plane(m, path), normal_plane(m, path, at_end=True)
    if p0 is not None and p1 is not None:
        s0 = sketch_on(m, p0, rect(4, 2, ("n", "z")))
        s1 = sketch_on(m, p1, rect(8, 2, ("n", "z")))
        loft(m, s0, s1, path)
except Exception:
    traceback.print_exc()

print("== C. ellipse section (semi-axes 2 x 1), Swept Boss; expected volume pi x 2 x 1 x 50.71 = 318.6 mm^3")
try:
    m = new_part()
    path = path_sketch(m)
    pl = normal_plane(m, path)
    if pl is not None:
        prof = sketch_on(m, pl, lambda sm: [sm.CreateEllipse(0, 0, 0, 0.002, 0, 0, 0, 0.001, 0)])
        sweep(m, prof, path)
except Exception:
    traceback.print_exc()
print("done")
