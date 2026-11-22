"""Seventh probe: features that do not lie along the Front plane normal.
Run on the Windows machine with SolidWorks:  python sw_probe7.py

Each case makes a new part with a 40 x 10 x 10 mm block (x 10..50, y 20..30, z 0..10 mm) and adds
a boss of radius 2 mm and height 5 mm by one method. From the volume and centre of mass before and
after, the probe prints where the boss went (its centre of mass) in model mm.

A. Sketch on the block's end face x = 50 (normal +x), selected with SelectByID2 "FACE":
   A1 SketchManager.InsertSketch(True), A2 ModelDoc2.InsertSketch2(True).
   Circle centre at sketch (0, 0), (5, 0) and (0, 5) mm; the three boss positions give the sketch
   origin and axes on the face. Also prints whether an active sketch exists and its transform.
   Expected if the sketch lies on the face: boss centres at x = 52.5 mm.
B. Separate body: circle on Front, extruded with Merge off (a second body), turned 90 deg about y
   and moved with InsertMoveCopyBody2, then joined with InsertCombineFeature (add).
   Expected: one body, boss centre near (52.5, 25, 5) mm.
"""
import math

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
RECT = [(10, 20), (50, 20), (50, 30), (10, 30), (10, 20)]


def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True


def extrude(model, depth, merge=True, reverse=False):
    f = model.FeatureByPositionReverse(0)
    model.ClearSelection2(True)
    f.Select2(False, 0)
    return model.FeatureManager.FeatureExtrusion3(
        True, reverse, False, 0, 0, depth, 0.0, False, False, False, False, 0.0, 0.0,
        False, False, False, False, merge, True, True, 0, 0.0, False)


def block():
    model = call0(app, "NewPart")
    sm = model.SketchManager
    c = sm.CreateCircleByRadius(0.005, 0.005, 0.0, 0.0005)
    for (x1, y1), (x2, y2) in zip(RECT[:-1], RECT[1:]):
        sm.CreateLine(x1 * MM, y1 * MM, 0.0, x2 * MM, y2 * MM, 0.0)
    c.ConstructionGeometry = True
    sm.InsertSketch(True)
    extrude(model, 0.010)
    return model


def mass(model):
    model.ForceRebuild3(False)
    mp = call0(model.Extension, "CreateMassProperty")
    return mp.Volume / MM ** 3, [v / MM for v in mp.CenterOfMass]


def boss_centre(before, after):
    (v0, c0), (v1, c1) = before, after
    dv = v1 - v0
    if abs(dv) < 1e-6:
        return "no volume added"
    return f"added {dv:.1f} mm^3 at " + str(tuple(round((c1[m] * v1 - c0[m] * v0) / dv, 2) for m in range(3)))


def bodies(model):
    b = model.GetBodies2(0, True)
    return 0 if b is None else len(b)


def case_a(label, open_sketch, a, b):
    model = block()
    before = mass(model)
    try:
        model.ClearSelection2(True)
        ok = model.Extension.SelectByID2("", "FACE", 50 * MM, 25 * MM, 5 * MM, False, 0, NOTHING, 0)
        open_sketch(model)
        sk = call0(model, "GetActiveSketch2")
        info = "no active sketch"
        if sk is not None:
            try:
                t = call0(sk, "ModelToSketchTransform")
                info = "transform " + str([round(v, 4) for v in t.ArrayData])
            except Exception as e:
                info = f"active sketch, transform failed ({type(e).__name__})"
        model.SketchManager.CreateCircleByRadius(a * MM, b * MM, 0.0, 2 * MM)
        model.SketchManager.InsertSketch(True)
        f = extrude(model, 0.005)
        print(f"{label} circle at ({a}, {b}): face selected {ok}; {info}; feature {'yes' if f else 'NO'}; "
              f"{boss_centre(before, mass(model))}")
    except Exception as e:
        print(f"{label} circle at ({a}, {b}) failed: {type(e).__name__}: {e}")


for label, opener in (("A1 SketchManager.InsertSketch", lambda m: m.SketchManager.InsertSketch(True)),
                      ("A2 ModelDoc2.InsertSketch2   ", lambda m: m.InsertSketch2(True))):
    for a, b in ((0, 0), (5, 0), (0, 5)):
        case_a(label, opener, a, b)

# ---- B: separate body, moved and combined
model = block()
before = mass(model)
try:
    sm = model.SketchManager
    sm.CreateCircleByRadius(0.0, 0.0, 0.0, 2 * MM)       # boss on Front, z 0..5, about the origin
    sm.InsertSketch(True)
    f = extrude(model, 0.005, merge=False)
    print(f"B  separate boss body: feature {'yes' if f else 'NO'}, bodies {bodies(model)}")
    bs = model.GetBodies2(0, True)
    tool = bs[-1]
    model.ClearSelection2(True)
    sd = call0(model.SelectionManager, "CreateSelectData")
    sd.Mark = 1
    tool.Select2(False, sd)
    # turn 90 deg about y (z -> +x), then move the base centre from the origin to (50, 25, 5) mm
    mv = model.FeatureManager.InsertMoveCopyBody2(50 * MM, 25 * MM, 5 * MM, 0.0, 0.0, 0.0, 0.0,
                                                  0.0, math.pi / 2, 0.0, False, 1)
    print(f"B  move/copy: feature {'yes' if mv else 'NO'}; {boss_centre(before, mass(model))} "
          f"(both bodies counted)")
    bs = model.GetBodies2(0, True)
    for label, args in (("main + tool arrays", (15903, bs[0], [bs[1]])),
                        ("selected bodies   ", (15903, NOTHING, NOTHING))):
        if bodies(model) < 2:
            break
        try:
            if args[1] is NOTHING:
                model.ClearSelection2(True)
                for bb in bs:
                    bb.Select2(True, NOTHING)
            cf = model.FeatureManager.InsertCombineFeature(*args)
            print(f"B  combine ({label}): feature {'yes' if cf else 'NO'}, bodies {bodies(model)}, "
                  f"{boss_centre(before, mass(model))}")
        except Exception as e:
            print(f"B  combine ({label}) failed: {type(e).__name__}: {e}")
except Exception as e:
    print(f"B failed: {type(e).__name__}: {e}")
print("done")
