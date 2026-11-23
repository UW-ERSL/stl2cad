"""Thirteenth probe: patterns that follow their seed. Run:  python sw_probe13.py

Test part: a 60 x 10 mm bar (x 0..60, y 0..10) with 9 trapezoidal teeth on top (root 4, top 2,
height 2.5, pitch 6, first tooth at x 3..7), extruded 10 mm. After each pattern the first
tooth's top width is changed from 2 to 1.2 through its dimension, and the volume shows whether
all teeth followed.

A. Linear Pattern FEATURE of one tooth boss along the bar's bottom edge (sw_probe8 never reached
   the pattern call: Feature.Select2 was given a SelectData, not an integer mark).
   A1 FeatureLinearPattern5, A2 FeatureLinearPattern4, A3 FeatureLinearPattern2.
B. Sketch linear pattern with its spacing and count dimensions on (probe 8 had them off).
C. One Equal relation on 5 lines selected together (fewer relations than pairwise Equal).
"""
import math
import time

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001


def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


def settle(seconds=0.3):
    end = time.time() + seconds
    while time.time() < end:
        pythoncom.PumpWaitingMessages()
        time.sleep(0.02)


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True


def new_part():
    model = call0(app, "NewPart")
    if model is None:
        model = app.NewDocument(app.GetUserPreferenceStringValue(8), 0, 0, 0)
    return model


def open_sketch(model):
    model.ClearSelection2(True)
    model.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
    if model.SelectionManager.GetSelectedObjectCount2(-1) == 0:
        f = call0(model, "FirstFeature")
        while f is not None:
            if call0(f, "GetTypeName2") == "RefPlane":
                f.Select2(False, 0)
                break
            f = call0(f, "GetNextFeature")
    sm = model.SketchManager
    sm.InsertSketch(True)
    settle()
    sm.AddToDB = True
    try:
        sm.InferenceMode = False
    except Exception:
        pass
    sm.DisplayWhenAdded = False
    return sm.CreateCircleByRadius(-40 * MM, -40 * MM, 0.0, 0.5 * MM)


def active_sketch(model):
    sk = call0(model, "GetActiveSketch2")
    return sk if sk is not None else model.SketchManager.ActiveSketch


def close_sketch(model, helper):
    sm = model.SketchManager
    sm.AddToDB = False
    sm.DisplayWhenAdded = True
    if helper is not None:
        helper.ConstructionGeometry = True
    model.ClearSelection2(True)
    sm.InsertSketch(True)
    f = model.FeatureByPositionReverse(0)
    sk = call0(f, "GetSpecificFeature2")
    try:
        st = call0(sk, "GetConstrainedStatus")
        st = {1: "unknown", 2: "under defined", 3: "fully defined", 4: "over defined"}.get(st, str(st))
    except Exception as e:
        st = f"status not available ({type(e).__name__})"
    rels = {}
    try:
        rm = call0(sk, "RelationManager")
        for r in (rm.GetRelations(0) or []):
            t = call0(r, "GetRelationType")
            rels[t] = rels.get(t, 0) + 1
    except Exception as e:
        rels = f"relations not read ({type(e).__name__}: {e})"
    return st, rels


def extrude(model, depth=0.010):
    f = model.FeatureByPositionReverse(0)
    model.ClearSelection2(True)
    f.Select2(False, 0)
    return model.FeatureManager.FeatureExtrusion3(
        True, False, False, 0, 0, depth, 0.0, False, False, False, False, 0.0, 0.0,
        False, False, False, False, True, True, True, 0, 0.0, False)


def volume(model):
    model.ForceRebuild3(False)
    mp = call0(model.Extension, "CreateMassProperty")
    return 0.0 if mp is None else mp.Volume / MM ** 3


def line(sm, a, b):
    return sm.CreateLine(a[0] * MM, a[1] * MM, 0.0, b[0] * MM, b[1] * MM, 0.0)


def arc(sm, c, a, b):
    return sm.CreateArc(c[0] * MM, c[1] * MM, 0.0, a[0] * MM, a[1] * MM, 0.0, b[0] * MM, b[1] * MM, 0.0, 1)


def pick(model, ents):
    model.ClearSelection2(True)
    ok = True
    for k, e in enumerate(ents):
        if e == "origin":
            ok = ok and model.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0, k > 0, 0, NOTHING, 0)
        else:
            ok = ok and e.Select4(k > 0, NOTHING)
    return ok


def rel(model, name, ents):
    if pick(model, ents):
        model.SketchAddConstraints(name)
    model.ClearSelection2(True)


def dim(model, ents, x, y):
    d = model.AddDimension2(x * MM, y * MM, 0.0) if pick(model, ents) else None
    model.ClearSelection2(True)
    return d is not None


ends = lambda g: (call0(g, "GetStartPoint2"), call0(g, "GetEndPoint2"))


def delete_helper(model, h):
    model.ClearSelection2(True)
    if h is not None and h.Select4(False, NOTHING):
        model.EditDelete()
    model.ClearSelection2(True)


def npoints(model):
    sk = active_sketch(model)
    return len(call0(sk, "GetSketchPoints2") or [])



def dimvalue(model, d, mm):
    """Set a display dimension's value (mm) and rebuild."""
    d.GetDimension2(0).SystemValue = mm * MM
    call0(model, "EditRebuild3")


H, R0, TOP, PITCH, N = 2.5, 4.0, 2.0, 6.0, 9
S = PITCH * MM


def tooth(x0, top=TOP):
    side = (R0 - top) / 2
    return [(x0, 10), (x0 + side, 10 + H), (x0 + side + top, 10 + H), (x0 + R0, 10)]


def tooth_area(top):
    return (R0 + top) / 2 * H


BAR = 60 * 10


def feature_case(label, make):
    try:
        model = new_part()
        h = open_sketch(model)
        sm = model.SketchManager
        for a, b in zip([(0, 0), (60, 0), (60, 10), (0, 10)], [(60, 0), (60, 10), (0, 10), (0, 0)]):
            line(sm, a, b)
        close_sketch(model, h)
        extrude(model)
        h = open_sketch(model)
        pts = tooth(3)
        tl = [line(sm, a, b) for a, b in zip(pts, pts[1:] + pts[:1])]
        delete_helper(model, h)
        model.ClearSelection2(True)
        tl[1].Select4(False, NOTHING)
        dtop = model.AddDimension2(5 * MM, 14 * MM, 0.0)
        close_sketch(model, None)
        boss = extrude(model)
        model.ClearSelection2(True)
        # direction: the bar's bottom front edge along x (mark 1); the tooth feature (mark 4)
        e_ok = model.Extension.SelectByRay(30 * MM, -1 * MM, -1 * MM, 0.0, 0.7071068, 0.7071068, 0.0005, 1, False, 1, 0)
        how = "ray"
        if not e_ok:
            e_ok = model.Extension.SelectByID2("", "EDGE", 30 * MM, 0.0, 0.0, False, 1, NOTHING, 0)
            how = "id"
        t_ok = boss.Select2(True, 4)
        f = make(model)
        v2 = volume(model)
        exp2 = (BAR + N * tooth_area(TOP)) * 10
        print(f"{label}: edge {e_ok} ({how}), tooth {t_ok}; pattern {'yes' if f else 'NO'}; "
              f"volume {v2:.3f} (expected {exp2:.3f})")
        if f and dtop is not None:
            dimvalue(model, dtop, 1.2)
            v3 = volume(model)
            print(f"{label}: top width -> 1.2: volume {v3:.3f} (all teeth follow: "
                  f"{(BAR + N * tooth_area(1.2)) * 10:.3f}; only the first: "
                  f"{exp2 - (tooth_area(TOP) - tooth_area(1.2)) * 10:.3f})")
    except Exception as ex:
        print(f"{label}: failed: {type(ex).__name__}: {ex}")


print("A. Linear Pattern feature of one tooth boss, 9 at 6 mm; then top width 2 -> 1.2")
feature_case("A1 FeatureLinearPattern5", lambda m: m.FeatureManager.FeatureLinearPattern5(
    N, S, 1, 0.0, False, False, "", "", False, False, False, False, False, False, False, False,
    False, False, 0.0, 0.0, False, False))
feature_case("A2 FeatureLinearPattern4", lambda m: m.FeatureManager.FeatureLinearPattern4(
    N, S, 1, 0.0, False, False, "", "", False, False, False, False, False, False, False, False,
    False, False, 0.0, 0.0, False, False))
feature_case("A3 FeatureLinearPattern2", lambda m: m.FeatureManager.FeatureLinearPattern2(
    N, S, 1, 0.0, False, False, "", "", False))

# ---- B: sketch pattern with its dimensions ----------------------------------------------------
print("B. sketch linear pattern with spacing/count dimensions on; seed top width 2 -> 1.2")
try:
    model = new_part()
    h = open_sketch(model)
    sm = model.SketchManager
    t0 = tooth(3)
    xe = 3 + N * PITCH                                    # where the last copy's root line ends
    for a, b in zip([(xe, 10), (60, 10), (60, 0), (0, 0), (0, 10)], [(60, 10), (60, 0), (0, 0), (0, 10), (3, 10)]):
        line(sm, a, b)
    seed = [line(sm, a, b) for a, b in zip(t0, t0[1:])] + [line(sm, t0[3], (3 + PITCH, 10))]
    delete_helper(model, h)
    model.ClearSelection2(True)
    seed[1].Select4(False, NOTHING)
    dtop = model.AddDimension2(5 * MM, 15 * MM, 0.0)
    model.ClearSelection2(True)
    for x in seed:
        x.Select4(True, NOTHING)
    ok = sm.CreateLinearSketchStepAndRepeat(N, 1, S, 0.0, 0.0, math.pi / 2, "", True, False, False, True, False)
    print(f"  pattern returned {ok}; points {npoints(model)}")
    close_sketch(model, None)
    f = extrude(model)
    v2 = volume(model)
    print(f"  extrusion {'yes' if f else 'NO'}; volume {v2:.3f} (expected {(BAR + N * tooth_area(TOP)) * 10:.3f})")
    if f and dtop is not None:
        dimvalue(model, dtop, 1.2)
        print(f"  top width -> 1.2: volume {volume(model):.3f} (all teeth follow: "
              f"{v2 - N * (tooth_area(TOP) - tooth_area(1.2)) * 10:.3f}; only the first: "
              f"{v2 - (tooth_area(TOP) - tooth_area(1.2)) * 10:.3f})")
except Exception as ex:
    print(f"B failed: {type(ex).__name__}: {ex}")

# ---- C: one Equal relation on many entities ---------------------------------------------------
print("C. Equal on 5 lines selected together")
try:
    model = new_part()
    h = open_sketch(model)
    sm = model.SketchManager
    ls = [line(sm, (0, 3 * k), (10, 3 * k)) for k in range(5)]
    delete_helper(model, h)
    rel(model, "sgSAMELENGTH", ls)
    st, rels = close_sketch(model, None)
    print(f"  relations by type: {rels} (type 14 = Equal: 1 if one relation holds all 5 lines)")
except Exception as ex:
    print(f"C failed: {type(ex).__name__}: {ex}")
print("done")
