"""Third probe on the part left by the failed knuckle1_frame_d3 run (same part open).
sw_probe_bar4b: Bar 4 does not loft even with Merge off, and not with a new start circle (no relations) and its
own end section. So the loft itself fails. This probe reads the sections' geometry and lofts new circles that
have no relations:
  geometry  each Bar 3 and Bar 4 section: circle centre (model mm, from the sketch transform), circumference/2pi,
            relations in the sketch
  T0  Bar 3's own two sections again, Merge off           (SolidWorks lofts at all in this state)
  T1  new circles on Bar 3's planes, Merge off            (new circles loft where the script's did)
  T2  new circles on Bar 4's start and end planes, Merge off
  T3  new circles on Bar 4's start, middle, end planes, Merge off
  T4  as T3, the skeleton spline Mark 4
  T5  new circles on Bar 4's start and middle planes, Merge off
Run:  python sw_probe_bar4c.py knuckle1_frame_d3_solidworks > probe_bar4c.log 2>&1
"""
import importlib
import math
import sys
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

S = importlib.import_module(sys.argv[1] if len(sys.argv) > 1 else "knuckle1_frame_d3_solidworks")
call0, MM, NOTHING = S.call0, 0.001, VARIANT(pythoncom.VT_DISPATCH, None)
app = win32com.client.Dispatch("SldWorks.Application")
model = app.ActiveDoc
fm, ext, sm = model.FeatureManager, model.Extension, model.SketchManager


def step(label, fn):
    try:
        r = fn()
        print(f"{label}: {r}")
        return r
    except Exception as ex:
        print(f"{label}: FAILED {type(ex).__name__}: {ex}")
        traceback.print_exc(file=sys.stdout)
        return None


def mass():
    model.ForceRebuild3(False)
    return round(call0(ext, "CreateMassProperty").Volume / MM ** 3, 3)


def sd(mark):
    s = call0(model.SelectionManager, "CreateSelectData")
    s.Mark = mark
    return s


def delete(f):
    model.ClearSelection2(True); f.Select2(False, 0); ext.DeleteSelection2(0)


def features():
    out, f = [], call0(model, "FirstFeature")
    while f is not None:
        out.append(f)
        f = call0(f, "GetNextFeature")
    return out


def to_model(T, q):
    """Sketch point (m) -> model point (mm), inverting ModelToSketchTransform (see to_sk)."""
    R = [[T[3 * i + j] for j in range(3)] for i in range(3)]
    y = [q[j] / T[12] - T[9 + j] for j in range(3)]
    return tuple(round(sum(y[j] * R[i][j] for j in range(3)) / MM, 4) for i in range(3))


fs = features()
names = [f.Name for f in fs]


def sections(bar, wheres):
    out = []
    for w in wheres:
        i = names.index(f"Bar {bar} plane {w}")
        out.append((fs[i], next(f for f in fs[i + 1:] if call0(f, "GetTypeName2") == "ProfileFeature")))
    return out


B3 = sections(3, ["start", "end"])
B4 = sections(4, ["start", "section 1", "end"])
for bar, secs, b in ((3, B3, S.BARS[2]), (4, B4, S.BARS[3])):
    for k, (p, f) in enumerate(secs):
        s = call0(f, "GetSpecificFeature2")
        T = list(call0(s, "ModelToSketchTransform").ArrayData)
        pts = [to_model(T, (o.X, o.Y, o.Z)) for o in (call0(s, "GetSketchPoints2") or [])]
        segs = list(call0(s, "GetSketchSegments") or [])
        L = [round(call0(g, "GetLength") / MM / (2 * math.pi), 4) for g in segs]
        try:
            nrel = call0(s, "RelationManager").GetRelationsCount(0)
        except Exception as ex:
            nrel = f"{type(ex).__name__}"
        i, r, t, _ = b["sec"][k]
        print(f"Bar {bar} {p.Name} / {f.Name}: centre {pts}, radius {L}, relations {nrel}; "
              f"expected centre {S.JOINTS[b['joints'][i]]}, radius {r}")
skel = call0(model.FeatureByName("Frame skeleton"), "GetSpecificFeature2")
b4 = S.BARS[3]
A, C = [tuple(c * MM for c in S.JOINTS[k]) for k in (b4["joints"][0], b4["joints"][-1])]
same = lambda u, v: max(abs(u[i] - v[i]) for i in range(3)) < 1e-7
spl = next((g for g in (call0(skel, "GetSketchSegments") or [])
            if (lambda e: e and ((same(e[0], A) and same(e[1], C)) or (same(e[0], C) and same(e[1], A))))(S.seg_ends(g))),
           None)


def new_circle(plane, p, r, name):
    for attempt in range(2):
        model.ClearSelection2(True)
        plane.Select2(False, 0)
        sm.InsertSketch(True)
        if call0(sm, "ActiveSketch") is not None:
            break
    act = call0(sm, "ActiveSketch")
    T = list(act.ModelToSketchTransform.ArrayData)
    Q = S.to_sk(T, p)
    sm.AddToDB = True
    sm.CreateCircleByRadius(Q[0], Q[1], 0.0, r * MM)
    sm.AddToDB = False
    sm.InsertSketch(True)
    f = model.FeatureByPositionReverse(0)
    f.Name = name
    print(f"  {name}: point at sketch z {Q[2] / MM:+.5f} mm")
    return f


def loft(label, secs, centerline=False):
    def run():
        model.ClearSelection2(True)
        ok = [s.Select2(i > 0, 1) for i, s in enumerate(secs)]
        if centerline:
            ok.append(spl.Select4(True, sd(4)))
        f = fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                      0.0, 0.0, 0, False, True, True, 0)
        if f is None:
            return dict(selected=ok, loft=None)
        n = len(model.GetBodies2(0, True) or [])
        delete(f)
        return dict(selected=ok, loft="made", bodies_with_it=n)
    step(label, run)


loft("T0 Bar 3's own sections, Merge off", [f for _, f in B3])
n3 = [new_circle(p, S.JOINTS[S.BARS[2]["joints"][i]], r, f"probe c Bar 3 {k}")
      for k, ((p, _), (i, r, _, _)) in enumerate(zip(B3, S.BARS[2]["sec"]))]
loft("T1 new circles on Bar 3's planes, Merge off", n3)
n4 = [new_circle(p, S.JOINTS[b4["joints"][i]], r, f"probe c Bar 4 {k}")
      for k, ((p, _), (i, r, _, _)) in enumerate(zip(B4, b4["sec"]))]
loft("T2 new circles on Bar 4's start and end planes, Merge off", [n4[0], n4[2]])
loft("T3 new circles on Bar 4's three planes, Merge off", n4)
if spl is not None:
    loft("T4 new circles on Bar 4's three planes, spline Mark 4, Merge off", n4, centerline=True)
loft("T5 new circles on Bar 4's start and middle planes, Merge off", [n4[0], n4[1]])
print("done")
