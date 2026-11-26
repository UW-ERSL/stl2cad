"""Second probe on the part left by the failed knuckle1_frame_d3 run (run after sw_probe_bar4.py; same part open).
sw_probe_bar4: every loft of Bar 4's sections failed, even start + end without a centerline. Suspected cause: the
merge, not the loft. Bar 4 starts at joint 1, where Bar 1 ends with a section of the same radius (1.69 mm) and the
same centre on another plane, so the two end faces' circles cross; Bars 1-3 met only sections of other radii.
  M1  start + end, no centerline, Merge off      (a loft body forms: the merge is the problem)
  M2  all three, no centerline, Merge off
  M3  all three, skeleton spline Mark 4, Merge off
  X1  start moved 0.5 r back along the bar (plane offset from the start plane) + end, Merge on
  X2  that start + middle + end, Merge on
  X3  that start + middle + end, skeleton spline Mark 4, Merge on
Run:  python sw_probe_bar4b.py knuckle1_frame_d3_solidworks > probe_bar4b.log 2>&1
"""
import importlib
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
BAR = 4
b = S.BARS[BAR - 1]


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
    mp = call0(ext, "CreateMassProperty")
    return round(mp.Volume / MM ** 3, 3)


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


fs = features()
pl = {f.Name: f for f in fs if f.Name.startswith(f"Bar {BAR} plane")}
names = [f"Bar {BAR} plane start", f"Bar {BAR} plane section 1", f"Bar {BAR} plane end"]
sk = []
for nm in names:                      # the sketch right after each section plane
    i = [f.Name for f in fs].index(nm)
    sk.append(next(f for f in fs[i + 1:] if call0(f, "GetTypeName2") == "ProfileFeature"))
print(f"sections: {[(n, f.Name) for n, f in zip(names, sk)]}")
skel = call0(model.FeatureByName("Frame skeleton"), "GetSpecificFeature2")
A, C = [tuple(c * MM for c in S.JOINTS[k]) for k in (b["joints"][0], b["joints"][-1])]
same = lambda u, v: max(abs(u[i] - v[i]) for i in range(3)) < 1e-7
spl = next((g for g in (call0(skel, "GetSketchSegments") or [])
            if (lambda e: e and ((same(e[0], A) and same(e[1], C)) or (same(e[0], C) and same(e[1], A))))(S.seg_ends(g))),
           None)
print(f"spline found: {spl is not None}")


def loft(label, secs, centerline=False, merge=True):
    def run():
        model.ClearSelection2(True)
        ok = [s.Select2(i > 0, 1) for i, s in enumerate(secs)]
        if centerline:
            ok.append(spl.Select4(True, sd(4)))
        v0 = mass()
        nb0 = len(model.GetBodies2(0, True) or [])
        f = fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                      0.0, 0.0, 0, merge, True, True, 0)
        if f is None:
            return dict(selected=ok, loft=None)
        v1 = mass()
        nb1 = len(model.GetBodies2(0, True) or [])
        delete(f)
        return dict(selected=ok, volume_added=round(v1 - v0, 3), bodies=(nb0, nb1))
    step(label, run)


loft("M1 start + end, no centerline, Merge off", [sk[0], sk[2]], merge=False)
loft("M2 all three, no centerline, Merge off", sk, merge=False)
if spl is not None:
    loft("M3 all three, spline Mark 4, Merge off", sk, centerline=True, merge=False)


# ---- start section moved back along the bar by 0.5 r ----
def moved_start():
    i0, r0, t0, _ = b["sec"][0]
    p0 = S.JOINTS[b["joints"][i0]]
    d = 0.5 * r0
    q = tuple(p0[k] - d * t0[k] for k in range(3))
    for flip in (0, 256):
        model.ClearSelection2(True)
        pl[names[0]].Select2(False, 0)
        p = fm.InsertRefPlane(8 + flip, d * MM, 0, 0, 0, 0)
        if p is None:
            continue
        p = model.FeatureByPositionReverse(0)
        model.ClearSelection2(True)
        p.Select2(False, 0)
        sm.InsertSketch(True)
        if call0(sm, "ActiveSketch") is None:
            model.ClearSelection2(True); p.Select2(False, 0); sm.InsertSketch(True)
        act = call0(sm, "ActiveSketch")
        T = list(act.ModelToSketchTransform.ArrayData)
        Q = S.to_sk(T, q)
        if abs(Q[2]) < 1e-6:
            sm.AddToDB = True
            sm.CreateCircleByRadius(Q[0], Q[1], 0.0, r0 * MM)
            sm.AddToDB = False
            sm.InsertSketch(True)
            f = model.FeatureByPositionReverse(0)
            p.Name, f.Name = f"Bar {BAR} probe plane", f"Bar {BAR} probe start"
            return f, dict(flip=flip, offset_mm=d)
        sm.InsertSketch(True)
        delete(model.FeatureByPositionReverse(0))
        delete(p)
    return None


got = step("start section 0.5 r back along the bar (plane offset from the start plane)", moved_start)
if got:
    s0 = got[0]
    loft("X1 moved start + end, Merge on", [s0, sk[2]])
    loft("X2 moved start + middle + end, Merge on", [s0, sk[1], sk[2]])
    if spl is not None:
        loft("X3 moved start + middle + end, spline Mark 4, Merge on", [s0, sk[1], sk[2]], centerline=True)
print("done")
