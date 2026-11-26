"""Probe on the part left open by knuckle1_frame_d3 after 'Bar 4: no loft gives the expected volume'.
Run with that part still open in SolidWorks, from the folder of the generated script:
    python sw_probe_bar4.py knuckle1_frame_d3_solidworks > probe_bar4.log 2>&1
and send probe_bar4.log. It reads the three section sketches of Bar 4 (the last three sketches in the tree) and
the bar's spline in 'Frame skeleton', then tries lofts and deletes each one:
  A  start + end sections, no centerline          (2 sections, as Bars 1-3 were built)
  B  start + middle, C  middle + end              (which pair fails)
  D  all three, no centerline                     (what the script tried)
  E  all three, skeleton spline selected Mark 4   (what the script tried first)
  F  all three, spline converted into its own 3D sketch (SketchUseEdge3), that segment Mark 4
  G  all three, that sketch feature Mark 4
  H  all three, skeleton spline as a guide curve (Mark 2)
  I  start + end, skeleton spline Mark 4          (centerline of a sketch with branches)
  J  start + end, converted spline Mark 4
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


def features():
    out, f = [], call0(model, "FirstFeature")
    while f is not None:
        out.append(f)
        f = call0(f, "GetNextFeature")
    return out


# ---- the sections and the spline ----
fs = features()
print(f"features: {len(fs)}; last 12: {[(f.Name, call0(f, 'GetTypeName2')) for f in fs[-12:]]}")
sk = [f for f in fs if call0(f, "GetTypeName2") == "ProfileFeature"][-3:]
print(f"Bar {BAR} sections (last three sketches): {[f.Name for f in sk]}")
for f in sk:
    s = call0(f, "GetSpecificFeature2")
    segs = list(call0(s, "GetSketchSegments") or [])
    pts = [(round(p.X / MM, 3), round(p.Y / MM, 3), round(p.Z / MM, 3)) for p in (call0(s, "GetSketchPoints2") or [])]
    print(f"  {f.Name}: status {call0(s, 'GetConstrainedStatus')}, segments {len(segs)}, points (sketch mm) {pts}")
print(f"expected section radii {[r for _, r, _, _ in b['sec']]}, joints {[S.JOINTS[k] for k in b['joints']]}")

skel = model.FeatureByName("Frame skeleton")
sks = call0(skel, "GetSpecificFeature2")
A, C = [tuple(c * MM for c in S.JOINTS[k]) for k in (b["joints"][0], b["joints"][-1])]
same = lambda u, v: max(abs(u[i] - v[i]) for i in range(3)) < 1e-7
spl = None
allsegs = list(call0(sks, "GetSketchSegments") or [])
print(f"skeleton segments: {len(allsegs)}")
for g in allsegs:
    e = S.seg_ends(g)
    if e and ((same(e[0], A) and same(e[1], C)) or (same(e[0], C) and same(e[1], A))):
        spl = g
print(f"Bar {BAR} spline found by its end points: {spl is not None}")
if spl is None and len(allsegs) >= BAR:
    spl = allsegs[BAR - 1]
    print(f"  taken by order: segment {BAR}, ends {S.seg_ends(spl)}")


# ---- converted spline in its own 3D sketch ----
def convert():
    model.ClearSelection2(True)
    ok = spl.Select4(False, NOTHING)
    sm.Insert3DSketch(False)
    if call0(sm, "ActiveSketch") is None:
        sm.Insert3DSketch(False)
    model.ClearSelection2(True)
    ok = ok and spl.Select4(False, NOTHING)
    r = sm.SketchUseEdge3(False, False)
    act = call0(sm, "ActiveSketch")
    n = len(list(call0(act, "GetSketchSegments") or [])) if act is not None else None
    sm.Insert3DSketch(True)
    f = model.FeatureByPositionReverse(0)
    f.Name = f"Bar {BAR} path"
    return f, dict(selected=ok, SketchUseEdge3=r, segments=n)


got = step("convert the spline into its own 3D sketch (SketchUseEdge3)", convert)
path_f = got[0] if got else None
path_seg = None
if path_f is not None:
    segs = list(call0(call0(path_f, "GetSpecificFeature2"), "GetSketchSegments") or [])
    path_seg = segs[0] if segs else None
    print(f"  path sketch segments: {len(segs)}")


# ---- lofts ----
def loft(label, secs, extra=None):
    def run():
        model.ClearSelection2(True)
        ok = [sk[k].Select2(i > 0, 1) for i, k in enumerate(secs)]
        if extra is not None:
            obj, mark, how = extra
            ok.append(obj.Select4(True, sd(mark)) if how == "Select4" else obj.Select2(True, mark))
        n = model.SelectionManager.GetSelectedObjectCount2(-1)
        v0 = mass()
        f = fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                      0.0, 0.0, 0, True, True, True, 0)
        if f is None:
            return dict(selected=ok, count=n, loft=None)
        v1 = mass()
        model.ClearSelection2(True); f.Select2(False, 0); ext.DeleteSelection2(0)
        return dict(selected=ok, count=n, volume_added=round(v1 - v0, 3))
    step(label, run)


loft("A start + end, no centerline", [0, 2])
loft("B start + middle, no centerline", [0, 1])
loft("C middle + end, no centerline", [1, 2])
loft("D all three, no centerline", [0, 1, 2])
if spl is not None:
    loft("E all three, skeleton spline Mark 4", [0, 1, 2], (spl, 4, "Select4"))
if path_seg is not None:
    loft("F all three, converted spline Mark 4", [0, 1, 2], (path_seg, 4, "Select4"))
if path_f is not None:
    loft("G all three, path sketch feature Mark 4", [0, 1, 2], (path_f, 4, "Select2"))
if spl is not None:
    loft("H all three, skeleton spline guide curve Mark 2", [0, 1, 2], (spl, 2, "Select4"))
    loft("I start + end, skeleton spline Mark 4", [0, 2], (spl, 4, "Select4"))
if path_seg is not None:
    loft("J start + end, converted spline Mark 4", [0, 2], (path_seg, 4, "Select4"))
print("done")
