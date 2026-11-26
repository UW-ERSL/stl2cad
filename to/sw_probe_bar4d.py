"""Fourth probe on the part left by the failed knuckle1_frame_d3 run (same part open, after sw_probe_bar4c).
sw_probe_bar4c: Bar 4's sections are right (centres and radii as expected), circles on Bar 3's planes loft, new
circles on Bar 4's planes do not. Bar 4 bends: its end planes are 96 deg apart (Bar 3's are parallel). Suspected
cause: the loft joins the circles' start points, which lie on each sketch's x axis; these axes point in unrelated
directions on Bar 4's planes, so the loft twists and intersects itself. A profile selected at a point puts the
loft's connector there. Each circle is selected at centre + r u, where u is the same direction (the normal to both
end planes, projected into each plane), so the connectors line up.
  P1  new circles at Bar 4's start and end (from sw_probe_bar4c), aligned, Merge off
  P2  the three new circles, aligned, Merge off
  P3  as P2, the skeleton spline Mark 4
  P4  the script's own three sections, aligned, spline Mark 4, Merge on (volume added)
  P5  as P4, no centerline
Each with the selection type SKETCHSEGMENT, then EXTSKETCHSEGMENT.
Run:  python sw_probe_bar4d.py knuckle1_frame_d3_solidworks > probe_bar4d.log 2>&1
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
fm, ext = model.FeatureManager, model.Extension
b4 = S.BARS[3]


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


unit = lambda v: (lambda n: tuple(x / n for x in v))(math.sqrt(sum(x * x for x in v)))
cross = lambda a, c: (a[1] * c[2] - a[2] * c[1], a[2] * c[0] - a[0] * c[2], a[0] * c[1] - a[1] * c[0])
dot = lambda a, c: sum(x * y for x, y in zip(a, c))
T = [unit(t) for _, _, t, _ in b4["sec"]]
P = [S.JOINTS[b4["joints"][i]] for i, _, _, _ in b4["sec"]]
Rr = [r for _, r, _, _ in b4["sec"]]
ref = unit(cross(T[0], T[-1]))
U = [unit(tuple(ref[k] - dot(ref, t) * t[k] for k in range(3))) for t in T]
print(f"end planes {math.degrees(math.acos(max(-1, min(1, dot(T[0], T[-1]))))):.1f} deg apart; u = {[tuple(round(x, 3) for x in u) for u in U]}")
at = [tuple(P[k][j] + Rr[k] * U[k][j] for j in range(3)) for k in range(3)]

fs = features()
by = {f.Name: f for f in fs}
new = [by.get(f"probe c Bar 4 {k}") for k in range(3)]
names = [f.Name for f in fs]
own = []
for w in ("start", "section 1", "end"):
    i = names.index(f"Bar 4 plane {w}")
    own.append(next(f for f in fs[i + 1:] if call0(f, "GetTypeName2") == "ProfileFeature"))
print(f"new circles: {[f.Name if f else None for f in new]}; own sections: {[f.Name for f in own]}")
skel = call0(model.FeatureByName("Frame skeleton"), "GetSpecificFeature2")
A, C = [tuple(c * MM for c in S.JOINTS[k]) for k in (b4["joints"][0], b4["joints"][-1])]
same = lambda u, v: max(abs(u[i] - v[i]) for i in range(3)) < 1e-7
spl = next((g for g in (call0(skel, "GetSketchSegments") or [])
            if (lambda e: e and ((same(e[0], A) and same(e[1], C)) or (same(e[0], C) and same(e[1], A))))(S.seg_ends(g))),
           None)


def loft(label, sk, idx, centerline=False, merge=False):
    for typ in ("SKETCHSEGMENT", "EXTSKETCHSEGMENT"):
        def run():
            model.ClearSelection2(True)
            ok = [ext.SelectByID2("", typ, *[c * MM for c in at[k]], n > 0, 1, NOTHING, 0) for n, k in enumerate(idx)]
            if centerline:
                ok.append(spl.Select4(True, sd(4)))
            smg = model.SelectionManager
            cnt = smg.GetSelectedObjectCount2(-1)
            typs = [smg.GetSelectedObjectType3(i + 1, -1) for i in range(cnt)]
            v0 = mass()
            f = fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                          0.0, 0.0, 0, merge, True, True, 0)
            if f is None:
                return dict(selected=ok, types=typs, loft=None)
            v1 = mass()
            delete(f)
            return dict(selected=ok, types=typs, loft="made", volume_added=round(v1 - v0, 3))
        r = step(f"{label} [{typ}]", run)
        if r and r.get("loft"):
            return


# the profile sketches selected at the aligned points must be the ones meant: hide the others by order of tests
if all(new):
    loft("P1 new start + end, aligned, Merge off", new, [0, 2])
    loft("P2 new three, aligned, Merge off", new, [0, 1, 2])
    if spl is not None:
        loft("P3 new three, aligned, spline Mark 4, Merge off", new, [0, 1, 2], centerline=True)
    for f in new:                      # the new circles lie on the own sections: removed before P4, P5
        delete(f)
if spl is not None:
    loft("P4 own three, aligned, spline Mark 4, Merge on", own, [0, 1, 2], centerline=True, merge=True)
loft("P5 own three, aligned, no centerline, Merge on", own, [0, 1, 2], merge=True)
print(f"expected volume of Bar 4 (CadQuery, own): "
      f"{next((s['own'] for s in S.STEPS if s.get('bar') == 3), None)}")
print("done")
