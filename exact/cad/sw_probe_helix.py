"""Probe: helical wires (stent) in SolidWorks. Run on the Windows machine with SolidWorks:
    python sw_probe_helix.py > probe_helix.log
and send the log. Each step prints what SolidWorks returned; nothing is assumed.

Stent numbers (mm): helix radius 1.5 about the z axis (Front Plane normal), pitch 16, wire diameter 0.12.
  A. one helix (InsertHelix, height and pitch, 16 x 16) from a circle on the Front Plane, variants of
     Clockwised and StartAngle; swept with a circular profile (CreateDefinition(swFmSweep),
     CircularProfile). Expected volume pi 0.06^2 x 18.5695 = 0.21002 mm^3. The planar end faces (centre,
     normal) show where the helix starts and which way it turns;
  B. two helices of opposite hand in one part, both swept: union 0.4158 +- 0.003 mm^3 (sum 0.42003);
  C. circular pattern of the two sweeps, 10 instances about the z axis (axis from the Top and Right
     Planes): union of 20 wires 16 high, 3.939 +- 0.005 mm^3;
  D. the same with helices 16.12 high and a cut above z = 16: 3.938 +- 0.005 mm^3 (the stent's flat top).
"""
import math
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
R, WIRE_D, PITCH = 1.5, 0.12, 16.0


def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


def sw_constant(app, enum, member):
    import os
    exe = app.GetExecutablePath
    exe = str(exe() if callable(exe) else exe)
    for f in (exe, os.path.dirname(exe)):
        for sub in ((), ("api",)):
            path = os.path.join(f, *sub, "swconst.tlb")
            if not os.path.exists(path):
                continue
            tl = pythoncom.LoadTypeLib(path)
            for i in range(tl.GetTypeInfoCount()):
                if tl.GetDocumentation(i)[0] != enum:
                    continue
                ti = tl.GetTypeInfo(i)
                for j in range(ti.GetTypeAttr().cVars):
                    vd = ti.GetVarDesc(j)
                    if ti.GetNames(vd.memid)[0] == member:
                        return int(vd.value)
    raise LookupError(f"{enum}.{member} not found")


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True
try:
    app.SetUserPreferenceToggle(10, False)
except Exception:
    pass
try:
    FM_SWEEP = sw_constant(app, "swFeatureNameID_e", "swFmSweep")
    print(f"swFmSweep = {FM_SWEEP}")
except Exception as ex:
    FM_SWEEP = None
    print(f"swFmSweep not read: {type(ex).__name__}: {ex}")


def new_part():
    m = call0(app, "NewPart")
    try:                                         # MMGS, so dimensions read in mm
        m.Extension.SetUserPreferenceInteger(sw_constant(app, "swUserPreferenceIntegerValue_e", "swUnitSystem"),
                                             0, sw_constant(app, "swUnitSystem_e", "swUnitSystem_MMGS"))
    except Exception:
        pass
    return m


def mass(model):
    model.ForceRebuild3(False)
    mp = call0(model.Extension, "CreateMassProperty")
    if mp is None:
        return None
    return round(mp.Volume / MM ** 3, 5), tuple(round(c / MM, 4) for c in mp.CenterOfMass)


def bodies(model):
    try:
        b = model.GetBodies2(0, True)
        return list(b) if b else []
    except Exception as ex:
        print(f"  GetBodies2: {type(ex).__name__}: {ex}")
        return []


def planar_faces(model):
    """Centre (from the face box) and normal of every planar face: the wire ends."""
    out = []
    for b in bodies(model):
        for f in (call0(b, "GetFaces") or []):
            try:
                s = call0(f, "GetSurface")
                if not call0(s, "IsPlane"):
                    continue
                bx = call0(f, "GetBox")
                c = [(bx[k] + bx[k + 3]) / 2 / MM for k in range(3)]
                n = call0(f, "Normal")
                out.append(((round(c[0], 4), round(c[1], 4), round(c[2], 4)),
                            tuple(round(v, 4) for v in n)))
            except Exception as ex:
                out.append(f"face: {type(ex).__name__}: {ex}")
    return out


def base_circle(model, plane="Front Plane"):
    model.ClearSelection2(True)
    ok = model.Extension.SelectByID2(plane, "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
    sm = model.SketchManager
    sm.InsertSketch(True)
    sm.AddToDB = True
    c = sm.CreateCircleByRadius(0.0, 0.0, 0.0, R * MM)
    sm.AddToDB = False
    sm.InsertSketch(True)
    f = model.FeatureByPositionReverse(0)
    print(f"  base circle on {plane} (selected {ok}): circle {c is not None}, feature {f.Name}")
    return f


def helix(model, sk, clockwised, start, height=16.0):
    model.ClearSelection2(True)
    sk.Select2(False, 0)
    before = model.FeatureByPositionReverse(0).Name
    try:
        r = model.InsertHelix(False, clockwised, False, False, 2, height * MM, PITCH * MM, height / PITCH,
                              0.0, start)
    except Exception as ex:
        print(f"  InsertHelix: FAILED {type(ex).__name__}: {ex}")
        return None
    f = model.FeatureByPositionReverse(0)
    made = f.Name != before
    print(f"  InsertHelix(clockwised={clockwised}, start={start:.4f}) returned {r}; new feature "
          f"{f.Name if made else None} ({call0(f, 'GetTypeName2') if made else ''})")
    return f if made else None


def sweep(model, path):
    fm = model.FeatureManager
    model.ClearSelection2(True)
    ok = path.Select2(True, 4)
    if FM_SWEEP is None:
        print("  no sweep: swFmSweep unknown")
        return None
    try:
        d = fm.CreateDefinition(FM_SWEEP)
        d.CircularProfile = True
        d.CircularProfileDiameter = WIRE_D * MM
        for name, val in (("TwistControlType", 0), ("PathAlignmentType", 0), ("Merge", True),
                          ("MergeSmoothFaces", True), ("AutoSelect", True)):
            try:
                setattr(d, name, val)
            except Exception as ex:
                print(f"    {name}: {type(ex).__name__}")
        f = fm.CreateFeature(d)
        print(f"  sweep (path selected {ok}): {'feature ' + f.Name if f is not None else 'None'}; "
              f"mass {mass(model)}; bodies {len(bodies(model))}")
        return f
    except Exception as ex:
        print(f"  sweep: FAILED {type(ex).__name__}: {ex}")
        return None


print("== A. one helix and one sweep; expected volume 0.21002 mm^3")
for cw, st in ((False, 0.0), (True, 0.0), (False, math.pi / 2)):
    print(f"-- clockwised {cw}, start angle {st:.4f} rad")
    try:
        m = new_part()
        sk = base_circle(m)
        h = helix(m, sk, cw, st)
        if h is not None:
            sweep(m, h)
            print(f"  planar faces (centre mm, normal): {planar_faces(m)}")
    except Exception:
        traceback.print_exc()

print("== B. two helices of opposite hand, both swept; expected union 0.4158 +- 0.003 mm^3 (sum 0.42003)")
mB, sweepsB = None, []
try:
    mB = new_part()
    for cw in (False, True):
        sk = base_circle(mB)
        h = helix(mB, sk, cw, 0.0)
        if h is not None:
            s = sweep(mB, h)
            if s is not None:
                sweepsB.append(s)
except Exception:
    traceback.print_exc()

print("== C. circular pattern of the two sweeps, 10 about z; expected union 3.939 +- 0.005 mm^3")
try:
    if mB is not None and len(sweepsB) == 2:
        m = mB
        m.ClearSelection2(True)
        ok1 = m.Extension.SelectByID2("Top Plane", "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
        ok2 = m.Extension.SelectByID2("Right Plane", "PLANE", 0, 0, 0, True, 0, NOTHING, 0)
        a = m.InsertAxis2(True)
        ax = m.FeatureByPositionReverse(0)
        print(f"  axis from Top and Right Planes (selected {ok1}, {ok2}): returned {a}; feature {ax.Name} "
              f"({call0(ax, 'GetTypeName2')})")
        fm = m.FeatureManager
        for label, mk in (
                ("FeatureCircularPattern4 (7 args)",
                 lambda: fm.FeatureCircularPattern4(10, 2 * math.pi, False, "NULL", False, True, False)),
                ("FeatureCircularPattern5 (14 args)",
                 lambda: fm.FeatureCircularPattern5(10, 2 * math.pi, False, "NULL", False, True, False, False,
                                                    False, False, 1, 0.0, "NULL", False))):
            for axis_mark, feat_mark in ((1, 4), (2, 1)):
                m.ClearSelection2(True)
                okA = ax.Select2(False, axis_mark)
                okF = all(s.Select2(True, feat_mark) for s in sweepsB)
                try:
                    f = mk()
                    print(f"  {label}, axis mark {axis_mark}, features mark {feat_mark} (selected {okA}, {okF}):"
                          f" {'feature ' + f.Name if f is not None else 'None'}; mass {mass(m)}; "
                          f"bodies {len(bodies(m))}")
                    if f is not None:
                        raise StopIteration
                except StopIteration:
                    raise
                except Exception as ex:
                    print(f"  {label}, marks {axis_mark}/{feat_mark}: FAILED {type(ex).__name__}: {ex}")
except StopIteration:
    pass
except Exception:
    traceback.print_exc()

print("== D. helices 16.12 high, pattern, cut above z = 16; expected 3.938 +- 0.005 mm^3")
try:
    m = new_part()
    sw = []
    for cw in (False, True):
        sk = base_circle(m)
        h = helix(m, sk, cw, 0.0, height=16.12)
        if h is not None:
            s = sweep(m, h)
            if s is not None:
                sw.append(s)
    m.ClearSelection2(True)
    m.Extension.SelectByID2("Top Plane", "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
    m.Extension.SelectByID2("Right Plane", "PLANE", 0, 0, 0, True, 0, NOTHING, 0)
    m.InsertAxis2(True)
    ax = m.FeatureByPositionReverse(0)
    fm = m.FeatureManager
    done = False
    for axis_mark, feat_mark in ((1, 4), (2, 1)):
        m.ClearSelection2(True)
        ax.Select2(False, axis_mark)
        for s in sw:
            s.Select2(True, feat_mark)
        try:
            f = fm.FeatureCircularPattern4(10, 2 * math.pi, False, "NULL", False, True, False)
        except Exception as ex:
            f = None
            print(f"  pattern: FAILED {type(ex).__name__}: {ex}")
        print(f"  pattern (marks {axis_mark}/{feat_mark}): {f is not None}; mass {mass(m)}")
        if f is not None:
            done = True
            break
    before = mass(m)
    print(f"  before the cut: {before}")
    m.ClearSelection2(True)
    m.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
    sm = m.SketchManager
    sm.InsertSketch(True)
    sm.AddToDB = True
    sm.CreateCircleByRadius(0.0, 0.0, 0.0, 3 * R * MM)
    sm.AddToDB = False
    sm.InsertSketch(True)
    sk = m.FeatureByPositionReverse(0)
    for rev, flip in ((False, False), (True, False), (False, True), (True, True)):
        m.ClearSelection2(True)
        sk.Select2(False, 0)
        f = fm.FeatureCut4(True, False, rev, 0, 0, 1.0 * MM, 0.0, False, False, False, False,
                           0.0, 0.0, False, False, False, False, False, True, True,
                           False, False, False, 3, 16.0 * MM, flip, False)
        print(f"  cut z 16..17 (reverse={rev}, flip_start={flip}): {f is not None}; mass {mass(m)}")
        if f is not None:
            v = mass(m)
            if v is not None and before is not None and 0 < before[0] - v[0] < 0.1:
                break
            m.ClearSelection2(True)
            f.Select2(False, 0)
            m.Extension.DeleteSelection2(0)
    print(f"  planar faces: {len(planar_faces(m))}")
except Exception:
    traceback.print_exc()
print("done")
