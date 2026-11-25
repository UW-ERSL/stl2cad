"""Probe: how to fix the z coordinate of a 3D sketch point (swframe's joints).
Run on the Windows machine with SolidWorks:  python sw_probe_frame_z.py > probe_frame_z.log 2> probe_frame_z.err
and send both files.

sw_probe_frame (run 3) found that a dimension from a 3D sketch point to the Front Plane, or to any plane parallel
to it, is refused, while dimensions to the Right and Top Planes, point-to-point and line lengths work. Each test
below uses a fresh 3D sketch with point C = (12, 22, 32) mm, dimensioned to the Right and Top Planes (x 12, y 22),
then tries one way to fix z, and prints the result and the sketch's constrained status after closing
(3 = fully defined, 2 = under defined).
  Z1  Front Plane selected as a feature object (Select2), not by name;
  Z2  relation point on a plane (sgCOINCIDENT, then sgONPLANE... names) with a plane offset 32 mm from the Front Plane;
  Z3  a second point E = (12, 22, 0) on the Front Plane by relation, a line E-C, the line's length (32);
  Z4  a second point E = (12, 22, 37) and the point-to-point dimension C-E (a distance along z only: 5);
  Z5  a plane at 45 deg (through the x axis between the Top and Front Planes, made from a 3D sketch line), and
      the point-to-plane dimension to it (|y - z| / sqrt 2 = 7.071);
  Z6  dimension to the Front Plane after switching the 3D sketch plane (ISketchManager.Sketch3DSetPlane...,
      if any of the candidate calls exists).
"""
import math
import time
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
C = (12.0, 22.0, 32.0)


def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


def settle(t=0.3):
    pythoncom.PumpWaitingMessages(); time.sleep(t); pythoncom.PumpWaitingMessages()


def step(label, fn):
    try:
        r = fn()
        print(f"{label}: {r}")
        return r
    except Exception as ex:
        print(f"{label}: FAILED {type(ex).__name__}: {ex}")
        traceback.print_exc()
        return None


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True
try:
    app.SetUserPreferenceToggle(10, False)
except Exception:
    pass
model = call0(app, "NewPart")
sm, fm, ext = model.SketchManager, model.FeatureManager, model.Extension


def open3d():
    """Open a 3D sketch (the first call after NewPart opened nothing in sw_probe_frame)."""
    for flag in (False, True, False):
        model.ClearSelection2(True)
        sm.Insert3DSketch(flag)
        settle()
        sk = call0(sm, "ActiveSketch") or call0(model, "GetActiveSketch2")
        if sk is not None:
            return sk
    return None


def close3d():
    sm.Insert3DSketch(True)
    settle()
    f = model.FeatureByPositionReverse(0)
    try:
        st = call0(call0(f, "GetSpecificFeature2"), "GetConstrainedStatus")
    except Exception as ex:
        st = f"n/a ({type(ex).__name__})"
    return f, st


def points_at(sk, q):
    return [o for o in call0(sk, "GetSketchPoints2") if math.dist((o.X / MM, o.Y / MM, o.Z / MM), q) < 1e-6]


def dim(objs, plane=None, plane_obj=None, txt=None):
    model.ClearSelection2(True)
    ok = [o.Select4(k > 0, NOTHING) for k, o in enumerate(objs)]
    if plane:
        ok.append(ext.SelectByID2(plane, "PLANE", 0, 0, 0, True, 0, NOTHING, 0))
    if plane_obj is not None:
        ok.append(plane_obj.Select2(True, 0))
    t = txt or (C[0] + 3, C[1] + 3, C[2] + 3)
    d = model.AddDimension2(*[c * MM for c in t])
    return ok, None if d is None else round(d.GetDimension2(0).SystemValue / MM, 4)


def rel(objs, name, plane=None, plane_obj=None):
    model.ClearSelection2(True)
    ok = [o.Select4(k > 0, NOTHING) for k, o in enumerate(objs)]
    if plane:
        ok.append(ext.SelectByID2(plane, "PLANE", 0, 0, 0, True, 0, NOTHING, 0))
    if plane_obj is not None:
        ok.append(plane_obj.Select2(True, 0))
    r = model.SketchAddConstraints(name)
    return ok, r


def with_point_c(label, extra):
    """Fresh 3D sketch with point C dimensioned to the Right and Top Planes; extra(sk, c) fixes z."""
    print(f"== {label}")
    sk = open3d()
    if sk is None:
        print("  no 3D sketch opened"); return
    sm.AddToDB = True
    sm.CreatePoint(*[c * MM for c in C])
    sm.AddToDB = False
    c = points_at(sk, C)[0]
    step("  x to the Right Plane (12)", lambda: dim([c], "Right Plane"))
    c = points_at(sk, C)[0]
    step("  y to the Top Plane (22)", lambda: dim([c], "Top Plane"))
    try:
        extra(sk)
    except Exception as ex:
        print(f"  FAILED {type(ex).__name__}: {ex}")
        traceback.print_exc()
    f, st = close3d()
    print(f"  constrained status after closing: {st} (3 = fully defined)")


# ---- planes made before the sketches (a plane cannot be made while a sketch is open) ----
def offset_front(d, name):
    model.ClearSelection2(True)
    ext.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
    for flip in (0, 256):
        f = fm.InsertRefPlane(8 | flip, d * MM, 0, 0, 0, 0)
        if f is None:
            continue
        f = model.FeatureByPositionReverse(0)
        f.Name = name
        return f
    return None


level = step("plane 32 mm from the Front Plane (named 'Level 32')", lambda: offset_front(32.0, "Level 32"))

# Z1
with_point_c("Z1: Front Plane selected as a feature object", lambda sk: step(
    "  z to the Front Plane (object) (32)", lambda: dim([points_at(sk, C)[0]], plane_obj=model.FeatureByName("Front Plane"))))


# Z2
def z2(sk):
    for nm in ("sgCOINCIDENT", "sgONPLANE", "sgATINTERSECT", "sgPOINTONPLANE"):
        c = points_at(sk, C)[0]
        step(f"  relation {nm}: point C and the plane 'Level 32' (by name)", lambda: rel([c], nm, plane="Level 32"))
    c = points_at(sk, C)[0]
    step("  relation sgCOINCIDENT: point C and the plane 'Level 32' (object)", lambda: rel([c], "sgCOINCIDENT", plane_obj=level))
with_point_c("Z2: point on a plane by relation", z2)


# Z3
def z3(sk):
    E = (C[0], C[1], 0.0)
    sm.AddToDB = True
    ln = sm.CreateLine(*[v * MM for v in E], *[v * MM for v in C])
    sm.AddToDB = False
    step("  line E-C made", lambda: ln is not None)
    e = points_at(sk, E)
    cs = points_at(sk, C)
    step("  points at E, at C", lambda: (len(e), len(cs)))
    if len(cs) > 1:
        step("  line end at C coincident with C", lambda: rel([cs[0], cs[1]], "sgCOINCIDENT"))
    e = points_at(sk, E)
    step("  E on the Front Plane (sgCOINCIDENT, by name)", lambda: rel([e[0]], "sgCOINCIDENT", plane="Front Plane"))
    step("  line along z (sgALONGZ3D)", lambda: rel([ln], "sgALONGZ3D"))
    step("  line length (32)", lambda: dim([ln], txt=(C[0] + 3, C[1], 16.0)))
with_point_c("Z3: line from a point on the Front Plane, along z, its length", z3)


# Z4
def z4(sk):
    E = (C[0], C[1], C[2] + 5)
    sm.AddToDB = True
    sm.CreatePoint(*[v * MM for v in E])
    sm.AddToDB = False
    e, c = points_at(sk, E), points_at(sk, C)
    step("  point-to-point C-E, along z only (5)", lambda: dim([c[0], e[0]], txt=(C[0] + 3, C[1] + 3, C[2] + 2.5)))
with_point_c("Z4: a distance along z only", z4)


# Z5: a 45 deg plane through the x axis, made from a 3D sketch line before the test sketch
def plane45():
    sk = open3d()
    sm.AddToDB = True
    a = sm.CreateLine(0, 0, 0, 0.1, 0, 0)                    # along x through the origin
    b = sm.CreateLine(0, 0, 0, 0, 0.05, 0.05)                # in the 45 deg plane
    sm.AddToDB = False
    f, st = close3d()
    model.ClearSelection2(True)
    sd0 = call0(model.SelectionManager, "CreateSelectData"); sd0.Mark = 0
    sd1 = call0(model.SelectionManager, "CreateSelectData"); sd1.Mark = 1
    ok = a.Select4(False, sd0) and b.Select4(True, sd1)
    p = fm.InsertRefPlane(4, 0, 4, 0, 0, 0) if ok else None   # coincident with both lines
    if p is None:
        return ok, None
    p = model.FeatureByPositionReverse(0)
    p.Name = "Plane 45"
    return ok, p.Name
r45 = step("plane at 45 deg through two 3D sketch lines", plane45)
if r45 and r45[1]:
    with_point_c("Z5: dimension to the 45 deg plane", lambda sk: step(
        "  C to 'Plane 45' (|y - z| / sqrt 2 = 7.0711)", lambda: dim([points_at(sk, C)[0]], "Plane 45")))


# Z6
def z6(sk):
    for nm in ("Set3DSketchPlane", "Sketch3DSetActivePlane", "SetActiveSketchPlane"):
        step(f"  SketchManager.{nm} exists", lambda: hasattr(sm, nm))
with_point_c("Z6: calls that switch the 3D sketch plane", z6)
print("done")
