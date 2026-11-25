"""Probe: the SolidWorks calls swframe.py uses beyond swsweep's (3D frame of a TO part).
Run on the Windows machine with SolidWorks:  python sw_probe_frame.py  > probe_frame.log
and send the log. Each step prints what SolidWorks returned; nothing is assumed.

Geometry (mm): joints A = (10, 20, 30), B = (40, 25, 35); member A-B (length 30.82).
Steps:
  1. 3D sketch (Insert3DSketch): points A, B and the line A-B (AddToDB); GetSketchPoints2 (are X, Y, Z model
     coordinates?), the sketch's ModelToSketchTransform;
  2. in the 3D sketch: a point-to-plane dimension (point A and the Right, Top, Front Plane: expect 10, 20, 30), a
     point-to-point dimension (A, B: expect 30.82, the true length), coincident relations (line ends and points);
     the constrained status after closing;
  3. plane normal to the 3D line through A (InsertRefPlane: perpendicular to the segment Mark 0, coincident with
     the point Mark 1), checked through a sketch's transform; the same at B;
  4. circles on the two planes (r 3 at A, 2 at B), centre pierced by the 3D line (sgATPIERCE), a diameter
     dimension; a Loft between them without a centerline (expected frustum volume
     pi/3 x 30.82 x (9 + 6 + 4) = 613.2 mm^3), then with the 3D line as centerline (Mark 4);
  5. plane parallel to the Right Plane through the 3D sketch point B (InsertRefPlane: parallel Mark 0, coincident
     Mark 1); a semicircle (r 4) closed by a line on it, its centre made coincident with B (an external point),
     revolved about the line (expected sphere 268.08 mm^3, centred at B).
"""
import math
import time
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
PARALLEL, PERP, COINC = 1, 2, 4    # swRefPlaneReferenceConstraints_e
A, B = (10.0, 20.0, 30.0), (40.0, 25.0, 35.0)
L = math.dist(A, B)
T_AB = tuple((B[k] - A[k]) / L for k in range(3))


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
        traceback.print_exc()
        return None


def sd(model, mark):
    s = call0(model.SelectionManager, "CreateSelectData")
    s.Mark = mark
    return s


def to_sk(T, p):
    x = [c * MM for c in p]
    return tuple((sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3))


def plane_check(model, plane, p, t):
    """(sketch z of p, sketch z of p + t) on the plane, mm: (0, +-1) when the plane is normal to t through p."""
    model.ClearSelection2(True)
    plane.Select2(False, 0)
    sm = model.SketchManager
    sm.InsertSketch(True)
    sk = call0(sm, "ActiveSketch") or call0(model, "GetActiveSketch2")
    T = list(sk.ModelToSketchTransform.ArrayData)
    z0 = to_sk(T, p)[2] / MM
    z1 = to_sk(T, tuple(p[k] + t[k] for k in range(3)))[2] / MM
    sm.InsertSketch(True)
    f = model.FeatureByPositionReverse(0)
    model.ClearSelection2(True); f.Select2(False, 0); model.Extension.DeleteSelection2(0)
    return round(z0, 6), round(z1, 6)


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True
model = call0(app, "NewPart")
sm, fm, ext = model.SketchManager, model.FeatureManager, model.Extension
try:
    app.SetUserPreferenceToggle(10, False)        # swInputDimValOnCreate: no Modify dialog
except Exception:
    pass

# ---- 1. 3D sketch ----
print("== 1. 3D sketch")
step("active document is the new part", lambda: call0(call0(app, "ActiveDoc"), "GetTitle") == call0(model, "GetTitle"))


def active():
    pythoncom.PumpWaitingMessages()
    time.sleep(0.3)
    pythoncom.PumpWaitingMessages()
    a = call0(sm, "ActiveSketch")
    b = call0(model, "GetActiveSketch2")
    return a or b


def command_3d():
    """swCommands_e member for a 3D sketch, read from swconst.tlb."""
    import os
    exe = app.GetExecutablePath
    exe = str(exe() if callable(exe) else exe)
    for folder in (exe, os.path.dirname(exe)):
        for path in (os.path.join(folder, "swconst.tlb"), os.path.join(folder, "api", "swconst.tlb")):
            if not os.path.exists(path):
                continue
            tl = pythoncom.LoadTypeLib(path)
            for i in range(tl.GetTypeInfoCount()):
                if tl.GetDocumentation(i)[0] != "swCommands_e":
                    continue
                ti = tl.GetTypeInfo(i)
                for j in range(ti.GetTypeAttr().cVars):
                    vd = ti.GetVarDesc(j)
                    nm = ti.GetNames(vd.memid)[0]
                    if nm in ("swCommands_3DSketch", "swCommands_3D_Sketch", "swCommands_Insert3DSketch"):
                        return nm, int(vd.value)
                names = [ti.GetNames(ti.GetVarDesc(j).memid)[0] for j in range(ti.GetTypeAttr().cVars)]
                return "not found; members with 3D: " + ", ".join(n for n in names if "3D" in n.upper())[:600], None
    return "swconst.tlb not found", None


# ways to open a 3D sketch, tried in turn until a sketch is active (each printed)
sk = None
opened_by = None
for label, fn in (("SketchManager.Insert3DSketch(True)", lambda: sm.Insert3DSketch(True)),
                  ("SketchManager.Insert3DSketch(False)", lambda: sm.Insert3DSketch(False)),
                  ("ModelDoc2.Insert3DSketch2(True)", lambda: model.Insert3DSketch2(True)),
                  ("ModelDoc2.Insert3DSketch()", lambda: model.Insert3DSketch()),
                  ("RunCommand(swCommands 3D sketch)", lambda: (lambda c: (c, app.RunCommand(c[1], "") if c[1] else None))(command_3d()))):
    model.ClearSelection2(True)
    step(label, fn)
    sk = step(f"  active sketch after {label}", active)
    if sk is not None:
        opened_by = label
        step("  Is3D", lambda: call0(sk, "Is3D"))
        break
print(f"3D sketch opened by: {opened_by}")
if sk is None:
    print("no way opened a 3D sketch; the rest of the probe needs one")
    raise SystemExit
step("ModelToSketchTransform", lambda: [round(v, 6) for v in sk.ModelToSketchTransform.ArrayData])
sm.AddToDB = True
pA = step("CreatePoint A", lambda: sm.CreatePoint(*[c * MM for c in A]))
pB = step("CreatePoint B", lambda: sm.CreatePoint(*[c * MM for c in B]))
ln = step("CreateLine A-B", lambda: sm.CreateLine(*[c * MM for c in A], *[c * MM for c in B]))
if ln is None:
    ln = step("CreateLine A-B (second try, after the points)", lambda: sm.CreateLine(*[c * MM for c in A], *[c * MM for c in B]))
sm.AddToDB = False
pts = step("GetSketchPoints2 (mm)", lambda: [(round(p.X / MM, 4), round(p.Y / MM, 4), round(p.Z / MM, 4)) for p in call0(sk, "GetSketchPoints2")])
objs = list(call0(sk, "GetSketchPoints2")) if pts else []

# ---- 2. dimensions and relations in the 3D sketch ----
print("== 2. dimensions and relations")
atA = [o for o in objs if math.dist((o.X / MM, o.Y / MM, o.Z / MM), A) < 1e-6]
atB = [o for o in objs if math.dist((o.X / MM, o.Y / MM, o.Z / MM), B) < 1e-6]
print(f"points at A: {len(atA)}, at B: {len(atB)}")


def rel(objs_, name):
    model.ClearSelection2(True)
    ok = all(o.Select4(k > 0, NOTHING) for k, o in enumerate(objs_))
    model.SketchAddConstraints(name)
    return ok


for grp, nm in ((atA, "A"), (atB, "B")):
    for o in grp[1:]:
        step(f"coincident at {nm}", lambda: rel([grp[0], o], "sgCOINCIDENT"))
for plane, expect in (("Right Plane", A[0]), ("Top Plane", A[1]), ("Front Plane", A[2])):
    def dim_plane():
        model.ClearSelection2(True)
        ok = atA[0].Select4(False, NOTHING) and ext.SelectByID2(plane, "PLANE", 0, 0, 0, True, 0, NOTHING, 0)
        d = model.AddDimension2((A[0] + 3) * MM, (A[1] + 3) * MM, (A[2] + 3) * MM) if ok else None
        return (ok, None if d is None else round(d.GetDimension2(0).SystemValue / MM, 4), f"expected {expect}")
    step(f"dimension A to {plane}", dim_plane)


def dim_ab():
    model.ClearSelection2(True)
    ok = atA[0].Select4(False, NOTHING) and atB[0].Select4(True, NOTHING)
    d = model.AddDimension2(25 * MM, 25 * MM, 35 * MM) if ok else None
    return (ok, None if d is None else round(d.GetDimension2(0).SystemValue / MM, 4), f"expected {L:.4f} (true length)")
step("dimension A to B", dim_ab)
step("close 3D sketch (the same call that opened it)", lambda: {
    "SketchManager.Insert3DSketch(True)": lambda: sm.Insert3DSketch(True),
    "SketchManager.Insert3DSketch(False)": lambda: sm.Insert3DSketch(True),
    "ModelDoc2.Insert3DSketch2(True)": lambda: model.Insert3DSketch2(True),
    "ModelDoc2.Insert3DSketch()": lambda: model.Insert3DSketch()}.get(opened_by, lambda: sm.InsertSketch(True))())
step("sketch still active after closing (expect None)", lambda: call0(sm, "ActiveSketch"))
sk3 = model.FeatureByPositionReverse(0)
step("3D sketch constrained status (3 = fully defined; B is free here)",
     lambda: call0(call0(sk3, "GetSpecificFeature2"), "GetConstrainedStatus"))

# ---- 3. planes normal to the 3D line ----
print("== 3. planes normal to the 3D line")


def normal_plane(p_obj, p):
    model.ClearSelection2(True)
    ok = ln.Select4(False, sd(model, 0)) and p_obj.Select4(True, sd(model, 1))
    f = fm.InsertRefPlane(PERP, 0, COINC, 0, 0, 0) if ok else None
    if f is None:
        return ok, None
    f = model.FeatureByPositionReverse(0)
    return f, plane_check(model, f, p, T_AB)


rA = step("plane at A by objects (expect (0, +-1))", lambda: normal_plane(pA, A))
rB = step("plane at B by objects (expect (0, +-1))", lambda: normal_plane(pB, B))
plA = rA[0] if rA and not isinstance(rA[0], bool) else None
plB = rB[0] if rB and not isinstance(rB[0], bool) else None

# ---- 4. circles and lofts ----
print("== 4. circles and lofts")


def circle(plane, p, r):
    model.ClearSelection2(True)
    plane.Select2(False, 0)
    sm.InsertSketch(True)
    s = call0(sm, "ActiveSketch") or call0(model, "GetActiveSketch2")
    T = list(s.ModelToSketchTransform.ArrayData)
    C = to_sk(T, p)
    sm.AddToDB = True
    c = sm.CreateCircleByRadius(C[0], C[1], 0.0, r * MM)
    sm.AddToDB = False
    cen = [o for o in call0(s, "GetSketchPoints2") if abs(o.X - C[0]) < 1e-9 and abs(o.Y - C[1]) < 1e-9]
    model.ClearSelection2(True)
    ok1 = bool(cen) and cen[0].Select4(False, NOTHING) and ln.Select4(True, NOTHING)
    model.SketchAddConstraints("sgATPIERCE")
    model.ClearSelection2(True)
    c.Select4(False, NOTHING)
    d = model.AddDimension2(C[0] + r * MM, C[1] + r * MM, 0.0)
    sm.InsertSketch(True)
    f = model.FeatureByPositionReverse(0)
    st = call0(call0(f, "GetSpecificFeature2"), "GetConstrainedStatus")
    return f, dict(centre_z=round(C[2] / MM, 6), pierce_selected=ok1,
                   diameter=None if d is None else round(d.GetDimension2(0).SystemValue / MM, 4), status=st)


cA = step("circle r 3 at A", lambda: circle(plA, A, 3.0)) if plA else None
cB = step("circle r 2 at B", lambda: circle(plB, B, 2.0)) if plB else None
if cA and cB:
    for variant in ("no centerline", "centerline Mark 4"):
        def loft():
            model.ClearSelection2(True)
            cA[0].Select2(False, 1); cB[0].Select2(True, 1)
            if variant != "no centerline":
                ln.Select4(True, sd(model, 4))
            f = fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                          0.0, 0.0, 0, True, True, True, 0)
            if f is None:
                return None
            v = mass(model)
            model.ClearSelection2(True); f.Select2(False, 0); ext.DeleteSelection2(0)
            return v, "expected 613.2 mm^3"
        step(f"loft, {variant}", loft)

# ---- 5. node: plane parallel to the Right Plane through B, sphere ----
print("== 5. node")


def node_plane():
    model.ClearSelection2(True)
    ok = model.FeatureByName("Right Plane").Select2(False, 0) and pB.Select4(True, sd(model, 1))
    f = fm.InsertRefPlane(PARALLEL, 0, COINC, 0, 0, 0) if ok else None
    if f is None:
        return ok, None
    f = model.FeatureByPositionReverse(0)
    return f, plane_check(model, f, B, (1.0, 0.0, 0.0))


rN = step("plane parallel to the Right Plane through B (expect (0, +-1))", node_plane)
if rN and not isinstance(rN[0], bool):
    def sphere():
        r = 4.0
        model.ClearSelection2(True)
        rN[0].Select2(False, 0)
        sm.InsertSketch(True)
        s = call0(sm, "ActiveSketch") or call0(model, "GetActiveSketch2")
        T = list(s.ModelToSketchTransform.ArrayData)
        C = to_sk(T, B)
        top, bot, side = to_sk(T, (B[0], B[1] + r, B[2])), to_sk(T, (B[0], B[1] - r, B[2])), to_sk(T, (B[0], B[1], B[2] + r))
        turn = 1 if (side[0] - top[0]) * (bot[1] - top[1]) - (side[1] - top[1]) * (bot[0] - top[0]) > 0 else -1
        sm.AddToDB = True
        arc = sm.CreateArc(C[0], C[1], 0.0, top[0], top[1], 0.0, bot[0], bot[1], 0.0, turn)
        line = sm.CreateLine(bot[0], bot[1], 0.0, top[0], top[1], 0.0)
        sm.AddToDB = False
        cen = [o for o in call0(s, "GetSketchPoints2") if abs(o.X - C[0]) < 1e-9 and abs(o.Y - C[1]) < 1e-9]
        model.ClearSelection2(True)
        ok = bool(cen) and cen[0].Select4(False, NOTHING) and pB.Select4(True, NOTHING)
        model.SketchAddConstraints("sgCOINCIDENT")
        sm.InsertSketch(True)
        f = model.FeatureByPositionReverse(0)
        out = dict(arc=arc is not None, line=line is not None, centre_to_B_selected=ok, sketch_axes=[
            "xyz"[max(range(3), key=lambda i: abs(T[3 * i + j]))] for j in (0, 1)])
        for mark in (16, 4):
            model.ClearSelection2(True)
            f.Select2(False, 0)
            line.Select4(True, sd(model, mark))
            rv = fm.FeatureRevolve2(True, True, False, False, False, False, 0, 0, 2 * math.pi, 0.0,
                                    False, False, 0.0, 0.0, 0, 0.0, 0.0, True, True, True)
            out[f"revolve Mark {mark}"] = None if rv is None else (mass(model), "expected 268.083 mm^3 at B")
            if rv is not None:
                break
        return out
    step("sphere at B", sphere)
print("done")
