"""Probe: the SolidWorks calls the bike hook's feature tree needs (bikehook_exact.py). Run on the Windows machine:
    python sw_probe_hook.py > probe_hook.log
and send the log. Each step prints what SolidWorks returned and compares the volume and centre of mass with
CadQuery's values (EXPECT); nothing is assumed.

  A. Cradle: J path (lines and R15 arcs) on the Top Plane (model xz); a plane parallel to the Front Plane at the
     leg's top end (z = Z_LEG); the tapered section (lines and tangent arcs) on it; Sweep (CreateDefinition
     swFmSweep, profile Mark 1, path Mark 4; else InsertProtrusionSwept4);
  B. Rib: its outline (lines and arcs) on the Top Plane, Extruded Boss mid-plane over 2 x RIB_W;
  C. Rib bottom fillet RIB_F on the edges between the rib's side faces and its bottom (CreateDefinition
     swFmFillet, edges found from the faces);
  D. End fillets END_R on every edge of the two planar end faces;
  E. Trunk: two of its sections (12 corners, tangent fillets) on planes parallel to the Front Plane; Loft
     (InsertProtrusionBlend2) with the sketches selected plainly and by a point at corner 0 (where SolidWorks
     starts the connector); then all 11 lofts between consecutive sections in one part.
"""
import math
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
Z_LEG, Z_LIP = @Z_LEG@, @Z_LIP@
RIB_W, RIB_Z, RIB_F, END_R = @RIB_W@, @RIB_Z@, @RIB_F@, @END_R@
EXPECT = @EXPECT@
# sketch entities in model mm: ('L', p0, p1) line, ('A', start, mid, end) arc
PATH = @PATH@
PROFILE = @PROFILE@
RIB = @RIB@
SECS = @SECS@


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
CONST = {}
for enum, member, default in (("swFeatureNameID_e", "swFmSweep", 17), ("swFeatureNameID_e", "swFmFillet", 1),
                              ("swSimpleFilletType_e", "swConstRadiusFillet", 0),
                              ("swEndConditions_e", "swEndCondMidPlane", 6)):
    try:
        CONST[member] = sw_constant(app, enum, member)
    except Exception as ex:
        CONST[member] = default
        print(f"{member} not read ({type(ex).__name__}); using {default}")
print("constants:", CONST)


def new_part():
    m = call0(app, "NewPart")
    try:
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
    return round(mp.Volume / MM ** 3, 3), tuple(round(c / MM, 4) for c in mp.CenterOfMass)


def report(model, what, key):
    got = mass(model)
    ev, ec = EXPECT[key]
    if got is None:
        print(f"  {what}: no mass properties; expected {ev:.3f} mm^3")
        return False
    v, c = got
    dc = math.dist(c, ec)
    ok = abs(v - ev) <= 1e-3 * ev and dc <= 0.05
    print(f"  {what}: volume {v:.3f} (expected {ev:.3f}, {100 * (v / ev - 1):+.3f} %), centre {c}, "
          f"{dc:.4f} mm from expected -> {'OK' if ok else 'DIFFERS'}")
    return ok


def to_sk(T, p):
    x = [c * MM for c in p]
    return tuple((sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3))


def open_sketch(model, plane):
    sm = model.SketchManager
    for attempt in range(3):
        model.ClearSelection2(True)
        if isinstance(plane, str):
            model.Extension.SelectByID2(plane, "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
        else:
            plane.Select2(False, 0)
        sm.InsertSketch(True)
        sk = call0(sm, "ActiveSketch")
        if sk is not None:
            break
    else:
        raise RuntimeError("no sketch opened")
    T = list(sk.ModelToSketchTransform.ArrayData)
    sm.AddToDB = True
    helper = sm.CreateCircleByRadius(-0.5, -0.5, 0.0, 0.5 * MM)
    return sm, T, helper


def close_sketch(model, sm, helper):
    sm.AddToDB = False
    try:
        model.ClearSelection2(True)
        helper.Select4(False, NOTHING)
        model.EditDelete()
    except Exception as ex:
        print(f"  helper not deleted: {type(ex).__name__}")
    sm.InsertSketch(True)
    return model.FeatureByPositionReverse(0)


def draw(sm, T, ents):
    n = 0
    for e in ents:
        if e[0] == "L":
            s = sm.CreateLine(*to_sk(T, e[1]), *to_sk(T, e[2]))
        else:
            s = sm.Create3PointArc(*to_sk(T, e[1]), *to_sk(T, e[3]), *to_sk(T, e[2]))
        n += s is not None
    return n


def sketch(model, plane, ents, what):
    sm, T, h = open_sketch(model, plane)
    n = draw(sm, T, ents)
    off = max(abs(to_sk(T, e[1])[2]) for e in ents) / MM
    f = close_sketch(model, sm, h)
    print(f"  {what}: {n} of {len(ents)} entities; points {off:.2e} mm off the sketch plane; feature {f.Name}")
    return f


def plane_at_z(model, z):
    """Plane parallel to the Front Plane (model xy) at model z (mm); side checked through a sketch transform."""
    for flip in (0, 256):
        model.ClearSelection2(True)
        model.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
        p = model.FeatureManager.InsertRefPlane(8 | flip, abs(z) * MM, 0, 0, 0, 0)
        if p is None:
            continue
        f = model.FeatureByPositionReverse(0)
        sm, T, h = open_sketch(model, f)
        off = to_sk(T, (0.0, 0.0, z))[2] / MM
        sm.AddToDB = False
        sm.InsertSketch(True)
        sk = model.FeatureByPositionReverse(0)
        model.ClearSelection2(True); sk.Select2(False, 0); model.Extension.DeleteSelection2(0)
        if abs(off) < 1e-6:
            print(f"  plane at z = {z} (flag {flip}): {f.Name}")
            return f
        model.ClearSelection2(True); f.Select2(False, 0); model.Extension.DeleteSelection2(0)
    raise RuntimeError(f"no plane at z = {z}")


def bodies(model):
    b = model.GetBodies2(0, True)
    return list(b) if b else []


def face_info(f):
    """(is plane, normal, box centre mm)."""
    bx = call0(f, "GetBox")
    c = tuple((bx[k] + bx[k + 3]) / 2 / MM for k in range(3))
    try:
        pl = bool(call0(call0(f, "GetSurface"), "IsPlane"))
    except Exception:
        pl = False
    n = tuple(call0(f, "Normal")) if pl else None
    return pl, n, c


def fillet(model, edges, r):
    fm = model.FeatureManager
    d = fm.CreateDefinition(CONST["swFmFillet"])
    d.Initialize(CONST["swConstRadiusFillet"])
    d.Edges = VARIANT(pythoncom.VT_ARRAY | pythoncom.VT_DISPATCH, edges)
    d.DefaultRadius = r * MM
    return fm.CreateFeature(d)


# ---------------------------------------------------------------------------------------------------------------
print("== A. cradle sweep; expected", EXPECT["sweep"])
m = None
try:
    m = new_part()
    path = sketch(m, "Top Plane", PATH, "path (Top Plane)")
    pl = plane_at_z(m, Z_LEG)
    prof = sketch(m, pl, PROFILE, "profile")
    fm = m.FeatureManager
    f = None
    try:
        m.ClearSelection2(True)
        ok = prof.Select2(False, 1) and path.Select2(True, 4)
        d = fm.CreateDefinition(CONST["swFmSweep"])
        try:
            d.CircularProfile = False
        except Exception:
            pass
        f = fm.CreateFeature(d)
        print(f"  CreateDefinition(swFmSweep), profile Mark 1, path Mark 4 (selected {ok}): "
              f"{'feature ' + f.Name if f is not None else 'None'}")
    except Exception as ex:
        print(f"  CreateDefinition sweep: FAILED {type(ex).__name__}: {ex}")
    if f is None:
        try:
            m.ClearSelection2(True)
            ok = prof.Select2(False, 1) and path.Select2(True, 4)
            f = fm.InsertProtrusionSwept4(False, False, 0, False, False, 0, 0, False, 0.0, 0.0, 0, 0, True, True,
                                          True, 0.0, True, False, 0.0, 0)
            print(f"  InsertProtrusionSwept4 (selected {ok}): {'feature ' + f.Name if f is not None else 'None'}")
        except Exception as ex:
            print(f"  InsertProtrusionSwept4: FAILED {type(ex).__name__}: {ex}")
    if f is not None:
        report(m, "sweep", "sweep")
except Exception:
    traceback.print_exc()

print("== B. rib, Extruded Boss mid-plane; expected", EXPECT["rib"])
try:
    rib = sketch(m, "Top Plane", RIB, "rib outline (Top Plane)")
    fm = m.FeatureManager
    f = None
    for rev in (False, True):
        m.ClearSelection2(True)
        rib.Select2(False, 0)
        f = fm.FeatureExtrusion3(True, False, rev, CONST["swEndCondMidPlane"], 0, 2 * RIB_W * MM, 0.0, False, False,
                                 False, False, 0.0, 0.0, False, False, False, False, True, True, True, 0, 0.0, False)
        print(f"  FeatureExtrusion3 mid-plane {2 * RIB_W} mm (reverse {rev}): {'feature ' + f.Name if f else 'None'}")
        if f is not None:
            break
    if f is not None:
        report(m, "cradle + rib", "rib")
except Exception:
    traceback.print_exc()

print("== C. rib bottom fillet; expected", EXPECT["ribfil"])
try:
    edges, seen = [], 0
    for b in bodies(m):
        for e in call0(b, "GetEdges") or []:
            fs = list(call0(e, "GetTwoAdjacentFaces2") or [])
            if len(fs) != 2:
                continue
            info = [face_info(f) for f in fs]
            side = [k for k, (pl, n, c) in enumerate(info) if pl and abs(abs(n[1]) - 1) < 1e-6 and abs(abs(c[1]) - RIB_W) < 1e-3]
            if len(side) != 1:
                continue
            seen += 1
            pl, n, c = info[1 - side[0]]
            # the other face: the rib's bottom (curved, or planar facing -z), not the cradle and not the rib's end
            if (not pl and c[2] < -85.0) or (pl and n[2] < -0.999):
                edges.append(e)
    print(f"  edges on the rib's side faces: {seen}; bottom edges chosen: {len(edges)} (expected 6)")
    f = fillet(m, edges, RIB_F)
    print(f"  fillet r {RIB_F}: {'feature ' + f.Name if f is not None else 'None'}")
    if f is not None:
        report(m, "rib bottom fillet", "ribfil")
except Exception:
    traceback.print_exc()

print("== D. end fillets; expected", EXPECT["ends"])
try:
    edges = []
    for b in bodies(m):
        for fc in call0(b, "GetFaces") or []:
            pl, n, c = face_info(fc)
            if pl and abs(abs(n[2]) - 1) < 1e-6 and min(abs(c[2] - Z_LEG), abs(c[2] - Z_LIP)) < 1e-3:
                edges += list(call0(fc, "GetEdges") or [])
    print(f"  edges of the two end faces: {len(edges)}")
    f = fillet(m, edges, END_R)
    print(f"  fillet r {END_R}: {'feature ' + f.Name if f is not None else 'None'}")
    if f is not None:
        report(m, "end fillets", "ends")
except Exception:
    traceback.print_exc()

print("== E. trunk lofts; one pair expected", EXPECT["pair"], "; all 11 expected", EXPECT["lofts"])
for how in ("plain", "point at corner 0"):
    try:
        m2 = new_part()
        sk = []
        for z, ents in SECS[:2]:
            sk.append(sketch(m2, plane_at_z(m2, z), ents, f"section z = {z}"))
        m2.ClearSelection2(True)
        if how == "plain":
            ok = sk[0].Select2(False, 1) and sk[1].Select2(True, 1)
        else:
            sd = call0(m2.SelectionManager, "CreateSelectData"); sd.Mark = 1
            ok = True
            for k, (z, ents) in enumerate(SECS[:2]):
                p = ents[0][1]
                ok &= m2.Extension.SelectByID2(sk[k].Name, "SKETCH", p[0] * MM, p[1] * MM, p[2] * MM, k > 0, 1,
                                               NOTHING, 0)
        f = m2.FeatureManager.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                                     0.0, 0.0, 0, True, True, True, 0)
        print(f"  loft, sketches selected {how} ({ok}): {'feature ' + f.Name if f is not None else 'None'}")
        if f is not None:
            report(m2, f"loft ({how})", "pair")
    except Exception:
        traceback.print_exc()
try:
    m3 = new_part()
    sk = [sketch(m3, plane_at_z(m3, z), ents, f"section z = {z}") for z, ents in SECS]
    made = 0
    for k in range(len(sk) - 1):
        m3.ClearSelection2(True)
        sk[k].Select2(False, 1); sk[k + 1].Select2(True, 1)
        f = m3.FeatureManager.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                                     0.0, 0.0, 0, True, True, True, 0)
        made += f is not None
        if f is None:
            print(f"  loft {k + 1}: None")
    print(f"  {made} of {len(sk) - 1} lofts made; sketches stay visible to the next loft? (each used by two lofts)")
    report(m3, "11 lofts", "lofts")
except Exception:
    traceback.print_exc()
print("done")
