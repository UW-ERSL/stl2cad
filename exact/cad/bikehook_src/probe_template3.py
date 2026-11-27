"""Probe 3: the bike hook's rib fillet and trunk lofts (after probe_hook2.log). Run on the Windows machine:
    python sw_probe_hook3.py > probe_hook3.log
and send the log. Each step prints what SolidWorks returned and compares the volume and centre of mass with
CadQuery's values (EXPECT); nothing is assumed.

  Probe 2 (probe_hook2.log): fillets by edges found from their faces picked a concave edge (the volume grew);
  lofts with a guide line at every corner match CadQuery (loft 4: +0.030 %) but loft 1 was refused.
  C. Rib bottom fillet r 1.5 (then 1.0) on edges picked by a point on each (SelectByID2 "EDGE"): the flat bottom
     line, a point on each bottom arc, both sides (6 edges), all in one fillet, then one side at a time.
  E. Each of the 11 trunk lofts in its own part, guide lines at all 12 corners; if refused, one guide, then none.
     Volume and centre of mass against CadQuery's ruled loft of the same pair.
"""
import math
import sys
import traceback

sys.stdout.reconfigure(line_buffering=True)   # the log fills as the probe runs

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
    try:
        segs = len(call0(call0(f, "GetSpecificFeature2"), "GetSketchSegments") or [])
    except Exception as ex:
        segs = f"? ({type(ex).__name__})"
    what += f", {segs} segments"
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




def line3d(m, p0, p1):
    sm = m.SketchManager
    m.ClearSelection2(True)
    for attempt in range(3):
        sm.Insert3DSketch(True)
        if call0(sm, "ActiveSketch") is not None or call0(m, "GetActiveSketch2") is not None:
            break
    sm.AddToDB = True
    ln = None
    for attempt in range(4):
        ln = sm.CreateLine(*(v * MM for v in p0), *(v * MM for v in p1))
        if ln is not None:
            break
    sm.AddToDB = False
    sm.Insert3DSketch(True)
    return m.FeatureByPositionReverse(0), ln


def corners(ents):
    return [e[1] for e in ents if e[0] == "A"]


def loft(m, s0, s1, guides=(), keep_tangency=True):
    m.ClearSelection2(True)
    ok = s0.Select2(False, 1) and s1.Select2(True, 1)
    for g in guides:
        sd = call0(m.SelectionManager, "CreateSelectData"); sd.Mark = 2
        ok &= g.Select4(True, sd)
    f = m.FeatureManager.InsertProtrusionBlend2(False, keep_tangency, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                                0.0, 0.0, 0, True, True, True, 0)
    return f, ok


def expect_pair(k):
    return EXPECT["pairs"][k]


def pair_test(k, variant):
    m = new_part()
    (z0, e0), (z1, e1) = SECS[k], SECS[k + 1]
    s0 = sketch(m, plane_at_z(m, z0), e0, f"section z = {z0}")
    s1 = sketch(m, plane_at_z(m, z1), e1, f"section z = {z1}")
    c0, c1 = corners(e0), corners(e1)
    guides = []
    if variant == "one guide":
        guides = [line3d(m, c0[0], c1[0])[0]]
    elif variant == "all guides":
        guides = [line3d(m, a, b)[0] for a, b in zip(c0, c1)]
    f, ok = loft(m, s0, s1, guides, keep_tangency=(variant != "no tangency"))
    print(f"  loft {k + 1} ({variant}; selected {ok}): {'feature ' + f.Name if f is not None else 'None'}")
    good = False
    if f is not None:
        ev, ec = EXPECT["pairs"][k]
        got = mass(m)
        if got:
            v, c = got
            good = abs(v - ev) <= 1e-3 * ev
            print(f"    volume {v:.3f} (expected {ev:.3f}, {100 * (v / ev - 1):+.3f} %), centre {c} "
                  f"({math.dist(c, ec):.4f} mm off) -> {'OK' if good else 'DIFFERS'}")
    return good



RIB_PTS = @RIB_PTS@          # a point on each rib bottom edge (model mm)


def build_cradle_rib():
    m = new_part()
    path = sketch(m, "Top Plane", PATH, "path")
    prof = sketch(m, plane_at_z(m, Z_LEG), PROFILE, "profile")
    m.ClearSelection2(True); prof.Select2(False, 1); path.Select2(True, 4)
    m.FeatureManager.CreateFeature(m.FeatureManager.CreateDefinition(CONST["swFmSweep"]))
    rib = sketch(m, "Top Plane", RIB, "rib")
    m.ClearSelection2(True); rib.Select2(False, 0)
    m.FeatureManager.FeatureExtrusion3(True, False, False, CONST["swEndCondMidPlane"], 0, 2 * RIB_W * MM, 0.0,
                                       False, False, False, False, 0.0, 0.0, False, False, False, False, True,
                                       True, True, 0, 0.0, False)
    report(m, "cradle + rib", "rib")
    return m


def pick_edges(m, pts):
    m.ClearSelection2(True)
    got = []
    for p in pts:
        ok = m.Extension.SelectByID2("", "EDGE", p[0] * MM, p[1] * MM, p[2] * MM, False, 0, NOTHING, 0)
        o = m.SelectionManager.GetSelectedObject6(1, -1) if ok else None
        got.append(o)
        m.ClearSelection2(True)
    return got


print("== C. rib bottom fillet, edges picked by points; expected r 1.5", EXPECT["ribfil"], "r 1.0", EXPECT["ribfil1"])
for r in (RIB_F, 1.0):
    for how in ("both sides", "one side at a time"):
        try:
            m = build_cradle_rib()
            groups = [RIB_PTS] if how == "both sides" else [[p for p in RIB_PTS if p[1] > 0], [p for p in RIB_PTS if p[1] < 0]]
            made = 0
            for g in groups:
                es = pick_edges(m, g)
                n = sum(e is not None for e in es)
                f = fillet(m, [e for e in es if e is not None], r) if n else None
                print(f"  {how}: {n} of {len(g)} edges picked; fillet r {r}: {'feature' if f is not None else 'None'}")
                made += f is not None
            report(m, f"{how}, r {r}", "ribfil" if r == RIB_F else "ribfil1")
        except Exception:
            traceback.print_exc()

print("== E. each trunk loft alone, guides at all corners (else one, else none)")
for k in range(len(SECS) - 1):
    for variant in ("all guides", "one guide", "no tangency"):
        try:
            if pair_test(k, variant):
                break
        except Exception:
            traceback.print_exc()
print("done")
