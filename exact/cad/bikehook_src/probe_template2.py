"""Probe 2: the bike hook's rib fillet and trunk lofts (after probe_hook.log). Run on the Windows machine:
    python sw_probe_hook2.py > probe_hook2.log
and send the log. Each step prints what SolidWorks returned and compares the volume and centre of mass with
CadQuery's values (EXPECT); nothing is assumed.

  Probe 1 (probe_hook.log): sweep, rib extrusion and end fillets agree with CadQuery; the rib bottom fillet
  (6 edges, r 1.5) was refused; 10 of 11 trunk lofts were refused and the one made was 5 % off (SolidWorks chose
  its own corner correspondence).
  C. Rib bottom fillet variants: each side's 3 edges, then edge by edge; r 1.5, then 1.0.
  E. Trunk lofts 1 (z -22/-28) and 4 (z -40/-46): no guide curve, KeepTangency off; one guide line (3D sketch)
     joining corner 0 of the two sections; a guide line at every corner (one 3D sketch each). Guide curves Mark 2.
     The best variant is then used for all 11 lofts.
  Sketches report their segment count (the helper circle must be gone: 24 per section).
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


print("== C. rib bottom fillet variants; expected", EXPECT["ribfil"])
for how in ("per side", "per edge"):
    for r in (RIB_F, 1.0):
        try:
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
            def bottom_edges():
                out = {1: [], -1: []}
                for b in bodies(m):
                    for e in call0(b, "GetEdges") or []:
                        fs = list(call0(e, "GetTwoAdjacentFaces2") or [])
                        if len(fs) != 2:
                            continue
                        info = [face_info(f) for f in fs]
                        side = [k for k, (pl, n, c) in enumerate(info)
                                if pl and abs(abs(n[1]) - 1) < 1e-6 and abs(abs(c[1]) - RIB_W) < 1e-3]
                        if len(side) != 1:
                            continue
                        pl, n, c = info[1 - side[0]]
                        if (not pl and c[2] < -85.0) or (pl and n[2] < -0.999):
                            out[1 if info[side[0]][2][1] > 0 else -1].append(e)
                return out
            E_ = bottom_edges()
            made = 0
            if how == "per side":
                for s in (1, -1):
                    f = fillet(m, bottom_edges()[s], r)
                    print(f"  side {s:+d} ({len(E_[s])} edges), r {r}: {'feature' if f is not None else 'None'}")
                    made += f is not None
            else:
                n = len(E_[1]) + len(E_[-1])
                for i in range(n):
                    es = bottom_edges(); flat = es[1] + es[-1]
                    if not flat:
                        break
                    f = fillet(m, [flat[0]], r)
                    made += f is not None
                    if f is None:
                        print(f"  edge {i}: None")
                        break
            print(f"  {how}, r {r}: {made} fillet feature(s)")
            report(m, f"{how}, r {r}", "ribfil" if r == RIB_F else "ribfil1")
        except Exception:
            traceback.print_exc()

print("== E. trunk loft variants")
best = None
for variant in ("no tangency", "one guide", "all guides"):
    res = []
    for k in (0, 3):
        try:
            res.append(pair_test(k, variant))
        except Exception:
            traceback.print_exc(); res.append(False)
    if all(res) and best is None:
        best = variant
print(f"  first variant right on both pairs: {best}")
if best is not None:
    try:
        m = new_part()
        sk = [sketch(m, plane_at_z(m, z), ents, f"section z = {z}") for z, ents in SECS]
        made = 0
        for k in range(len(sk) - 1):
            guides = []
            if best != "no tangency":
                c0, c1 = corners(SECS[k][1]), corners(SECS[k + 1][1])
                pairs = [(c0[0], c1[0])] if best == "one guide" else list(zip(c0, c1))
                guides = [line3d(m, a, b)[0] for a, b in pairs]
            f, ok = loft(m, sk[k], sk[k + 1], guides, keep_tangency=(best != "no tangency"))
            made += f is not None
            if f is None:
                print(f"  loft {k + 1}: None")
        print(f"  {made} of {len(sk) - 1} lofts")
        report(m, "11 lofts", "lofts")
    except Exception:
        traceback.print_exc()
print("done")
