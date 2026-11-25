"""Probe: curved bars for the 3D frame (swframe): a 3D spline centre line, lofts through several sections along it,
an ellipse section, and a fillet where two bars meet.
Run on the Windows machine with SolidWorks:  python sw_probe_bars.py > probe_bars.log 2> probe_bars.err
and send both files. Each step prints what SolidWorks returned; nothing is assumed.

Geometry (mm): spline through P0 = (0, 0, 0), P1 = (15, 8, 4), P2 = (30, 10, 12), P3 = (45, 6, 20).
  1. 3D sketch: the spline through P0..P3 (CreateSpline2, then CreateSpline), its points, a dimension from P1 to the
     Right Plane (15) and to the Top Plane (8);
  2. planes normal to the spline at P0, P2 and P3 (InsertRefPlane perpendicular to the spline Mark 0, coincident with
     the spline's point Mark 1), checked through a sketch's transform;
  3. circles r 3 at P0, 2.2 at P2, 2 at P3 (centre pierced by the spline), then lofts with the spline as centerline:
     a. P0 and P3 sections only (CadQuery, sweeping along its own spline through the points: 1014.1 mm^3);
     b. P0, P2 and P3 sections (956.2 mm^3);
     c. P0, an ellipse at P2 (half-axes 2.8 and 1.6, major axis 30 deg from the sketch's x axis), P3 (920.6 mm^3);
     the spline SolidWorks makes may differ a little from CadQuery's, so a volume within a few % is the expected result;
  4. a second bar from P2 to Q = (30, 25, 12): a 3D sketch line, circles r 1.8 at P2 and r 1.5 at Q on planes normal to
     it, a straight loft; then the edges where the two bars meet (edges between a face of each loft) found from the
     body, and a fillet of 1 mm on them (FeatureFillet3, as swgen; marks 1 then 0); volume before and after.
"""
import math
import time
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
EMPTY = VARIANT(pythoncom.VT_EMPTY, None)
MM = 0.001
PERP, COINC = 2, 4
P = [(0.0, 0.0, 0.0), (15.0, 8.0, 4.0), (30.0, 10.0, 12.0), (45.0, 6.0, 20.0)]
Q = (30.0, 25.0, 12.0)


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


def mass():
    model.ForceRebuild3(False)
    mp = call0(model.Extension, "CreateMassProperty")
    if mp is None:
        return 0.0, (0, 0, 0)
    return round(mp.Volume / MM ** 3, 3), tuple(round(c / MM, 3) for c in mp.CenterOfMass)


def to_sk(T, p):
    x = [c * MM for c in p]
    return tuple((sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3))


def sd(mark):
    s = call0(model.SelectionManager, "CreateSelectData")
    s.Mark = mark
    return s


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True
try:
    app.SetUserPreferenceToggle(10, False)
except Exception:
    pass
model = call0(app, "NewPart")
sm, fm, ext = model.SketchManager, model.FeatureManager, model.Extension


def open3d():
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
    return model.FeatureByPositionReverse(0)


def points_at(sk, q):
    return [o for o in call0(sk, "GetSketchPoints2") if math.dist((o.X / MM, o.Y / MM, o.Z / MM), q) < 1e-6]


def sketch_on(plane, p, label):
    """Sketch on the plane; (sketch, T) when p is on it (sketch z 0), else None. The first sketch call may open
    nothing (seen after NewPart): repeated."""
    for attempt in range(3):
        model.ClearSelection2(True)
        plane.Select2(False, 0)
        sm.InsertSketch(True)
        settle()
        s = call0(sm, "ActiveSketch") or call0(model, "GetActiveSketch2")
        if s is not None:
            break
    if s is None:
        print(f"  {label}: no sketch opened"); return None
    T = list(s.ModelToSketchTransform.ArrayData)
    z = to_sk(T, p)[2] / MM
    if abs(z) > 1e-3:
        print(f"  {label}: sketch {z:+.4f} mm off the point")
    return s, T


def close_sketch():
    sm.InsertSketch(True)
    settle(0.1)
    return model.FeatureByPositionReverse(0)


# ---- 1. spline ----
print("== 1. 3D spline through P0..P3")
sk = step("open 3D sketch", open3d)
pts = []
for p in P:
    pts += [c * MM for c in p]
arr = VARIANT(pythoncom.VT_ARRAY | pythoncom.VT_R8, pts)
sm.AddToDB = True
spl = step("CreateSpline2(points, natural ends True)", lambda: sm.CreateSpline2(arr, True))
if spl is None:
    spl = step("CreateSpline(points)", lambda: sm.CreateSpline(arr))
sm.AddToDB = False
step("points of the sketch (mm)", lambda: [(round(o.X / MM, 3), round(o.Y / MM, 3), round(o.Z / MM, 3))
                                           for o in call0(sk, "GetSketchPoints2")])


def dim_plane(o, plane, txt):
    model.ClearSelection2(True)
    ok = o.Select4(False, NOTHING) and ext.SelectByID2(plane, "PLANE", 0, 0, 0, True, 0, NOTHING, 0)
    d = model.AddDimension2(*[c * MM for c in txt])
    return ok, None if d is None else round(d.GetDimension2(0).SystemValue / MM, 4)


p1 = points_at(sk, P[1])
step("points at P1", lambda: len(p1))
if p1:
    step("P1 to the Right Plane (15)", lambda: dim_plane(points_at(sk, P[1])[0], "Right Plane", (18, 11, 7)))
    step("P1 to the Top Plane (8)", lambda: dim_plane(points_at(sk, P[1])[0], "Top Plane", (18, 11, 7)))
skf = step("close the 3D sketch", close3d)
spline_pts = {}
try:
    s1 = call0(skf, "GetSpecificFeature2")
    for k, p in enumerate(P):
        o = [x for x in call0(s1, "GetSketchPoints2") if math.dist((x.X / MM, x.Y / MM, x.Z / MM), p) < 1e-6]
        spline_pts[k] = o[0] if o else None
    segs = list(call0(s1, "GetSketchSegments") or [])
    print(f"segments in the sketch: {len(segs)}; spline points found: {[k for k, o in spline_pts.items() if o]}")
    if spl is None and segs:
        spl = segs[0]
except Exception as ex:
    print(f"reading the sketch: FAILED {type(ex).__name__}: {ex}")

# ---- 2. planes normal to the spline ----
print("== 2. planes normal to the spline")


def tangent(k):
    a, b = (P[max(k - 1, 0)], P[min(k + 1, 3)])
    d = [b[i] - a[i] for i in range(3)]
    L = math.sqrt(sum(x * x for x in d))
    return tuple(x / L for x in d)


def normal_plane(k):
    model.ClearSelection2(True)
    ok = spl.Select4(False, sd(0)) and spline_pts[k].Select4(True, sd(1))
    f = fm.InsertRefPlane(PERP, 0, COINC, 0, 0, 0) if ok else None
    if f is None:
        return ok, None
    f = model.FeatureByPositionReverse(0)
    f.Name = f"Normal P{k}"
    got = sketch_on(f, P[k], f"plane at P{k}")
    if got is None:
        return f, "no sketch"
    s, T = got
    z0 = to_sk(T, P[k])[2] / MM
    t = tangent(k)
    z1 = to_sk(T, tuple(P[k][i] + t[i] for i in range(3)))[2] / MM
    sm.AddToDB = True
    sm.CreatePoint(-0.0005, -0.0005, 0.0)         # not empty: an empty sketch leaves no feature
    sm.AddToDB = False
    g = close_sketch()
    model.ClearSelection2(True); g.Select2(False, 0); ext.DeleteSelection2(0)
    # z1: the sketch z of P + the chord direction (about +-1 if the plane is normal to the spline there)
    return f, (round(z0, 4), round(z1, 3))


planes = {}
for k in (0, 2, 3):
    r = step(f"plane normal to the spline at P{k} (expect (0, about +-1))", lambda: normal_plane(k))
    if r and not isinstance(r[0], bool):
        planes[k] = r[0]

# ---- 3. sections and lofts ----
print("== 3. sections and lofts along the spline")


def section(k, r, ellipse=None):
    got = sketch_on(planes[k], P[k], f"section at P{k}")
    if got is None:
        return None
    s, T = got
    C = to_sk(T, P[k])
    sm.AddToDB = True
    if ellipse is None:
        e = sm.CreateCircleByRadius(C[0], C[1], 0.0, r * MM)
    else:
        a, b, ang = ellipse
        ca, sa = math.cos(math.radians(ang)), math.sin(math.radians(ang))
        e = sm.CreateEllipse(C[0], C[1], 0.0, C[0] + a * MM * ca, C[1] + a * MM * sa, 0.0,
                             C[0] - b * MM * sa, C[1] + b * MM * ca, 0.0)
    sm.AddToDB = False
    cen = [o for o in call0(s, "GetSketchPoints2") if abs(o.X - C[0]) < 1e-9 and abs(o.Y - C[1]) < 1e-9]
    model.ClearSelection2(True)
    ok = bool(cen) and cen[0].Select4(False, NOTHING) and spl.Select4(True, NOTHING)
    model.SketchAddConstraints("sgATPIERCE")
    f = close_sketch()
    st = call0(call0(f, "GetSpecificFeature2"), "GetConstrainedStatus")
    print(f"  section at P{k}: {'ellipse' if ellipse else 'circle'} made {e is not None}, centre pierced {ok}, status {st}")
    return f


def loft(secs, label):
    model.ClearSelection2(True)
    for i, s_ in enumerate(secs):
        s_.Select2(i > 0, 1)
    spl.Select4(True, sd(4))
    f = fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                  0.0, 0.0, 0, True, True, True, 0)
    if f is None:
        return None
    v = mass()
    f.Name = label
    return f, v


if all(k in planes for k in (0, 2, 3)):
    s0, s2, s3 = section(0, 3.0), section(2, 2.2), section(3, 2.0)
    for secs, label, expect in (([s0, s3], "Loft 2 sections", 1014.1), ([s0, s2, s3], "Loft 3 sections", 956.2)):
        r = step(f"{label} along the spline (expect about {expect} mm^3)", lambda: loft(secs, label))
        if r:
            model.ClearSelection2(True); r[0].Select2(False, 0); ext.DeleteSelection2(0)
    s2e = section(2, None, ellipse=(2.8, 1.6, 30.0))
    r = step("Loft with the ellipse at P2 (expect about 920.6 mm^3)", lambda: loft([s0, s2e, s3], "Loft ellipse"))
    if r:
        model.ClearSelection2(True); r[0].Select2(False, 0); ext.DeleteSelection2(0)
    bar1 = step("Loft 3 sections again (bar 1, kept)", lambda: loft([s0, s2, s3], "Bar 1"))
else:
    bar1 = None
    print("planes missing; sections and lofts skipped")

# ---- 4. second bar and the junction fillet ----
print("== 4. second bar from P2 to Q, fillet where the bars meet")


def bar2():
    s = open3d()
    sm.AddToDB = True
    a = sm.CreatePoint(*[c * MM for c in P[2]])
    b = sm.CreatePoint(*[c * MM for c in Q])
    ln = sm.CreateLine(*[c * MM for c in P[2]], *[c * MM for c in Q])
    sm.AddToDB = False
    close3d()
    t = (0.0, 1.0, 0.0)
    out = []
    for p, pt, r in ((P[2], a, 1.8), (Q, b, 1.5)):
        model.ClearSelection2(True)
        ok = ln.Select4(False, sd(0)) and pt.Select4(True, sd(1))
        f = fm.InsertRefPlane(PERP, 0, COINC, 0, 0, 0) if ok else None
        if f is None:
            return f"no plane at {p}"
        f = model.FeatureByPositionReverse(0)
        got = sketch_on(f, p, f"bar 2 section at {p}")
        s_, T = got
        C = to_sk(T, p)
        sm.AddToDB = True
        sm.CreateCircleByRadius(C[0], C[1], 0.0, r * MM)
        sm.AddToDB = False
        out.append(close_sketch())
    v0 = mass()
    model.ClearSelection2(True)
    out[0].Select2(False, 1); out[1].Select2(True, 1)
    f = fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                  0.0, 0.0, 0, True, True, True, 0)
    if f is None:
        return "no loft"
    f.Name = "Bar 2"
    return v0, mass()


if bar1:
    step("bar 2 (volume before, after)", bar2)

    def junction_edges():
        """Edges between a face of Bar 1 and a face of Bar 2 (late binding: methods called through call0)."""
        bodies = model.GetBodies2(0, True)
        body = bodies[0]
        out, seen = [], {}
        for e in call0(body, "GetEdges") or []:
            try:
                fs = call0(e, "GetTwoAdjacentFaces2")
                names = {call0(call0(f, "GetFeature"), "Name") for f in fs}
            except Exception as ex:
                seen[type(ex).__name__] = seen.get(type(ex).__name__, 0) + 1
                continue
            key = tuple(sorted(names))
            seen[key] = seen.get(key, 0) + 1
            if names == {"Bar 1", "Bar 2"}:
                out.append(e)
        print(f"  edges by the features of their two faces: {seen}")
        return out

    edges = step("edges between Bar 1 and Bar 2 (count)", lambda: len(junction_edges()))
    if edges:
        for mark in (1, 0):
            def fillet():
                es = junction_edges()
                model.ClearSelection2(True)
                ok = [e.Select4(i > 0, sd(mark)) for i, e in enumerate(es)]
                v0 = mass()
                f = fm.FeatureFillet3(195, 1.0 * MM, 1.0 * MM, 0, 0, 0, 0, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY)
                if f is None:
                    return ok, "no fillet"
                return ok, v0, mass()
            r = step(f"fillet 1 mm on the junction edges (mark {mark}): selected, volume before, after", fillet)
            if r and len(r) == 3:
                break
print("done")
