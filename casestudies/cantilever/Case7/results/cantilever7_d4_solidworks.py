"""cantilever7_d4_solidworks: SolidWorks feature tree of a planar-truss TO (stl2cad, to/swtruss.py), member-driven sketch.

Feature tree:
  Truss          sketch on the Front Plane: the truss section (lines, tangent arcs) and one construction centerline
                 per member between its joints. Joints dimensioned from the origin; each member edge parallel to
                 its centerline at "w_k" / 2; each chord's inner edge at "w_k" from the outline; one corner radius
                 "R" (all fillet arcs equal); outline edges on the origin or dimensioned from it.
  Truss face 1   Extruded Boss, depth "t", from the Front Plane.
  Truss face 2   Extruded Boss of the same sketch from an offset H, reversed, depth "t" (else a copy of the sketch).
  Web n          rectangle on "Web plane" (z = "t"), Extruded Boss H - 2 "t" towards the other face.
Global variables (Equation Manager): t, R and the member widths w1 ... (edit them there).
After every feature the volume and centre of mass are compared with CadQuery's values (STEPS).

Needs: Windows, SolidWorks, Python with pywin32.   Run: python cantilever7_d4_solidworks.py [C:\\path\\part.SLDPRT]
"""
import math
import sys
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

sys.stdout.reconfigure(line_buffering=True)
NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
TOL = 2e-3                    # a feature's volume may differ from CadQuery's by this fraction
H = 20.000000                 # domain depth (mm)
BOX = (0.0, 0.0, 119.99999731779099, 39.99999910593033)
T_DEPTH = 3.232000            # global variable "t"
R_CORNER = 1.3            # global variable "R"
WIDTHS = {'w1': 3.2286, 'w2': 4.2525, 'w3': 5.4243}
JOINTS = {'J1': (0.75, 1.5827), 'J2': (2.0, 37.341), 'J3': (5.5256, 37.1913), 'J4': (36.5, 37.2924), 'J5': (37.5, 1.6132), 'J6': (53.4698, 2.2013), 'J7': (79.25, 1.5827), 'J8': (83.4716, 36.6064), 'J9': (119.9438, 1.704)}
MEMBERS = [{'a': 'J8', 'c': 'J9', 'w': 'w2'}, {'a': 'J1', 'c': 'J5', 'w': 'w1', 'edge': ((19.25, 0.0), (0.0, 1.0)), 'ext': (1.25, 0.0)}, {'a': 'J5', 'c': 'J6', 'w': 'w1', 'edge': ((19.25, 0.0), (0.0, 1.0)), 'ext': (0.0, 0.0)}, {'a': 'J2', 'c': 'J3', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (2.5, 0.0)}, {'a': 'J3', 'c': 'J4', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (0.0, 0.0)}, {'a': 'J4', 'c': 'J8', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (0.0, 0.0)}, {'a': 'J3', 'c': 'J6', 'w': 'w2'}, {'a': 'J9', 'c': 'J7', 'w': 'w1', 'edge': ((97.75, 0.0), (0.0, 1.0)), 'ext': (0.5562, 0.0)}, {'a': 'J6', 'c': 'J8', 'w': 'w2'}]
CHORD_SIDE = {1: 0, 2: 0, 3: 2, 4: 2, 5: 2, 7: 0}            # outline side of each chord (0 y0, 1 x1, 2 y1, 3 x0)
WEBS = [(0.0, 0.0, 91.875, 1.7886), (116.875, 0.0, 120.0, 2.5993)]
N_PARAMS = 19             # joints 8, widths 3, R and t, chord ends 1, webs 5
LOOPS = [
    [('L', (79.25, 3.2286), (112.03654450143272, 3.2286), ('C', 1)), ('A', (112.03654450143272, 3.2286), (113.24298513893078, 4.044346979159362), (112.93534838323785, 5.467829248933409), ('A',)), ('L', (112.93534838323785, 5.467829248933409), (84.59422035961109, 32.58912675036824), ('M', 0, -1)), ('A', (84.59422035961109, 32.58912675036824), (83.63677326560713, 32.948574125411334), (82.71561928808629, 32.50429639364939), ('A',)), ('L', (82.71561928808629, 32.50429639364939), (54.37135057376806, 0.0), ('M', 8, -1)), ('L', (54.37135057376806, 0.0), (0.0, 0.0), ('D', 0)), ('L', (0.0, 0.0), (0.0, 3.2286), ('D', 3)), ('L', (0.0, 3.2286), (44.46883983274119, 3.2286), ('C', 1)), ('A', (44.46883983274119, 3.2286), (45.704784136848794, 4.125559459676862), (45.23520238679786, 5.5786897274708185), ('A',)), ('L', (45.23520238679786, 5.5786897274708185), (5.845184901230553, 34.32578972747081), ('M', 6, -1)), ('A', (5.845184901230553, 34.32578972747081), (5.481862887497022, 34.511644304107605), (5.078822347173883, 34.5757), ('A',)), ('L', (5.078822347173883, 34.5757), (0.0, 34.5757), ('C', 3)), ('L', (0.0, 34.5757), (0.0, 39.99999910593033), ('D', 3)), ('L', (0.0, 39.99999910593033), (83.00070468233181, 39.99999910593033), ('D', 2)), ('L', (83.00070468233181, 39.99999910593033), (119.99999731779099, 4.593193065169143), ('M', 0, 1)), ('L', (119.99999731779099, 4.593193065169143), (119.99999731779099, 0.0), ('D', 1)), ('L', (119.99999731779099, 0.0), (79.25, 0.0), ('D', 0)), ('L', (79.25, 0.0), (79.25, 3.2286), ('E', 1))],
    [('L', (53.93189503949535, 5.966384962531873), (77.00100697948878, 32.42130110778544), ('M', 8, 1)), ('A', (77.00100697948878, 32.42130110778544), (77.20457756013448, 33.8138827942814), (76.02120978976913, 34.5757), ('A',)), ('L', (76.02120978976913, 34.5757), (16.702889771717224, 34.5757), ('C', 3)), ('A', (16.702889771717224, 34.5757), (15.46694546760962, 33.67874054032314), (15.936527217660554, 32.22561027252918), ('A',)), ('L', (15.936527217660554, 32.22561027252918), (52.18573529571903, 5.770694127275615), ('M', 6, 1)), ('A', (52.18573529571903, 5.770694127275615), (53.09688152176937, 5.5288714396330745), (53.93189503949535, 5.966384962531873), ('A',))],
]
STEPS = {
    'Truss face 1': (4453.8769, (55.8384, 21.2335, 1.616)),
    'Truss face 2': (8907.7538, (55.8384, 21.2335, 10.0)),
    'Web 1': (11132.0924, (53.86, 17.1694, 10.0)),
    'Web 2': (11242.0428, (54.4916, 17.0142, 10.0)),
}


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
    app.SetUserPreferenceToggle(10, False)        # no dimension input box
except Exception:
    pass
CONST = {}
for enum, member, default in (("swEndConditions_e", "swEndCondMidPlane", 6), ("swEndConditions_e", "swEndCondBlind", 0)):
    try:
        CONST[member] = sw_constant(app, enum, member)
    except Exception:
        CONST[member] = default
LOG = dict(rels=0, rels_failed=[], dims=0, dims_failed=[], eqs=0, eqs_failed=[])


def new_part():
    m = call0(app, "NewPart")
    try:
        m.Extension.SetUserPreferenceInteger(sw_constant(app, "swUserPreferenceIntegerValue_e", "swUnitSystem"),
                                             0, sw_constant(app, "swUnitSystem_e", "swUnitSystem_MMGS"))
    except Exception:
        pass
    return m


def mass(m):
    m.ForceRebuild3(False)
    mp = call0(m.Extension, "CreateMassProperty")
    if mp is None:
        return 0.0, (0.0, 0.0, 0.0)
    return mp.Volume / MM ** 3, tuple(c / MM for c in mp.CenterOfMass)


def check(m, step, f):
    ev, ec = STEPS[step]
    if f is None:
        print(f"{step}: NOT MADE")
        return False
    v, c = mass(m)
    ok = abs(v - ev) <= TOL * ev
    print(f"{step}: volume {v:.3f} mm^3 (expected {ev:.3f}); centre {math.dist(c, ec):.4f} mm off -> "
          f"{'OK' if ok else 'DIFFERS'}")
    return ok


def to_sk(T, p):
    x = [c * MM for c in p]
    return tuple((sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3))


def open_sketch(m, plane):
    sm = m.SketchManager
    for attempt in range(3):              # the first sketch call after NewPart can open nothing
        m.ClearSelection2(True)
        if isinstance(plane, str):
            m.Extension.SelectByID2(plane, "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
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
    helper = sm.CreateCircleByRadius(-0.5, -0.5, 0.0, 0.5 * MM)   # the first entity of a new sketch can be
    m.ClearSelection2(True)                                        # lost: a helper, deleted again
    if helper is not None and helper.Select4(False, NOTHING):
        m.EditDelete()
    return sm, T


def close_sketch(m, sm, name):
    sm.AddToDB = False
    m.ClearSelection2(True)
    sm.InsertSketch(True)
    f = m.FeatureByPositionReverse(0)
    f.Name = name
    return f


def arc_centre(s, mid, t):
    (ax, ay), (bx, by), (cx, cy) = s, mid, t
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    return ux, uy


def draw_loops(sm, T, loops):
    """Lines and arcs of the section (model mm); returns [[(entry, object)]] per loop."""
    out = []
    for loop in loops:
        objs = []
        for e in loop:
            if e[0] == "L":
                g = sm.CreateLine(*to_sk(T, (*e[1], 0.0)), *to_sk(T, (*e[2], 0.0)))
            else:
                c = arc_centre(e[1], e[2], e[3])
                cross = (e[2][0] - e[1][0]) * (e[3][1] - e[1][1]) - (e[2][1] - e[1][1]) * (e[3][0] - e[1][0])
                g = sm.CreateArc(*to_sk(T, (*c, 0.0)), *to_sk(T, (*e[1], 0.0)), *to_sk(T, (*e[3], 0.0)),
                                 1 if cross > 0 else -1)          # +1: counter-clockwise
            objs.append((e, g))
        out.append(objs)
    return out


def sel(m, objs):
    m.ClearSelection2(True)
    ok = True
    for k, o in enumerate(objs):
        if o == "ORIGIN":
            ok = ok and m.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0, k > 0, 0, NOTHING, 0)
        else:
            ok = ok and o is not None and o.Select4(k > 0, NOTHING)
    return ok


def rel(m, objs, r, what):
    try:
        if sel(m, objs):
            m.SketchAddConstraints(r); LOG["rels"] += 1
            return True
    except Exception as ex:
        what += f" ({type(ex).__name__})"
    LOG["rels_failed"].append(what)
    return False


def dim(m, objs, kind, text, what):
    try:
        d = None
        if sel(m, objs):
            tx, ty = text
            d = {"any": m.AddDimension2, "H": m.AddHorizontalDimension2, "V": m.AddVerticalDimension2}[kind](
                tx * MM, ty * MM, 0.0)
        if d is not None:
            LOG["dims"] += 1
            return call0(d.GetDimension2(0), "FullName")
    except Exception as ex:
        what += f" ({type(ex).__name__})"
    LOG["dims_failed"].append(what)
    return None


def point_at(objs, p, tol=1e-6):
    for e, g in objs:
        for nm in ("GetStartPoint2", "GetEndPoint2"):
            q = call0(g, nm)
            if q is not None and math.hypot(q.X / MM - p[0], q.Y / MM - p[1]) < tol:
                return q
    return None


def equation(m, text):
    try:
        em = call0(m, "GetEquationMgr")
        em.Add2(-1, text, True); LOG["eqs"] += 1
    except Exception as ex:
        LOG["eqs_failed"].append(f"{text} ({type(ex).__name__})")


def truss_sketch(m, dimension=True):
    sm, T = open_sketch(m, "Front Plane")
    E = draw_loops(sm, T, LOOPS)
    flat = [x for loop in E for x in loop]
    cl = {}
    for i, mb in enumerate(MEMBERS):          # construction centerlines (not chords)
        if mb.get("edge") is None:
            a, c = JOINTS[mb["a"]], JOINTS[mb["c"]]
            cl[i] = sm.CreateCenterLine(*to_sk(T, (*a, 0.0)), *to_sk(T, (*c, 0.0)))
    sm.AddToDB = False
    links = []
    if dimension:
        # joints: centerline ends coincident; one end per joint dimensioned from the origin
        for j, p in JOINTS.items():
            pts = [point_at([(None, cl[i])], p) for i in cl if j in (MEMBERS[i]["a"], MEMBERS[i]["c"])]
            pts = [q for q in pts if q is not None]
            for q in pts[1:]:
                rel(m, [pts[0], q], "sgCOINCIDENT", f"joint {j}")
            if pts:
                if abs(p[0]) > 1e-9:
                    dim(m, [pts[0], "ORIGIN"], "H", (p[0] / 2, p[1] + 4), f"joint {j} x")
                else:
                    rel(m, [pts[0], "ORIGIN"], "sgVERTICALPOINTS2D", f"joint {j} on x = 0")
                if abs(p[1]) > 1e-9:
                    dim(m, [pts[0], "ORIGIN"], "V", (p[0] + 4, p[1] / 2), f"joint {j} y")
                else:
                    rel(m, [pts[0], "ORIGIN"], "sgHORIZONTALPOINTS2D", f"joint {j} on y = 0")
        # outline edges: horizontal / vertical; the first piece of each side on the origin or dimensioned,
        # the others collinear with it
        first = {}
        for e, g in flat:
            if e[0] == "L" and e[3][0] == "D":
                side = e[3][1]
                rel(m, [g], "sgHORIZONTAL2D" if side in (0, 2) else "sgVERTICAL2D", f"outline side {side}")
                if side in first:
                    rel(m, [first[side], g], "sgCOLINEAR", f"outline side {side} piece")
                    continue
                first[side] = g
                if side == 0 and abs(BOX[1]) < 1e-9 or side == 3 and abs(BOX[0]) < 1e-9:
                    rel(m, [g, "ORIGIN"], "sgCOINCIDENT", f"outline side {side} on the origin")
                else:
                    q = call0(g, "GetStartPoint2")
                    dim(m, [q, "ORIGIN"], "V" if side in (0, 2) else "H", (e[1][0] + 4, e[1][1] + 4),
                        f"outline side {side}")
        # member edges: parallel to the centerline at w/2; chord inner edges parallel to the outline at w
        for e, g in flat:
            if e[0] != "L":
                continue
            tg = e[3]
            if tg[0] == "M":
                i = tg[1]
                rel(m, [g, cl[i]], "sgPARALLEL", f"member {i} edge")
                mid = ((e[1][0] + e[2][0]) / 2, (e[1][1] + e[2][1]) / 2)
                nm = dim(m, [g, cl[i]], "any", (mid[0] + 2, mid[1] + 2), f"member {i} edge offset")
                if nm:
                    links.append((nm, '"%s" / 2' % MEMBERS[i]["w"]))
            elif tg[0] == "C":
                i = tg[1]
                side = CHORD_SIDE[i]
                if side in first:
                    rel(m, [g, first[side]], "sgPARALLEL", f"chord {i} edge")
                    mid = ((e[1][0] + e[2][0]) / 2, (e[1][1] + e[2][1]) / 2)
                    nm = dim(m, [g, first[side]], "any", (mid[0] + 2, mid[1] + 2), f"chord {i} width")
                    if nm:
                        links.append((nm, '"%s"' % MEMBERS[i]["w"]))
            elif tg[0] == "E":
                i = tg[1]
                side = CHORD_SIDE[i]
                rel(m, [g], "sgVERTICAL2D" if side in (0, 2) else "sgHORIZONTAL2D", f"chord {i} end")
                q = call0(g, "GetStartPoint2")
                dim(m, [q, "ORIGIN"], "H" if side in (0, 2) else "V", (e[1][0] + 2, e[1][1] - 4), f"chord {i} end")
        # fillet arcs: tangent to their neighbours, equal, one radius
        arc0 = None
        for loop in E:
            for k, (e, g) in enumerate(loop):
                if e[0] != "A":
                    continue
                for nb in (loop[k - 1], loop[(k + 1) % len(loop)]):
                    if nb[0][0] == "L":
                        rel(m, [g, nb[1]], "sgTANGENT", "fillet tangency")
                if arc0 is None:
                    arc0 = g
                    nm = dim(m, [g], "any", (e[2][0] + 3, e[2][1] + 3), "corner radius")
                    if nm:
                        links.append((nm, '"R"'))
                else:
                    rel(m, [arc0, g], "sgSAMELENGTH", "equal radius")
    m.ClearSelection2(True)
    sm.InsertSketch(True)
    f = m.FeatureByPositionReverse(0)
    f.Name = "Truss" if dimension else "Truss copy"
    return f, links


def boss(m, sk, depth, offset=None, reverse=False):
    m.ClearSelection2(True)
    sk.Select2(False, 0)
    t0, off = (3, offset * MM) if offset is not None else (0, 0.0)   # 3: start offset
    return m.FeatureManager.FeatureExtrusion3(True, False, reverse, CONST["swEndCondBlind"], 0, depth * MM, 0.0,
                                              False, False, False, False, 0.0, 0.0, False, False, False, False,
                                              True, True, True, t0, off, False)


def plane_at(m, z, name):
    """Plane parallel to the Front Plane at z (mm), its side checked through a sketch transform (as the bike
    hook's plane_at_z, which SolidWorks 2025 built)."""
    for flip in (0, 256):
        m.ClearSelection2(True)
        m.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
        if m.FeatureManager.InsertRefPlane(8 | flip, z * MM, 0, 0, 0, 0) is None:
            continue
        f = m.FeatureByPositionReverse(0)
        sm, T = open_sketch(m, f)
        off = to_sk(T, (0.0, 0.0, z))[2] / MM
        sm.AddToDB = False
        sm.InsertSketch(True)
        sk = m.FeatureByPositionReverse(0)
        m.ClearSelection2(True); sk.Select2(False, 0); m.Extension.DeleteSelection2(0)
        if abs(off) < 1e-6:
            f.Name = name
            return f
        m.ClearSelection2(True); f.Select2(False, 0); m.Extension.DeleteSelection2(0)
    raise RuntimeError(f"no plane at z = {z}")


def web(m, plane, k, box):
    x0, y0, x1, y1 = box
    sm, T = open_sketch(m, plane)
    z = T_DEPTH
    c = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    ln = [sm.CreateLine(*to_sk(T, (*c[i], z)), *to_sk(T, (*c[(i + 1) % 4], z))) for i in range(4)]
    sm.AddToDB = False
    pts = [call0(g, "GetStartPoint2") for g in ln]
    for i, g in enumerate(ln):
        rel(m, [g], "sgHORIZONTAL2D" if i % 2 == 0 else "sgVERTICAL2D", f"web {k} side {i}")
    dim(m, [ln[0]], "any", ((x0 + x1) / 2, y0 - 3), f"web {k} length")
    dim(m, [ln[3]], "any", (x0 - 3, (y0 + y1) / 2), f"web {k} height")
    # position: a web at the domain's end is placed by that end (on the origin, or dimensioned from it)
    q, qx = (pts[1], x1) if abs(x1 - BOX[2]) < 1e-9 and abs(x0 - BOX[0]) > 1e-9 else (pts[0], x0)
    if abs(qx) > 1e-9:
        dim(m, [q, "ORIGIN"], "H", (qx / 2, y0 - 6), f"web {k} position")
    else:
        rel(m, [q, "ORIGIN"], "sgVERTICALPOINTS2D", f"web {k} on x = 0")
    rel(m, [pts[0], "ORIGIN"], "sgHORIZONTALPOINTS2D", f"web {k} on y = 0")
    m.ClearSelection2(True)
    sm.InsertSketch(True)
    sk = m.FeatureByPositionReverse(0)
    sk.Name = f"Web {k} sketch"
    f = None
    for rev in (False, True):                 # towards +z from the plane at z = t; the side checked by volume
        f = boss(m, sk, H - 2 * T_DEPTH, reverse=rev)
        if f is None:
            continue
        f.Name = f"Web {k}"
        v, c = mass(m)
        if abs(v - STEPS[f"Web {k}"][0]) <= TOL * v:
            return f
        m.ClearSelection2(True); f.Select2(False, 0); m.Extension.DeleteSelection2(0)
    return None


def status(f):
    try:
        st = call0(call0(f, "GetSpecificFeature2"), "GetConstrainedStatus")
        return {1: "unknown", 2: "under defined", 3: "fully defined", 4: "over defined", 5: "no solution",
                6: "invalid solution"}.get(st, str(st))
    except Exception as ex:
        return f"not available ({type(ex).__name__})"


def build():
    m = new_part()
    ok = True
    sk, links = truss_sketch(m)
    print(f"Truss sketch: {status(sk)}; {LOG['rels']} relations, {LOG['dims']} dimensions")
    f1 = boss(m, sk, T_DEPTH)
    if f1 is not None:
        f1.Name = "Truss face 1"
    ok &= check(m, "Truss face 1", f1)
    f2 = boss(m, sk, T_DEPTH, offset=H, reverse=True)            # the same sketch, from z = H towards the part
    if f2 is None or not check(m, "Truss face 2", f2):
        if f2 is not None:
            m.ClearSelection2(True); f2.Select2(False, 0); m.Extension.DeleteSelection2(0)
        print("Truss face 2: the shared sketch was not used; a copy of the sketch instead")
        sk2, _ = truss_sketch(m, dimension=False)
        f2 = boss(m, sk2, T_DEPTH, offset=H, reverse=True)
        ok &= check(m, "Truss face 2", f2)
    if f2 is not None:
        f2.Name = "Truss face 2"
    wp = plane_at(m, T_DEPTH, "Web plane") if WEBS else None
    for k, box in enumerate(WEBS, 1):
        ok &= check(m, f"Web {k}", web(m, wp, k, box))
    # global variables and links
    equation(m, f'"t" = {T_DEPTH:g}mm')
    equation(m, f'"R" = {R_CORNER:g}mm')
    for w, v in WIDTHS.items():
        equation(m, f'"{w}" = {v:g}mm')
    for nm, expr in links:
        equation(m, f'"{nm}" = {expr}')
    for nm in ("D1@Truss face 1", "D1@Truss face 2"):
        equation(m, f'"{nm}" = "t"')
    if WEBS:
        equation(m, '"D1@Web plane" = "t"')
    for k in range(1, len(WEBS) + 1):
        equation(m, f'"D1@Web {k}" = {H:g}mm - 2 * "t"')
    try:
        m.EditRebuild3()
    except Exception:
        pass
    v, c = mass(m)
    ev, ec = STEPS[list(STEPS)[-1]]
    print(f"after the equations: volume {v:.3f} mm^3 (expected {ev:.3f})")
    print(f"relations: {LOG['rels']} added" + (f"; failed: {LOG['rels_failed']}" if LOG['rels_failed'] else ""))
    print(f"dimensions: {LOG['dims']} added" + (f"; failed: {LOG['dims_failed']}" if LOG['dims_failed'] else ""))
    print(f"equations: {LOG['eqs']} added" + (f"; failed: {LOG['eqs_failed']}" if LOG['eqs_failed'] else ""))
    print(f"independent parameters: {N_PARAMS}")
    print("Truss sketch: " + status(sk))
    print("all features " + ("agree" if ok else "do NOT all agree"))
    return m


if __name__ == "__main__":
    try:
        m = build()
        if len(sys.argv) > 1:
            m.SaveAs3(sys.argv[1], 0, 2)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
