"""cantilever7_d4_solidworks: SolidWorks feature tree of a planar-truss TO (stl2cad, to/swtruss.py), member-driven sketch.

Feature tree (one sketch per structure; the face trusses "Truss" on the Front Plane, face layers "Base" on the
plane of their face):
  <name>         sketch: the section (lines, tangent arcs) and one construction centerline per member between its
                 joints. Joints dimensioned from the origin; each member edge parallel to its centerline at
                 "w" / 2; each chord's inner edge at "w" from the outline; one corner radius "R" (all fillet arcs
                 equal); outline edges on the origin or dimensioned from it (the trusses' side on a layer: "tb").
  <name> face n  Extruded Boss of the sketch, depth "t" ("tb" for a layer), towards the part; the second face of the
                 trusses from an offset H, reusing the same sketch (else a copy of it).
  Web n          rectangle on "Web plane" (z = "t"), Extruded Boss H - 2 "t" towards the other face.
Global variables (Equation Manager): t, tb, R and the member widths w1 ..., b1 ... (edit them there).
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
R_CORNER = 1.1            # global variable "R"
WEBS = []
N_PARAMS = 31             # Truss: joints 10, widths 3, depth 1, chord ends 2; Base: joints 8, widths 3, depth 1, chord ends 2; R 1, webs 0
STRUCTS = [
    {'name': 'Truss', 'axis': 2, 'frame': (0, 1), 'plane': 'Front Plane', 'box': (0.0, 2.4046, 120.0, 40.0), 'faces': [(0.0, 1), (20.0, -1)], 't': 3.1579, 'tvar': 't', 'joints': {'J1': (2.0, 37.2878), 'J2': (5.1222, 37.2878), 'J3': (36.25, 3.7212), 'J4': (36.5, 37.2878), 'J5': (52.1341, 3.7212), 'J6': (54.5333, 3.7212), 'J7': (83.567, 36.6375), 'J8': (113.75, 3.7212), 'J9': (117.8736, 3.7212), 'J10': (120.0, 3.7212)}, 'widths': {'w1': 2.6324, 'w2': 4.2525, 'w3': 5.4243}, 'members': [{'a': 'J7', 'c': 'J9', 'w': 'w2'}, {'a': 'J5', 'c': 'J2', 'w': 'w2'}, {'a': 'J9', 'c': 'J10', 'w': 'w1', 'edge': ((117.5, 2.4046), (0.0, 1.0)), 'ext': (0.0, 0.5)}, {'a': 'J1', 'c': 'J2', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (2.5, 0.0)}, {'a': 'J2', 'c': 'J4', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (5.6222, 0.0)}, {'a': 'J4', 'c': 'J7', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (0.0, 0.0)}, {'a': 'J3', 'c': 'J5', 'w': 'w1', 'edge': ((40.0, 2.4046), (0.0, 1.0)), 'ext': (0.0, 0.0)}, {'a': 'J8', 'c': 'J9', 'w': 'w1', 'edge': ((117.5, 2.4046), (0.0, 1.0)), 'ext': (0.0, 2.6264)}, {'a': 'J5', 'c': 'J6', 'w': 'w1', 'edge': ((40.0, 2.4046), (0.0, 1.0)), 'ext': (0.0, 0.0)}, {'a': 'J6', 'c': 'J7', 'w': 'w2'}], 'links': {0: 'tb'}, 'loops': [[('L', (113.41156782499391, 5.055733330433753), (84.52258580871518, 32.773968268971394), ('M', 0, -1)), ('A', (84.52258580871518, 32.773968268971394), (83.71523257686836, 33.07928049607934), (82.93607085022813, 32.70787532549851), ('A',)), ('L', (82.93607085022813, 32.70787532549851), (56.2071795465905, 2.4045999999999834), ('M', 9, -1)), ('L', (56.2071795465905, 2.4045999999999834), (36.25, 2.4045999999999994), ('D', 0)), ('L', (36.25, 2.4045999999999994), (36.25, 5.037), ('E', 2)), ('L', (36.25, 5.037), (43.19852712273649, 5.036999999999999), ('C', 2)), ('A', (43.19852712273649, 5.036999999999999), (44.246083692248014, 5.8014031679630245), (43.837721061113925, 7.032226847866813), ('A',)), ('L', (43.837721061113925, 7.032226847866813), (5.548343211406177, 34.37092684786681), ('M', 1, 1)), ('A', (5.548343211406177, 34.37092684786681), (5.244746105065721, 34.52325656951152), (4.909149273028746, 34.5757), ('A',)), ('L', (4.909149273028746, 34.5757), (0.0, 34.5757), ('C', 3)), ('L', (0.0, 34.5757), (0.0, 40.0), ('D', 3)), ('L', (0.0, 40.0), (83.13360821624222, 40.0), ('D', 2)), ('L', (83.13360821624222, 40.0), (119.573355569, 5.037), ('M', 0, 1)), ('L', (119.573355569, 5.037), (120.0, 5.037), ('C', 2)), ('L', (120.0, 5.037), (120.0, 2.4046), ('D', 1)), ('L', (120.0, 2.4046), (113.75, 2.4045999999999994), ('D', 0)), ('L', (113.75, 2.4045999999999994), (113.75, 4.261998835823597), ('E', 2)), ('A', (113.75, 4.261998835823597), (113.6618608124375, 4.693435617117619), (113.41156782499391, 5.055733330433753), ('A',))], [('L', (53.93701906461931, 6.259503026695434), (77.30115086438352, 32.74805844886272), ('M', 9, 1)), ('A', (77.30115086438352, 32.74805844886272), (77.478801782517, 33.9282451876603), (76.47620373089038, 34.5757), ('A',)), ('L', (76.47620373089038, 34.5757), (16.013369094038524, 34.5757), ('C', 3)), ('A', (16.013369094038524, 34.5757), (14.965812524527003, 33.81129683203697), (15.374175155661092, 32.58047315213319), ('A',)), ('L', (15.374175155661092, 32.58047315213319), (52.472877992748735, 6.091917729965894), ('M', 1, -1)), ('A', (52.472877992748735, 6.091917729965894), (53.237160980582566, 5.894280117020617), (53.93701906461931, 6.259503026695434), ('A',))]], 'chord_side': {2: 0, 3: 2, 4: 2, 5: 2, 6: 0, 7: 0, 8: 0}},
    {'name': 'Base', 'axis': 1, 'frame': (0, 2), 'plane': 'Top Plane', 'box': (0.0, 0.0, 120.0, 20.0), 'faces': [(0.0, 1)], 't': 2.4046, 'tvar': 'tb', 'joints': {'B1': (10.25, 9.875), 'B2': (43.25, 9.9375), 'B3': (68.6669, 10.0329), 'B4': (118.4583, 3.2623), 'B5': (118.4583, 16.5776)}, 'widths': {'b1': 3.0833, 'b2': 9.0, 'b3': 20.25}, 'members': [{'a': 'B4', 'c': 'B3', 'w': 'b1'}, {'a': 'B2', 'c': 'B3', 'w': 'b2'}, {'a': 'B3', 'c': 'B5', 'w': 'b1'}, {'a': 'B4', 'c': 'B5', 'w': 'b1', 'edge': ((120.0, 10.0), (-1.0, -0.0)), 'ext': (0.0, 0.0)}, {'a': 'B1', 'c': 'B2', 'w': 'b3', 'edge': ((26.75, 20.0), (0.0, -1.0)), 'ext': (10.75, 12.25)}], 'links': {}, 'loops': [[('L', (120.0, 1.4968231467041562), (74.6085346557402, 7.6691230872530385), ('M', 0, 1)), ('A', (74.6085346557402, 7.6691230872530385), (74.33272618364276, 7.6717285920999725), (74.06494008713014, 7.605639951239633), ('A',)), ('L', (74.06494008713014, 7.605639951239633), (68.68379021788478, 5.532931697914503), ('B',)), ('L', (68.68379021788478, 5.532931697914503), (56.59587128010077, 5.487560804185265), ('M', 1, -1)), ('A', (56.59587128010077, 5.487560804185265), (55.820724184935415, 5.1639249166285595), (55.5, 4.387568552543467), ('A',)), ('L', (55.5, 4.387568552543467), (55.5, 0.0), ('E', 4)), ('L', (55.5, 0.0), (0.0, 0.0), ('D', 0)), ('L', (0.0, 0.0), (0.0, 20.0), ('D', 3)), ('L', (0.0, 20.0), (55.5, 20.0), ('D', 2)), ('L', (55.5, 20.0), (55.5, 15.5876474454736), ('E', 4)), ('A', (55.5, 15.5876474454736), (55.82364363591491, 14.808371630409013), (56.60412871989923, 14.4876551938318), ('A',)), ('L', (56.60412871989923, 14.4876551938318), (68.65000978211522, 14.532868302085497), ('M', 1, 1)), ('L', (68.65000978211522, 14.532868302085497), (74.0558369744691, 12.432632707508585), ('B',)), ('A', (74.0558369744691, 12.432632707508585), (74.32257986087943, 12.36586995470986), (74.59754685035062, 12.367348911415114), ('A',)), ('L', (74.59754685035062, 12.367348911415114), (120.0, 18.335155302247873), ('M', 2, 1)), ('L', (120.0, 18.335155302247873), (120.0, 1.4968231467041562), ('D', 1))], [('L', (115.67334644823431, 14.656628494949217), (88.45586052454145, 11.079097432862337), ('M', 2, -1)), ('A', (88.45586052454145, 11.079097432862337), (87.49921680800703, 9.990929928554229), (88.45100082288366, 8.89850926468785), ('A',)), ('L', (88.45100082288366, 8.89850926468785), (115.66848674657652, 5.197494426885438), ('M', 0, -1)), ('A', (115.66848674657652, 5.197494426885438), (116.54022105056943, 5.458900226341508), (116.91669999999999, 6.287463615215837), ('A',)), ('L', (116.91669999999999, 6.287463615215837), (116.91669999999999, 13.566009515105126), ('C', 3)), ('A', (116.91669999999999, 13.566009515105126), (116.54206580187439, 14.392958397123309), (115.67334644823431, 14.656628494949217), ('A',))]], 'chord_side': {3: 1, 4: 2}},
]
STEPS = {
    'Truss face 1': (3500.7075, (55.036, 26.1101, 1.5789)),
    'Truss face 2': (7001.4150, (55.036, 26.1101, 10.0)),
    'Base': (10816.2178, (51.7681, 17.3254, 9.9977)),
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


def model_pt(S, p, w=0.0):
    """Model point (mm) of the structure's local point p = (u, v) at a distance w along its normal."""
    q = [0.0, 0.0, 0.0]
    u, v = S["frame"]
    q[u], q[v], q[S["axis"]] = p[0], p[1], w
    return tuple(q)


def draw_loops(sm, T, S, loops):
    """Lines and arcs of the section (local mm); returns [[(entry, object)]] per loop."""
    P = lambda p: to_sk(T, model_pt(S, p))
    out = []
    for loop in loops:
        objs = []
        for e in loop:
            if e[0] == "L":
                g = sm.CreateLine(*P(e[1]), *P(e[2]))
            else:
                c = arc_centre(e[1], e[2], e[3])
                a, b, d = (P(q)[:2] for q in (e[1], e[2], e[3]))      # direction in the sketch's own axes
                cross = (b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0])
                g = sm.CreateArc(*P(c), *P(e[1]), *P(e[3]), 1 if cross > 0 else -1)   # +1: counter-clockwise
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


def point_at(objs, xy, tol=1e-9):
    """An end point of one of objs at the sketch point xy (m)."""
    for e, g in objs:
        for nm in ("GetStartPoint2", "GetEndPoint2"):
            q = call0(g, nm)
            if q is not None and math.hypot(q.X - xy[0], q.Y - xy[1]) < tol:
                return q
    return None


def equation(m, text):
    try:
        em = call0(m, "GetEquationMgr")
        em.Add2(-1, text, True); LOG["eqs"] += 1
    except Exception as ex:
        LOG["eqs_failed"].append(f"{text} ({type(ex).__name__})")


def truss_sketch(m, S, dimension=True):
    """The structure's sketch on its plane: the section, construction centerlines, relations and dimensions.
    Returns the sketch feature and [(dimension name, expression)] to link to the global variables."""
    sm, T = open_sketch(m, S["plane"])
    E = draw_loops(sm, T, S, S["loops"])
    flat = [x for loop in E for x in loop]
    skp = lambda p: to_sk(T, model_pt(S, p))[:2]                     # sketch point (m) of a local point
    txt = lambda p: tuple(c / MM for c in skp(p))                    # (mm, for dimension text)
    o, ex = skp((0.0, 0.0)), skp((1.0, 0.0))
    swap = abs(ex[1] - o[1]) > abs(ex[0] - o[0])                     # local u along the sketch's y axis
    hv = lambda k: {"H": "V", "V": "H", "sgHORIZONTAL2D": "sgVERTICAL2D", "sgVERTICAL2D": "sgHORIZONTAL2D",
                    "sgHORIZONTALPOINTS2D": "sgVERTICALPOINTS2D",
                    "sgVERTICALPOINTS2D": "sgHORIZONTALPOINTS2D"}[k] if swap else k
    J, MB, box = S["joints"], S["members"], S["box"]
    cl = {}
    for i, mb in enumerate(MB):               # construction centerlines (not chords)
        if mb.get("edge") is None:
            cl[i] = sm.CreateCenterLine(*to_sk(T, model_pt(S, J[mb["a"]])), *to_sk(T, model_pt(S, J[mb["c"]])))
    sm.AddToDB = False
    links = []
    if dimension:
        # joints: centerline ends coincident; one end per joint dimensioned from the origin
        for j, p in J.items():
            pts = [point_at([(None, cl[i])], skp(p)) for i in cl if j in (MB[i]["a"], MB[i]["c"])]
            pts = [q for q in pts if q is not None]
            for q in pts[1:]:
                rel(m, [pts[0], q], "sgCOINCIDENT", f"{S['name']} joint {j}")
            if pts:
                for k, kind, relk in ((0, "H", "sgVERTICALPOINTS2D"), (1, "V", "sgHORIZONTALPOINTS2D")):
                    if abs(p[k]) > 1e-9:
                        q2 = (p[0] / 2, p[1] + 4) if k == 0 else (p[0] + 4, p[1] / 2)
                        dim(m, [pts[0], "ORIGIN"], hv(kind), txt(q2), f"{S['name']} joint {j} {'uv'[k]}")
                    else:
                        rel(m, [pts[0], "ORIGIN"], hv(relk), f"{S['name']} joint {j} on the axis")
        # outline edges: along the local axes; the first piece of each side on the origin or dimensioned (linked
        # to a layer's thickness where the side was cut by that layer), the others collinear with it
        first = {}
        for e, g in flat:
            if e[0] == "L" and e[3][0] == "D":
                side = e[3][1]
                rel(m, [g], hv("sgHORIZONTAL2D" if side in (0, 2) else "sgVERTICAL2D"), f"{S['name']} outline {side}")
                if side in first:
                    rel(m, [first[side], g], "sgCOLINEAR", f"{S['name']} outline {side} piece")
                    continue
                first[side] = g
                if side == 0 and abs(box[1]) < 1e-9 or side == 3 and abs(box[0]) < 1e-9:
                    rel(m, [g, "ORIGIN"], "sgCOINCIDENT", f"{S['name']} outline {side} on the origin")
                else:
                    q = call0(g, "GetStartPoint2")
                    nm = dim(m, [q, "ORIGIN"], hv("V" if side in (0, 2) else "H"), txt((e[1][0] + 4, e[1][1] + 4)),
                             f"{S['name']} outline {side}")
                    if nm and side in S["links"]:
                        links.append((nm, '"%s"' % S["links"][side]))
        # member edges: parallel to the centerline at w/2; chord inner edges parallel to the outline at w
        for e, g in flat:
            if e[0] != "L":
                continue
            tg = e[3]
            mid = ((e[1][0] + e[2][0]) / 2 + 2, (e[1][1] + e[2][1]) / 2 + 2)
            if tg[0] == "M":
                i = tg[1]
                rel(m, [g, cl[i]], "sgPARALLEL", f"{S['name']} member {i} edge")
                nm = dim(m, [g, cl[i]], "any", txt(mid), f"{S['name']} member {i} edge offset")
                if nm:
                    links.append((nm, '"%s" / 2' % MB[i]["w"]))
            elif tg[0] == "C":
                i = tg[1]
                side = S["chord_side"][i]
                if side in first:
                    rel(m, [g, first[side]], "sgPARALLEL", f"{S['name']} chord {i} edge")
                    nm = dim(m, [g, first[side]], "any", txt(mid), f"{S['name']} chord {i} width")
                    if nm:
                        links.append((nm, '"%s"' % MB[i]["w"]))
            elif tg[0] == "E":
                i = tg[1]
                side = S["chord_side"][i]
                rel(m, [g], hv("sgVERTICAL2D" if side in (0, 2) else "sgHORIZONTAL2D"), f"{S['name']} chord {i} end")
                q = call0(g, "GetStartPoint2")
                dim(m, [q, "ORIGIN"], hv("H" if side in (0, 2) else "V"), txt((e[1][0] + 2, e[1][1] - 4)),
                    f"{S['name']} chord {i} end")
        # fillet arcs: tangent to their neighbours, equal, one radius
        arc0 = None
        for loop in E:
            for k, (e, g) in enumerate(loop):
                if e[0] != "A":
                    continue
                for nb in (loop[k - 1], loop[(k + 1) % len(loop)]):
                    if nb[0][0] == "L":
                        rel(m, [g, nb[1]], "sgTANGENT", f"{S['name']} fillet tangency")
                if arc0 is None:
                    arc0 = g
                    nm = dim(m, [g], "any", txt((e[2][0] + 3, e[2][1] + 3)), f"{S['name']} corner radius")
                    if nm:
                        links.append((nm, '"R"'))
                else:
                    rel(m, [arc0, g], "sgSAMELENGTH", f"{S['name']} equal radius")
    m.ClearSelection2(True)
    sm.InsertSketch(True)
    f = m.FeatureByPositionReverse(0)
    f.Name = S["name"] if dimension else S["name"] + " copy"
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
    S = STRUCTS[0]
    x0, y0, x1, y1 = box
    sm, T = open_sketch(m, plane)
    z = S["faces"][0][0] + S["t"]
    c = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    ln = [sm.CreateLine(*to_sk(T, (*c[i], z)), *to_sk(T, (*c[(i + 1) % 4], z))) for i in range(4)]
    sm.AddToDB = False
    pts = [call0(g, "GetStartPoint2") for g in ln]
    for i, g in enumerate(ln):
        rel(m, [g], "sgHORIZONTAL2D" if i % 2 == 0 else "sgVERTICAL2D", f"web {k} side {i}")
    dim(m, [ln[0]], "any", ((x0 + x1) / 2, y0 - 3), f"web {k} length")
    dim(m, [ln[3]], "any", (x0 - 3, (y0 + y1) / 2), f"web {k} height")
    # position: a web at the domain's end is placed by that end (on the origin, or dimensioned from it)
    BOX = S["box"]
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
        f = boss(m, sk, abs(S["faces"][1][0] - S["faces"][0][0]) - 2 * S["t"], reverse=rev)
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


def boss_checked(m, sk, S, i, name):
    """Extruded Boss of the structure's face i: blind t from its start, towards the part; the direction is
    checked by the volume (the reverse is tried when it differs)."""
    w0, d = S["faces"][i]
    off = None if abs(w0) < 1e-9 else abs(w0)
    for rev in ((d < 0, d >= 0) if off is None else (True, False)):
        f = boss(m, sk, S["t"], offset=off, reverse=rev)
        if f is None:
            continue
        f.Name = name
        v, c = mass(m)
        if abs(v - STEPS[name][0]) <= TOL * STEPS[name][0]:
            return f
        m.ClearSelection2(True); f.Select2(False, 0); m.Extension.DeleteSelection2(0)
    return None


def build():
    m = new_part()
    ok = True
    links, depth = [], []
    for S in STRUCTS:
        sk, lk = truss_sketch(m, S)
        links += lk
        print(f"{S['name']} sketch: {status(sk)}; {LOG['rels']} relations, {LOG['dims']} dimensions so far")
        for i in range(len(S["faces"])):
            name = f"{S['name']} face {i + 1}" if len(S["faces"]) > 1 else S["name"]
            f = boss_checked(m, sk, S, i, name)
            if f is None and i > 0:           # the shared sketch refused: a copy of it (geometry only)
                print(f"{name}: the shared sketch was not used; a copy of the sketch instead")
                sk2, _ = truss_sketch(m, S, dimension=False)
                f = boss_checked(m, sk2, S, i, name)
            ok &= check(m, name, f)
            depth.append((name, S["tvar"]))
    S0 = STRUCTS[0]
    wp = plane_at(m, S0["faces"][0][0] + S0["t"], "Web plane") if WEBS else None
    for k, box in enumerate(WEBS, 1):
        ok &= check(m, f"Web {k}", web(m, wp, k, box))
    # global variables and links
    for S in STRUCTS:
        equation(m, f'"{S["tvar"]}" = {S["t"]:g}mm')
    equation(m, f'"R" = {R_CORNER:g}mm')
    for S in STRUCTS:
        for w, v in S["widths"].items():
            equation(m, f'"{w}" = {v:g}mm')
    for nm, expr in links:
        equation(m, f'"{nm}" = {expr}')
    for nm, tv in depth:
        equation(m, f'"D1@{nm}" = "{tv}"')
    if WEBS:
        H0 = abs(S0["faces"][1][0] - S0["faces"][0][0])
        equation(m, f'"D1@Web plane" = "{S0["tvar"]}"')
        for k in range(1, len(WEBS) + 1):
            equation(m, f'"D1@Web {k}" = {H0:g}mm - 2 * "{S0["tvar"]}"')
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
