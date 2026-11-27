"""Probe: a truss joint in SolidWorks (organic-truss, the 6-strut node at (49.42, 187.24, 9.82) mm).
Run on the Windows machine with SolidWorks:  python sw_probe_struts.py > probe_struts.log  and send the log.
  A. one strut: a 3D sketch line and a circular-profile sweep along it (path: the sketch, then the line);
  B. the node's 6 struts (3D sketch + sweep each, merged); expected union (CadQuery) 701.722 mm^3;
  C. fillets r 2.2 on the node's 21 re-entrant edges, picked by rays from the CadQuery model (all at once,
     then smaller radii, then edge by edge); volume before and after;
  D. a tapered strut: planes normal to a 3D sketch line at its ends, a circle on each, loft along the line;
     expected cone volume 262.5975 mm^3.
"""
DATA_JSON = r'''{"node": [49.42156807831479, 187.2385495357149, 9.818582264538595], "struts": [{"p0": [49.42156807831479, 187.2385495357149, 9.818582264538595], "p1": [82.88909834846423, 172.47102339391864, 10.062089020179592], "r": 1.37}, {"p0": [49.42156807831479, 187.2385495357149, 9.818582264538595], "p1": [17.46240858526515, 192.73319158300293, 4.2384324912460904], "r": 1.0}, {"p0": [49.42156807831479, 187.2385495357149, 9.818582264538595], "p1": [19.116551667343465, 201.44099428817896, 9.028962548471942], "r": 1.37}, {"p0": [44.67037379797362, 177.9815421014386, 20.003186208838297], "p1": [49.42156807831479, 187.2385495357149, 9.818582264538595], "r": 1.0}, {"p0": [76.64361291185674, 158.1775042075479, 7.837228743002764], "p1": [49.42156807831479, 187.2385495357149, 9.818582264538595], "r": 1.0}, {"p0": [49.42156807831479, 187.2385495357149, 9.818582264538595], "p1": [43.93547295791189, 174.41844215604377, 5.505881588062544], "r": 1.0}], "union": 701.7219403610064, "single": [215.70186798616558, 103.37305667133049, 197.39716718963257, 45.741437051764805, 125.25108117205332, 45.85560800772336], "edges": [[[49.96593708167102, 188.41952561802609, 9.387076924469111], [0.3973207952809854, 0.8619660479458515, -0.31488206971165417], 0.8782547950650578], [[49.57404415992012, 187.52512985575856, 8.487512701320291], [0.11125049106638832, 0.20908302932993958, -0.9715490801208672], 1.8938504245864487], [[49.60345524046702, 186.31375497354986, 8.690825732097146], [0.2083438883765522, -0.2843902955902595, -0.9357964436511947], 1.152542877548174], [[52.12124168815524, 185.09404958750002, 8.784226173727241], [0.05236093978006176, -0.12640584388470166, -0.9905957271345083], 5.106585063189127], [[52.01399359923084, 184.8554305952106, 10.609893032362768], [-0.0685655000426298, -0.3710010472187724, 0.9260977244149121], 5.643183792561727], [[49.54626307713396, 185.72710611827878, 10.141879624262897], [0.1071471317584955, -0.8274215461512373, 0.551264979044219], 0.09665688563668007], [[49.944098677480824, 186.29197869590718, 11.027574848621938], [0.48859693030932555, -0.3757684627024757, 0.7874459360048296], 2.360909039611035], [[49.92015411168968, 188.3570762551687, 10.432818032497645], [0.36390817304098805, 0.816394067094746, 0.44841004538978546], 1.2738440420808852], [[44.03862098859546, 189.09740967267422, 8.453550430395971], [0.03954483406082649, 0.35543396762277235, -0.9338644980718753], 7.809344570410393], [[48.03797767298997, 186.74684784610068, 8.893081529394442], [-0.15753853747982247, -0.7170902165753755, -0.6789427299121772], 2.145931816631724], [[47.9982947153208, 186.40649103768936, 9.611102293661009], [-0.921122659017829, -0.3877952888789915, 0.033880096918159], 0.3908167511414175], [[48.42201550927949, 186.68436182191232, 8.792815197632253], [-0.5041819219619147, -0.09077560231444146, -0.8588133555035344], 1.6404240805371675], [[48.26960126405018, 187.04442580554695, 10.994235830397344], [-0.6363511025980042, 0.14519028537777467, 0.7576127343532608], 2.493918616435299], [[48.11474676383302, 186.34789054038527, 9.955981536280829], [-0.804123360556369, -0.44807441159467293, 0.39065962509788393], 0.3513460161400766], [[45.440806321616336, 187.65030157962445, 10.108859955478747], [-0.3518086981025505, -0.6265085480705348, 0.695498151783266], 5.061071795547067], [[48.37047711682128, 185.74553490245614, 10.234843397110167], [-0.8478452117051668, -0.5039046706284315, 0.1650411461044574], 1.2512131663018722], [[49.22873707849764, 185.5639745048935, 10.23636962161766], [0.4457471328309734, -0.8745483177416798, 0.1909835948665866], 0.7018345605666276], [[49.742419906642375, 185.47391821596887, 9.509161615967871], [0.15858220677745202, -0.9836662002762396, -0.08516155310759975], 1.3644870404109506], [[49.532889772818976, 187.53663320347377, 10.907601480516051], [0.8404152570158211, 0.15135358232355034, 0.5203789858294239], 0.8010543456006992], [[49.54300396646069, 187.5317201271142, 10.907854344043086], [-0.441814640992606, 0.7411149819587348, 0.5055179586531807], 0.8015765067743704], [[49.70295394850002, 187.8919726517652, 10.771955485332303], [-0.3401113315048258, -0.7258596785959491, 0.5978729038605113], 0.3608738792108104]], "cone": {"p0": [176.06247255524573, 81.5237135932405, 19.56351719200165], "p1": [189.87610377145825, 101.5035220862702, 20.37337171236201], "r0": 2.3, "r1": 1.37, "volume": 262.59748243974474}}'''
import json
import math
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
EMPTY = VARIANT(pythoncom.VT_EMPTY, None)
MM = 0.001
DATA = json.loads(DATA_JSON)


def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True
try:
    app.SetUserPreferenceToggle(10, False)
except Exception:
    pass
FM_SWEEP = 17


def new_part():
    return call0(app, "NewPart")


def mass(m):
    m.ForceRebuild3(False)
    mp = call0(m.Extension, "CreateMassProperty")
    return None if mp is None else round(mp.Volume / MM ** 3, 4)


def nbodies(m):
    try:
        return len(m.GetBodies2(0, True) or [])
    except Exception:
        return None


def last(m):
    return m.FeatureByPositionReverse(0)


def line3d(m, p0, p1):
    """A 3D sketch with one line (model mm); returns (sketch feature, line). The first sketch after NewPart
    may not open, and the first line may not be made: both are retried (and reported)."""
    sm = m.SketchManager
    m.ClearSelection2(True)
    for attempt in range(3):                    # an odd number of toggles leaves the sketch open
        sm.Insert3DSketch(True)
        if call0(sm, "ActiveSketch") is not None or call0(m, "GetActiveSketch2") is not None:
            break
        print(f"    3D sketch not reported open (try {attempt + 1})")
    sm.AddToDB = True
    ln = None
    for attempt in range(4):
        ln = sm.CreateLine(*(v * MM for v in p0), *(v * MM for v in p1))
        if ln is not None:
            break
        print(f"    CreateLine returned None (try {attempt + 1})")
    sm.AddToDB = False
    sm.Insert3DSketch(True)
    return last(m), ln


def sweep(m, sel_path, d_mm):
    fm = m.FeatureManager
    m.ClearSelection2(True)
    ok = sel_path()
    d = fm.CreateDefinition(FM_SWEEP)
    d.CircularProfile = True
    d.CircularProfileDiameter = d_mm * MM
    try:
        d.Merge = True
    except Exception:
        pass
    f = fm.CreateFeature(d)
    return ok, f


def strut(m, s, how="sketch"):
    sk, ln = line3d(m, s["p0"], s["p1"])
    if how == "sketch":
        sel = lambda: sk.Select2(True, 4)
    else:
        def sel():
            sd = call0(m.SelectionManager, "CreateSelectData")
            sd.Mark = 4
            return ln.Select4(True, sd)
    ok, f = sweep(m, sel, 2 * s["r"])
    return sk, ln, ok, f


print("== A. one strut")
s = DATA["struts"][0]
L = math.dist(s["p0"], s["p1"])
exp = math.pi * s["r"] ** 2 * L
for how in ("sketch", "line"):
    try:
        m = new_part()
        sk, ln, ok, f = strut(m, s, how)
        print(f"  path = {how}: line {ln is not None}, selected {ok}, sweep {'feature' if f is not None else None}; "
              f"volume {mass(m)} (expected {exp:.4f})")
    except Exception:
        traceback.print_exc()

print("== B. the node's 6 struts; expected union %.4f" % DATA["union"])
mB = None
try:
    mB = new_part()
    for k, s in enumerate(DATA["struts"]):
        sk, ln, ok, f = strut(mB, s, "sketch")
        print(f"  strut {k} (r {s['r']}): sweep {'feature' if f is not None else None}; volume {mass(mB)}, bodies {nbodies(mB)}")
except Exception:
    traceback.print_exc()


def pick(m, p, d, mark):
    ext, sm = m.Extension, m.SelectionManager
    n0 = sm.GetSelectedObjectCount2(-1)
    q = [p[i] + 0.3 * d[i] for i in range(3)]
    ok = ext.SelectByRay(*(v * MM for v in q), -d[0], -d[1], -d[2], 0.05 * MM, 1, True, mark, 0)
    return bool(ok) and sm.GetSelectedObjectCount2(-1) > n0


def fillet(m, edges, r, mark=1):
    fm = m.FeatureManager
    m.ClearSelection2(True)
    picked = sum(pick(m, p, d, mark) for p, d, _ in edges)
    f = None
    for how, mk in (("FeatureFillet3", lambda: fm.FeatureFillet3(195, r * MM, r * MM, 0, 0, 0, 0, EMPTY, EMPTY, EMPTY,
                                                                   EMPTY, EMPTY, EMPTY, EMPTY)),
                    ("FeatureFillet", lambda: fm.FeatureFillet(195, r * MM, 0, 0, EMPTY, EMPTY, EMPTY))):
        try:
            f = mk()
        except Exception as ex:
            print(f"    {how}: {type(ex).__name__}: {ex}")
            f = None
        if f is not None:
            return picked, how, f
        m.ClearSelection2(True)
        picked = sum(pick(m, p, d, mark) for p, d, _ in edges)
    return picked, None, None


print("== C. fillets on the node's %d re-entrant edges" % len(DATA["edges"]))
try:
    if mB is not None:
        v0 = mass(mB)
        done = False
        for r in (2.2, 1.5, 1.0, 0.5):
            picked, how, f = fillet(mB, DATA["edges"], r)
            print(f"  all edges, r {r}: picked {picked} of {len(DATA['edges'])}; {how or 'no feature'}; volume {mass(mB)} (before {v0})")
            if f is not None:
                done = True
                break
        if not done:
            ok = 0
            for k, e in enumerate(DATA["edges"]):
                picked, how, f = fillet(mB, [e], 2.2)
                ok += f is not None
            print(f"  edge by edge, r 2.2: {ok} of {len(DATA['edges'])} filleted; volume {mass(mB)}")
except Exception:
    traceback.print_exc()

print("== D. tapered strut: loft between circles on planes normal to a 3D sketch line; expected %.4f" % DATA["cone"]["volume"])
try:
    c = DATA["cone"]
    m = new_part()
    sk, ln = line3d(m, c["p0"], c["p1"])
    fm = m.FeatureManager
    planes = []
    for end, pt in ((0, call0(ln, "GetStartPoint2")), (1, call0(ln, "GetEndPoint2"))):
        m.ClearSelection2(True)
        sd0 = call0(m.SelectionManager, "CreateSelectData"); sd0.Mark = 0
        sd1 = call0(m.SelectionManager, "CreateSelectData"); sd1.Mark = 1
        ok = ln.Select4(False, sd0) and pt.Select4(True, sd1)
        p = fm.InsertRefPlane(2, 0, 4, 0, 0, 0)
        print(f"  plane at end {end}: selected {ok}; {'plane' if p is not None else None}")
        planes.append(last(m) if p is not None else None)
    secs = []
    for pl, ctr, r in zip(planes, (c["p0"], c["p1"]), (c["r0"], c["r1"])):
        if pl is None:
            continue
        m.ClearSelection2(True)
        pl.Select2(False, 0)
        sm = m.SketchManager
        sm.InsertSketch(True)
        T = list(call0(sm, "ActiveSketch").ModelToSketchTransform.ArrayData)
        x = [v * MM for v in ctr]
        u = [(sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3)]
        sm.AddToDB = True
        circ = sm.CreateCircleByRadius(u[0], u[1], 0.0, r * MM)
        sm.AddToDB = False
        sm.InsertSketch(True)
        secs.append(last(m))
        print(f"  section r {r}: centre in sketch ({u[0] / MM:.4f}, {u[1] / MM:.4f}, z {u[2] / MM:.4f}); circle {circ is not None}")
    if len(secs) == 2:
        m.ClearSelection2(True)
        secs[0].Select2(False, 1); secs[1].Select2(True, 1); sk.Select2(True, 4)
        f = fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False, 0.0, 0.0, 0,
                                      True, True, True, 0)
        print(f"  loft: {'feature' if f is not None else None}; volume {mass(m)} (expected {c['volume']:.4f})")
except Exception:
    traceback.print_exc()
print("done")
