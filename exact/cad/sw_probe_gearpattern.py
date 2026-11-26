"""Probe: which dimensions make the gear's circular sketch pattern fully defined in SolidWorks.

Run from exact/cad (it imports Gear_solidworks.py from there):
    python sw_probe_gearpattern.py
Each variant opens a new part, draws the first sketch of the gear (seed period, its planned
relations and dimensions, the circular pattern) and prints the sketch status. The parts stay open.
"""
import math
import types

import Gear_solidworks as G
from Gear_solidworks import SW, FEATURES, NOTHING, MM, call0, set_units

ST = FEATURES[0]
LOOP = next(lp for lp in ST["loops"] if isinstance(lp, dict) and "centre" in lp)

# (label, (RadiusDim, AngleDim, CreateNumOfInstancesDim), dimension the pattern centre to the origin)
VARIANTS = [
    ("no pattern dims", (False, False, False), False),
    ("radius+angle+count", (True, True, True), False),
    ("angle+count", (False, True, True), False),
    ("no pattern dims, centre dims", (False, False, False), True),
    ("angle+count, centre dims", (False, True, True), True),
    ("radius+angle+count, centre dims", (True, True, True), True),
]


def seed_box_centre(ents):
    xs, ys = [], []
    for kind, geo, g in ents:
        xs += [geo[0][0], geo[1][0]]
        ys += [geo[0][1], geo[1][1]]
    return (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2


def points_near(sw, x, y, tol=1e-6):
    sk = sw.model.SketchManager.ActiveSketch
    pts = call0(sk, "GetSketchPoints2") or []
    return [p for p in pts if math.hypot(p.X - x, p.Y - y) < tol]


def make_pattern(flags, centre_dims, log):
    def circular_pattern(self, ents, loop):
        bc = seed_box_centre(ents)              # the gear period is lines only
        cx, cy = self.to_sketch(*loop["centre"])
        R, A = math.hypot(cx - bc[0], cy - bc[1]), math.atan2(cy - bc[1], cx - bc[0]) % (2 * math.pi)
        n = loop["count"]
        p0, p1, q1 = self.to_sketch(0.0, 0.0), self.to_sketch(1.0, 0.0), self.to_sketch(0.0, 1.0)
        turn = 1.0 if (p1[0] - p0[0]) * (q1[1] - p0[1]) - (p1[1] - p0[1]) * (q1[0] - p0[0]) > 0 else -1.0
        log.append(f"points at the centre before: {len(points_near(self, cx, cy))}")
        self.model.ClearSelection2(True)
        ok = all(g.Select4(True, NOTHING) for _, _, g in ents)
        r = self.model.SketchManager.CreateCircularSketchStepAndRepeat(
            R, A, n, turn * 2 * math.pi / n, True, "", *flags) if ok else False
        self.model.ClearSelection2(True)
        pts = points_near(self, cx, cy)
        log.append(f"pattern returned {r}; points at the centre after: {len(pts)}")
        if not pts:                              # where did the pattern put its centre point?
            try:
                sk = self.model.SketchManager.ActiveSketch
                ends = set()
                for g in call0(sk, "GetSketchSegments") or []:
                    for nm in ("GetStartPoint2", "GetEndPoint2"):
                        q = call0(g, nm)
                        ends.add((round(q.X / MM, 4), round(q.Y / MM, 4)))
                lone = [(round(p.X / MM, 4), round(p.Y / MM, 4)) for p in call0(sk, "GetSketchPoints2") or []]
                lone = [p for p in lone if p not in ends]
                log.append(f"points that are no segment end (mm): {lone[:10]}; centre is "
                           f"({cx / MM:.4f}, {cy / MM:.4f})")
            except Exception as ex:
                log.append(f"sketch points not readable: {type(ex).__name__}")
        if not (centre_dims and pts):
            return
        ext = self.model.Extension
        for kind, tx, ty in (("H", cx / 2, cy - 30 * MM), ("V", cx - 30 * MM, cy / 2)):
            self.model.ClearSelection2(True)
            ok1 = pts[0].Select4(False, NOTHING)
            ok2 = ext.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0.0, 0.0, 0.0, True, 0, NOTHING, 0)
            d = None
            if ok1 and ok2:
                d = (self.model.AddHorizontalDimension2(tx, ty, 0.0) if kind == "H"
                     else self.model.AddVerticalDimension2(tx, ty, 0.0))
            log.append(f"centre {kind} dimension: selections {ok1}, {ok2}; dimension "
                       f"{'added' if d is not None else 'failed'}")
        self.model.ClearSelection2(True)
    return circular_pattern


results = []
for label, flags, centre_dims in VARIANTS:
    print(f"\n=== {label} ===")
    log = []
    sw = SW()
    set_units(sw)
    sw.circular_pattern = types.MethodType(make_pattern(flags, centre_dims, log), sw)
    sw.report_free = lambda *a: None
    try:
        sw.sketch(ST["loops"], True, plans=ST["plans"])
        status = sw.status[-1] if sw.status else "none"
    except Exception as ex:
        status = f"error {type(ex).__name__}: {ex}"
    for line in log:
        print("  " + line)
    print(f"  dimensions added {sw.dims_ok}, failed {len(sw.dims_failed)}; relations added {sw.rels_ok}, "
          f"failed {len(sw.rels_failed)}")
    print(f"  sketch status: {status}")
    results.append((label, status))

print("\n----- summary -----")
for label, status in results:
    print(f"{label:36s} {status}")
