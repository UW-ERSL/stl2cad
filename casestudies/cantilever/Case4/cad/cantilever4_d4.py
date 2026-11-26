"""cantilever4_d4: topology-optimized part as plates of the original part.

Original features first, then 4 plate(s); units mm.
Run:  python cantilever4_d4.py [out.step]
"""

import math
import sys

from cadquery import Edge, Solid, Vector, Vertex
from OCP.BRepOffset import BRepOffset_Analyse
from OCP.ChFiDS import ChFiDS_TypeOfConcavity
from cadquery import func as cf

# ---- parameters ----------------------------------------------------------
z0 = 0.0    # level along x
z1 = 120.0    # level along x

DIR = (1, 0, 0)


def P(a, b):
    """Point (y, z) in the sketch plane (x = 0)."""
    return (0, a, b)


def prism(outline, lo, hi):
    """Outline extruded from level lo to level hi."""
    body = cf.extrude(cf.face(outline), [(hi - lo) * c for c in DIR])
    return body.translate(Vector(*[lo * c for c in DIR]))


def hole(a, b, d, lo, hi):
    return prism(Edge.makeCircle(d / 2, P(a, b), DIR), lo, hi)


# ---- 1. extrusion: x 0.0 to 120.0 ------------------------------
outline_1 = cf.wire(
    cf.segment(P(0.0, 20.0), P(40.0, 20.0)),
    cf.segment(P(40.0, 20.0), P(40.0, 0.0)),
    cf.segment(P(40.0, 0.0), P(0.0, 0.0)),
    cf.segment(P(0.0, 0.0), P(0.0, 20.0)),
)
part = prism(outline_1, z0, z1)


# ---- topology-optimized part: plates -----------------------------------------
# Each plate is the original part cut to a slab [lo, hi] along its axis, minus cuts
# sketched in the plate plane (u, v) = the other two axes in increasing order.
# A cut is a closed polygon whose vertices carry tangent fillet radii (x, y, r)
# (r = 0: sharp), or an arc notch closed by a polyline outside the original outline.
original = part


def to3(axis, u, v, w):
    """Point (u, v) of a plate plane at level w along axis, in model coordinates."""
    p = [0.0, 0.0, 0.0]
    a, b = [i for i in range(3) if i != axis]
    p[a], p[b], p[axis] = u, v, w
    return Vector(*p)


def rounded_polygon(verts, axis, w):
    """Closed polygon with a tangent fillet of radius r at each vertex (r = 0: sharp)."""
    V = [Vector(x, y, 0) for x, y, _ in verts]
    R = [r for _, _, r in verts]
    n = len(V)
    ends = []
    for k in range(n):
        a, b, c = V[k - 1], V[k], V[(k + 1) % n]
        u = (b - a).normalized()
        q = (c - b).normalized()
        if R[k] <= 0:
            ends.append((b, b, None))
            continue
        th = math.acos(max(-1.0, min(1.0, u.dot(q))))
        t = R[k] * math.tan(th / 2)
        bis = (q - u).normalized()
        centre = b + bis * (R[k] / math.cos(th / 2))
        ends.append((b - u * t, b + q * t, centre - bis * R[k]))
    P = lambda p: to3(axis, p.x, p.y, w)
    edges = []
    for k in range(n):
        t1, t2, mid = ends[k]
        if mid is not None:
            edges.append(Edge.makeThreePointArc(P(t1), P(mid), P(t2)))
        nxt = ends[(k + 1) % n][0]
        if (nxt - t2).Length > 1e-6:
            edges.append(cf.segment(P(t2), P(nxt)))
    return cf.wire(*edges)


def arc_notch(cx, cy, r, a0, a1, path, axis, w):
    pt = lambda a: to3(axis, cx + r * math.cos(a), cy + r * math.sin(a), w)
    edges = [Edge.makeThreePointArc(pt(a0), pt((a0 + a1) / 2), pt(a1))]
    poly = [pt(a1)] + [to3(axis, x, y, w) for x, y in path] + [pt(a0)]
    edges += [cf.segment(p, q) for p, q in zip(poly[:-1], poly[1:])]
    return cf.wire(*edges)


def slab(axis, lo, hi, box):
    """Box spanning the original part, cut to [lo, hi] along axis."""
    (x0, y0, z0), (x1, y1, z1) = box
    lo_, hi_ = [x0 - 1, y0 - 1, z0 - 1], [x1 + 1, y1 + 1, z1 + 1]
    lo_[axis], hi_[axis] = lo, hi
    return Solid.makeBox(hi_[0] - lo_[0], hi_[1] - lo_[1], hi_[2] - lo_[2], Vector(*lo_))


def through(wire, axis, lo, hi):
    d = [0.0, 0.0, 0.0]
    d[axis] = hi - lo + 2
    return cf.extrude(cf.face(wire), d)


def concave_edges(solid, angle=0.01):
    """Re-entrant edges, from the two face normals at each edge (OCC BRepOffset_Analyse:
    convex, concave or tangent from the sign of (n1 x n2) . t along the edge)."""
    an = BRepOffset_Analyse(solid.wrapped, angle)
    return [e for e in solid.Edges()
            if any(i.Type() == ChFiDS_TypeOfConcavity.ChFiDS_Concave for i in an.Type(e.wrapped))]


def fillet_reentrant(solid, r):
    """Fillet every re-entrant edge with radius r (all at once; if that fails, one edge at a time)."""
    edges = concave_edges(solid)
    if not edges:
        return solid, 0, 0
    try:
        return solid.fillet(r, edges), len(edges), 0
    except Exception:
        pass
    done = skipped = 0
    for e0 in edges:
        mid = Vertex.makeVertex(*e0.positionAt(0.5).toTuple())
        cand = [e for e in concave_edges(solid) if e.distance(mid) < 1e-3]
        if not cand:
            continue
        try:
            solid = solid.fillet(r, cand[:1]); done += 1
        except Exception:
            skipped += 1
    return solid, done, skipped


def ring_cylinder(axis, lo, hi, cx, cy, r):
    """Cylinder of radius r about the axis through (cx, cy) (in-plane coordinates), from lo to hi."""
    u, v = [i for i in range(3) if i != axis]
    p = [0.0, 0.0, 0.0]; p[axis], p[u], p[v] = lo, cx, cy
    d = [0.0, 0.0, 0.0]; d[axis] = 1.0
    return Solid.makeCylinder(r, hi - lo, Vector(*p), Vector(*d))


def plate(axis, lo, hi, cuts, box, ring=None):
    """The original between lo and hi along axis, minus the cuts; with ring = (cx, cy, r): only the part
    of it within radius r of the axis through (cx, cy) (a hub around one of the original's holes)."""
    body = original.intersect(slab(axis, lo, hi, box))
    if ring is not None:
        body = body.intersect(ring_cylinder(axis, lo, hi, *ring))
    for cut in cuts:
        if cut[0] == "poly":
            body = body - through(rounded_polygon(cut[1], axis, lo - 1), axis, lo, hi)
        else:
            (cx, cy, r, a0, a1), path = cut[1], cut[2]
            body = body - through(arc_notch(cx, cy, r, a0, a1, path, axis, lo - 1), axis, lo, hi)
    return body

BOX = ((0.0, 0.0, 0.0), (120.0, 40.0, 20.0))

PLATES = [
    dict(  # normal to z, 20.0 mm thick, 4 cut(s)
        axis=2, lo=0.0000, hi=20.0000, ring=None,
        cuts=[
            ("poly", [
                (  80.9058,    7.3650,   1.500),
                (  97.6320,   12.6142,   0.250),
                (  98.3050,   11.7148,   0.250),
                (  82.4929,    4.3657,   0.250),
                ( 110.3800,    3.1956,   0.000),
                ( 110.3407,    2.6682,   0.250),
                ( 118.3987,   -0.1058,   0.250),
                ( 107.1814,   14.4537,   1.500),
                ( 104.3719,   16.0966,   1.500),
                ( 103.5419,   20.7779,   1.500),
                (  97.9408,   25.4093,   1.500),
                (  92.0113,   25.4234,   1.500),
                (  87.6116,   27.6912,   1.500),
                (  71.5721,   30.5535,  25.809),
                (  58.7706,   29.1941,   1.500),
                (  54.1926,   26.5115,   1.500),
                (  55.6989,   13.5278,   1.500),
                (  58.9206,   11.0568,   1.500),
                (  66.4165,   10.2246,  17.645),
                (  72.2472,   11.1539,  15.686),
                (  88.5898,   18.6912,   0.250),
                (  71.8597,    7.0016,   0.250),
            ]),
            ("poly", [
                (  30.6867,   13.1075,   1.500),
                (  43.4488,   15.9345,  10.507),
                (  49.3377,   19.8899,   1.500),
                (  43.4354,   24.0759,  12.290),
                (  30.9180,   26.7447,   1.500),
                (  23.4134,   19.2498,   1.500),
            ]),
            ("poly", [
                (-1089.8246,   39.9800,   0.000),
                (-1089.3246,   40.0000,   0.000),
                (  89.7124,   39.7351,  22.283),
                ( 114.7749,   23.8195,   0.000),
                ( 113.6126,   22.1053,   1.500),
                ( 119.6353,   12.7334,   1.500),
                ( 120.0000,   -1.2682,   0.000),
                ( 120.0200,   -1.7682,   0.000),
                ( 120.0200,   -4.0000,   0.000),
                ( 124.0000,   -4.0000,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  -4.0000,   44.0000,   0.000),
                (  -4.0000,   39.9800,   0.000),
            ]),
            ("poly", [
                (  -0.0441,    0.0000,   0.000),
                (  -0.0000,   14.7127,   0.000),
                (   0.0108,   14.0926,   0.000),
                (   4.0673,   13.7118,   1.500),
                (  19.4271,   18.0897,   1.500),
                (  18.0741,   22.3661,   1.500),
                (   5.0923,   26.1083,   1.500),
                (   0.1061,   25.4684,   1.500),
                (  -0.0000,   20.5316,   0.000),
                (  -0.0584,   40.0000,   0.000),
                (  -4.0000,   40.0000,   0.000),
                (  -4.0000,    0.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 18.5 mm thick, 4 cut(s)
        axis=2, lo=0.7500, hi=19.2500, ring=None,
        cuts=[
            ("poly", [
                (  85.9693,    5.0751,   0.250),
                ( 114.4663,    3.6634,   0.250),
                ( 106.3214,   14.6555,   1.500),
                ( 102.9124,   14.5150,   1.500),
                ( 102.3960,   20.4094,   1.500),
                (  96.5319,   24.7461,   1.500),
                (  91.5349,   23.9045,   1.500),
                (  84.7385,   28.2066,  12.435),
                (  55.3773,   29.7322,   1.500),
                (  57.5134,   14.1525,   1.500),
                (  63.2526,   11.0677,   1.500),
                (  72.7421,   11.9156,   1.500),
                (  85.6444,   18.6890,   0.250),
                (  86.9603,   16.8677,   0.250),
                (  74.3390,    8.3131,   0.250),
                (  80.5358,    7.7881,  12.231),
                ( 101.8947,   14.9048,   0.250),
                ( 102.4796,   13.2043,   0.250),
            ]),
            ("poly", [
                (  44.8776,   22.9108,   1.500),
                (  31.5737,   25.5828,   1.500),
                (  26.8938,   20.1661,   1.500),
                (  30.7860,   13.9762,   1.500),
                (  48.1028,   17.8827,   1.500),
            ]),
            ("poly", [
                (  69.7510,   40.0600,   0.000),
                (  89.7510,   40.0000,  17.206),
                ( 114.3972,   25.3037,   0.000),
                ( 113.5045,   23.6908,   1.500),
                ( 121.5236,   10.4483,   0.000),
                ( 124.0000,   10.4483,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  69.7510,   44.0000,   0.000),
            ]),
            ("poly", [
                (  -1.5422,   13.4198,   0.000),
                (  18.9852,   18.1180,   1.500),
                (  16.4339,   22.1588,   1.500),
                (  -1.5315,   26.9558,   0.000),
                (  -4.0000,   26.9558,   0.000),
                (  -4.0000,   13.4198,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 16.0 mm thick, 8 cut(s)
        axis=2, lo=2.0000, hi=18.0000, ring=None,
        cuts=[
            ("poly", [
                (  87.6576,    5.2545,   0.250),
                ( 113.4736,    4.4393,   1.500),
                ( 106.3587,   14.5188,   1.500),
            ]),
            ("poly", [
                (  97.6196,   22.8202,   1.500),
                (  93.3650,   22.8382,   1.500),
                (  76.3045,    8.8060,   0.250),
                (  84.1770,    9.4326,   1.500),
                (  99.5260,   15.4629,   1.500),
                ( 101.1293,   19.0539,   1.500),
            ]),
            ("poly", [
                (  80.9768,   27.3350,   0.250),
                (  80.9647,   28.1308,   0.000),
                (  62.3514,   28.2541,   1.500),
                (  64.5839,   18.6722,   1.500),
                (  62.7375,   13.0369,   0.250),
                (  70.3892,   11.8035,   1.500),
                (  87.6355,   21.0088,   1.500),
                (  87.6946,   25.3578,   1.500),
            ]),
            ("poly", [
                (  43.0365,   23.8924,   1.500),
                (  37.8072,   23.5696,   1.500),
                (  35.7390,   19.7785,   1.500),
                (  37.8218,   16.3217,   1.500),
                (  41.6984,   16.2962,   1.500),
                (  47.2130,   19.4812,   1.500),
            ]),
            ("poly", [
                (  33.8011,   14.7974,   1.500),
                (  35.7412,   16.7282,   0.250),
                (  33.7185,   17.5949,   1.500),
                (  30.7863,   15.4890,   0.250),
            ]),
            ("poly", [
                (  33.5464,   25.0498,   1.500),
                (  30.8018,   24.1537,   0.250),
                (  33.5777,   22.2795,   1.500),
                (  36.3542,   23.5417,   0.250),
            ]),
            ("poly", [
                (  -1.5234,   10.6741,   0.000),
                (   0.5190,   14.0803,   1.500),
                (  23.7353,   19.8452,   1.500),
                (   0.1339,   25.9942,   1.500),
                (  -1.5113,   32.2002,   0.000),
                (  -4.0000,   32.2002,   0.000),
                (  -4.0000,   10.6741,   0.000),
            ]),
            ("poly", [
                (  69.6203,   40.0600,   0.000),
                (  89.6203,   40.0000,  16.354),
                ( 111.2658,   27.6462,  29.115),
                ( 121.5130,   10.3531,   0.000),
                ( 124.0000,   10.3531,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  69.6203,   44.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 12.0 mm thick, 8 cut(s)
        axis=2, lo=4.0000, hi=16.0000, ring=None,
        cuts=[
            ("poly", [
                (  88.1008,    4.8908,   0.250),
                ( 113.0145,    5.3324,   1.500),
                ( 106.4831,   14.4631,   1.500),
            ]),
            ("poly", [
                (  82.9260,    9.7750,   1.500),
                (  98.8307,   17.0805,   1.500),
                (  94.1439,   21.9893,   1.500),
                (  77.5797,    8.4604,   0.250),
            ]),
            ("poly", [
                (  84.2180,   23.4219,   0.250),
                (  76.5166,   32.7748,   0.250),
                (  81.3277,   22.1752,   1.500),
            ]),
            ("poly", [
                (  77.6969,   20.5345,   0.250),
                (  75.0563,   19.0159,   1.500),
                (  75.2857,   16.3032,   0.250),
                (  78.3987,   17.7660,   1.500),
            ]),
            ("poly", [
                (  37.6783,   14.2558,   1.500),
                (  47.0878,   19.4922,   0.250),
                (  44.6156,   22.1663,   1.500),
                (  39.4802,   24.2782,   1.500),
                (  36.4482,   21.1650,   1.500),
            ]),
            ("poly", [
                (  33.0275,   24.4447,   0.250),
                (  33.9455,   23.3025,   0.250),
                (  34.4393,   23.9690,   0.250),
            ]),
            ("poly", [
                (  -1.5010,   13.2780,   0.000),
                (  18.4391,   19.9327,   1.500),
                (  -1.5113,   26.8106,   0.000),
                (  -4.0000,   26.8106,   0.000),
                (  -4.0000,   13.2780,   0.000),
            ]),
            ("poly", [
                (  69.5475,   40.0600,   0.000),
                (  89.5475,   40.0000,  22.285),
                ( 114.9460,   24.8175,   0.000),
                ( 113.9205,   22.9863,   1.500),
                ( 121.5076,    9.4275,   0.000),
                ( 124.0000,    9.4275,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  69.5475,   44.0000,   0.000),
            ]),
        ],
    ),
]

part = None
for p in PLATES:
    b = plate(p["axis"], p["lo"], p["hi"], p["cuts"], BOX, p.get("ring"))
    part = b if part is None else part + b

# ---- no re-entrant (concave) edges: fillet them ---------------------------------
R_FILLET = 0.00    # fillet radius on re-entrant edges (mm)
part, n_fillet, n_skip = fillet_reentrant(part.clean(), R_FILLET)
print("re-entrant edges filleted: %d%s" % (n_fillet, ", skipped: %d" % n_skip if n_skip else ""))

part = part.clean()

if __name__ == "__main__":
    part.exportStep(sys.argv[1] if len(sys.argv) > 1 else "Cantilever_exact.step")
