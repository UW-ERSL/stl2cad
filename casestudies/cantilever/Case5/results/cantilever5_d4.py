"""cantilever5_d4: topology-optimized part as plates of the original part.

Original features first, then 5 plate(s); units mm.
Run:  python cantilever5_d4.py [out.step]
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
    dict(  # normal to z, 20.0 mm thick, 1 cut(s)
        axis=2, lo=0.0000, hi=20.0000, ring=None,
        cuts=[
            ("poly", [
                (  -1.5422,    5.2635,   0.000),
                (  13.1748,    6.5182,  82.961),
                (  92.4443,    0.2782,   0.070),
                (  90.7308,   -0.0000,   1.598),
                (  70.7308,   -0.0600,   0.000),
                (  70.7308,   -4.0000,   0.000),
                ( 124.0000,   -4.0000,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  -4.0000,   44.0000,   0.000),
                (  -4.0000,    5.2635,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 20.0 mm thick, 1 cut(s)
        axis=2, lo=0.0000, hi=20.0000, ring=None,
        cuts=[
            ("poly", [
                (  60.0000,   40.0200,   0.000),
                (  63.8617,   40.0000,   1.000),
                (  71.1034,   37.3189,   1.000),
                (   0.7600,   33.0425,   1.000),
                (  -1.5141,   39.5756,   0.000),
                (  -4.0000,   39.5756,   0.000),
                (  -4.0000,   -4.0000,   0.000),
                ( 124.0000,   -4.0000,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 14.0 mm thick, 2 cut(s)
        axis=2, lo=3.0000, hi=17.0000, ring=None,
        cuts=[
            ("poly", [
                (  60.0000,   40.0254,   0.000),
                (  68.4826,   40.0000,  22.453),
                ( 102.2565,   28.4529,  38.532),
                ( 121.5182,    3.1256,   0.000),
                ( 124.0000,    3.1256,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
            ("poly", [
                (  -1.5197,    1.6664,   0.000),
                (   2.2954,    7.5191,   1.000),
                ( 118.0906,    0.5230,   1.000),
                ( 100.0944,   23.8977,  27.818),
                (  70.4488,   34.9710,   0.000),
                (  69.9618,   36.2360,   0.000),
                (   9.5923,   32.9318,   0.000),
                (   9.6103,   32.4961,   0.000),
                (  -1.5071,   34.7419,   0.000),
                (  -4.0000,   34.7419,   0.000),
                (  -4.0000,    1.6664,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 5.5 mm thick, 4 cut(s)
        axis=2, lo=7.2500, hi=12.7500, ring=None,
        cuts=[
            ("poly", [
                ( 105.3794,   16.6675,   1.000),
                (  79.2107,    2.3957,   1.000),
                ( 116.3647,    2.9278,   1.000),
            ]),
            ("poly", [
                (  30.1704,   25.2444,   1.000),
                (  28.7584,   17.0153,   1.000),
                (  22.8131,   10.5545,   1.000),
                (  28.3692,    8.7672,   1.000),
                (  38.2249,    9.5500,   1.000),
                (  40.7509,    6.8529,   1.000),
                (  43.6786,    6.2693,   1.000),
                (  78.0776,    6.2858,   1.000),
                ( 100.1476,   18.6943,   1.000),
                (  98.2696,   21.5071,   1.000),
                (  92.6646,   23.0434,   1.000),
                (  89.0584,   25.9674,   1.000),
                (  65.0535,   34.0621,   1.000),
                (  26.5450,   31.8573,   1.000),
                (  24.8407,   28.2400,   1.000),
            ]),
            ("poly", [
                (  -1.5355,    3.7894,   0.000),
                (  23.3910,   16.6043,   1.000),
                (  24.5147,   22.1383,   1.000),
                (  -1.5030,   36.5835,   0.000),
                (  -4.0000,   36.5835,   0.000),
                (  -4.0000,    3.7894,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0275,   0.000),
                (  69.1605,   40.0000,  44.412),
                ( 103.8914,   27.3722,  35.003),
                ( 121.5255,    2.4223,   0.000),
                ( 124.0000,    2.4223,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 3.5 mm thick, 9 cut(s)
        axis=2, lo=8.2500, hi=11.7500, ring=None,
        cuts=[
            ("poly", [
                (  79.9940,    2.3628,   1.000),
                ( 116.4499,    2.9079,   1.000),
                ( 105.6870,   16.2505,   1.000),
            ]),
            ("poly", [
                (  81.1140,   11.9715,   1.000),
                (  81.9705,    8.6320,   1.000),
                (  98.8248,   17.9400,   1.000),
                (  95.8649,   21.7582,   1.000),
            ]),
            ("poly", [
                (  84.2909,   18.0920,   0.000),
                (  84.6258,   17.8536,   0.000),
                (  90.9603,   22.1681,   1.000),
                (  87.8286,   26.9299,   1.000),
                (  76.9449,   14.8795,   1.000),
                (  79.1434,   12.0699,   1.000),
            ]),
            ("poly", [
                (  71.9122,   26.4930,   1.000),
                (  78.4016,   20.0160,   1.000),
                (  83.7130,   27.5030,   1.000),
                (  76.6915,   30.4363,   1.000),
            ]),
            ("poly", [
                (  71.6481,   12.5441,   1.000),
                (  73.3599,    8.9556,   1.000),
                (  76.5722,   11.6045,   1.000),
                (  74.5149,   14.5488,   1.000),
            ]),
            ("poly", [
                (  61.0862,   31.5615,   1.000),
                (  60.7495,   29.5427,   1.000),
                (  62.6320,   23.0533,   1.000),
                (  63.8273,   21.4425,   1.000),
                (  66.5445,   19.2608,   1.000),
                (  68.0836,   13.7739,   1.000),
                (  73.9632,   18.7445,   1.000),
                (  76.0033,   21.9932,   1.000),
                (  72.3540,   22.1871,   1.000),
                (  69.5371,   20.9393,   1.000),
                (  69.5580,   26.4033,   1.000),
                (  74.1778,   29.6030,   1.000),
                (  71.6719,   31.6066,   1.000),
                (  64.8260,   32.9263,   1.000),
                (  62.8382,   32.7147,   1.000),
            ]),
            ("poly", [
                (  64.0175,   11.1576,   1.000),
                (  61.1214,   12.0047,   1.000),
                (  58.9872,    9.4460,   1.000),
                (  62.2988,    8.4550,   1.000),
            ]),
            ("poly", [
                (  -1.5215,    1.6030,   0.000),
                (   0.4136,    5.3132,   1.000),
                (  10.7391,    9.9150,  19.375),
                (  20.2529,   19.1950,   1.000),
                (  18.2107,   24.2036,   1.000),
                (  -1.5175,   37.0670,   0.000),
                (  -4.0000,   37.0670,   0.000),
                (  -4.0000,    1.6030,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0234,   0.000),
                (  67.8018,   40.0000,   1.000),
                ( 101.4208,   29.0177,  37.639),
                ( 121.5291,    3.7075,   0.000),
                ( 124.0000,    3.7075,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
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
