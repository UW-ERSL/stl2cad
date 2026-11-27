"""cantilever2_d4: topology-optimized part as plates of the original part.

Original features first, then 4 plate(s); units mm.
Run:  python cantilever2_d4.py [out.step]
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
    dict(  # normal to z, 20.0 mm thick, 2 cut(s)
        axis=2, lo=0.0000, hi=20.0000, ring=None,
        cuts=[
            ("poly", [
                (  -0.8210,   -1.3028,   0.000),
                (  -0.0031,    8.8326,   1.250),
                (  23.3446,   10.3261,   0.000),
                (  23.3933,    9.7318,   0.250),
                (  80.6717,    2.5653,   0.000),
                (  80.6814,    2.6798,   0.000),
                ( 113.5923,    2.0550,   0.000),
                ( 113.6135,    1.3253,   0.250),
                ( 115.5660,    2.0607,   1.250),
                ( 115.8132,    4.3209,   1.250),
                ( 104.9288,   21.8430,  22.237),
                (  92.5603,   29.8591,  37.967),
                (  78.7020,   34.0475,  41.685),
                (  61.8527,   34.8329,  82.477),
                (  18.5539,   29.6701,  68.796),
                (   0.3069,   31.4333,   1.250),
                (  -1.5097,   37.1602,   0.000),
                (  -4.0000,   37.1602,   0.000),
                (  -4.0000,   -4.0000,   0.000),
                (  -0.8210,   -4.0000,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0391,   0.000),
                (  73.0205,   40.0000,   1.250),
                (  96.7548,   31.5668,   1.250),
                (  98.2532,   32.5198,   0.000),
                ( 121.5331,    6.5256,   0.000),
                ( 124.0000,    6.5256,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 16.5 mm thick, 2 cut(s)
        axis=2, lo=1.7500, hi=18.2500, ring=None,
        cuts=[
            ("poly", [
                (  60.0000,   40.0473,   0.000),
                (  75.7764,   40.0000,  35.695),
                ( 105.5874,   27.7821,  35.567),
                ( 121.5018,    4.8504,   0.000),
                ( 124.0000,    4.8504,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
            ("poly", [
                (  -1.5418,    8.7148,   0.000),
                (  17.2468,   12.7373,   1.250),
                (  88.9107,    3.0495,   1.250),
                ( 116.5511,    3.2594,   0.250),
                (  96.4592,   27.8395,  33.421),
                (  64.1470,   34.6834,   0.000),
                (  63.7399,   33.2365,   1.250),
                (  11.2961,   27.0995,   0.000),
                (  11.1206,   28.5164,   1.250),
                (  -1.5038,   31.3024,   0.000),
                (  -4.0000,   31.3024,   0.000),
                (  -4.0000,    8.7148,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 8.0 mm thick, 5 cut(s)
        axis=2, lo=6.0000, hi=14.0000, ring=None,
        cuts=[
            ("poly", [
                (  84.2935,    2.8065,   0.250),
                ( 115.3408,    3.9847,   1.250),
                ( 105.6850,   16.7917,   3.000),
            ]),
            ("poly", [
                (  54.2336,    9.4596,   3.000),
                (  83.3078,    7.7655,  15.609),
                ( 100.8508,   17.6399,   1.250),
                ( 100.4277,   21.5874,   1.250),
                (  80.3239,   31.4983,  31.873),
                (  54.7084,   29.7157,   3.000),
                (  51.7328,   25.2525,   1.250),
            ]),
            ("poly", [
                (  34.1517,   13.5403,   1.250),
                (  46.6035,   16.8621,   3.000),
                (  50.4671,   24.7099,   3.000),
                (  30.1423,   26.7296,   3.000),
                (  29.6658,   14.7051,   1.250),
            ]),
            ("poly", [
                (  -1.5456,    7.6122,   0.000),
                (  20.1505,   16.8337,   1.250),
                (  22.0333,   22.2027,   3.000),
                (  -1.5143,   32.5029,   0.000),
                (  -4.0000,   32.5029,   0.000),
                (  -4.0000,    7.6122,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0475,   0.000),
                (  75.8497,   40.0000,  38.878),
                ( 106.2653,   27.3380,  32.952),
                ( 121.5277,    4.6270,   0.000),
                ( 124.0000,    4.6270,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 5.5 mm thick, 5 cut(s)
        axis=2, lo=7.2500, hi=12.7500, ring=None,
        cuts=[
            ("poly", [
                (  85.2753,    2.8824,   0.250),
                ( 115.1065,    4.2294,   1.250),
                ( 106.3514,   16.0043,   3.000),
            ]),
            ("poly", [
                (  92.6382,   21.7350,   1.250),
                (  83.5240,   13.2708,   1.250),
                (  82.9737,   11.0663,   0.250),
                (  92.1696,   14.6969,   3.000),
                (  96.4051,   17.8361,   1.250),
                (  96.9777,   20.3980,   1.250),
            ]),
            ("poly", [
                (  82.7427,   27.0722,  12.887),
                (  77.0520,   27.4419,   1.250),
                (  77.7679,   17.3181,   0.250),
                (  87.8986,   19.9283,   3.000),
                (  90.4660,   23.3519,   1.250),
                (  87.8440,   25.7003,   1.250),
            ]),
            ("poly", [
                (  -1.5190,    7.1762,   0.000),
                (  19.4316,   18.1728,   3.000),
                (  15.4798,   24.4352,  10.215),
                (  -1.5262,   32.3741,   0.000),
                (  -4.0000,   32.3741,   0.000),
                (  -4.0000,    7.1762,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0460,   0.000),
                (  75.3297,   40.0000,  34.147),
                ( 105.9585,   27.8899,  33.992),
                ( 121.5041,    4.7161,   0.000),
                ( 124.0000,    4.7161,   0.000),
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
