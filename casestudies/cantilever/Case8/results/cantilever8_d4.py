"""cantilever8_d4: topology-optimized part as plates of the original part.

Original features first, then 6 plate(s); units mm.
Run:  python cantilever8_d4.py [out.step]
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
                (  -1.4861,   -0.2826,   0.000),
                (   1.1027,    3.1557,   1.000),
                (  53.6696,    1.1754,   0.250),
                (  54.0718,    0.2811,   0.250),
                ( 141.9002,   -0.0000,   0.000),
                ( 121.9002,   -0.0600,   0.000),
                ( 124.0000,   -0.0600,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  -4.0000,   44.0000,   0.000),
                (  -4.0000,   -0.2826,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 20.0 mm thick, 1 cut(s)
        axis=2, lo=0.0000, hi=20.0000, ring=None,
        cuts=[
            ("arc", (17.4090, 186.4568, 149.5325, -1.360133, -1.687879), [
                (  -1.5000,   37.9480),
                (  -4.0000,   37.9480),
                (  -4.0000,   -4.0000),
                ( 124.0000,   -4.0000),
                ( 124.0000,   44.0000),
                (  48.6776,   44.0000),
                (  48.6776,   41.5000),
            ]),
        ],
    ),
    dict(  # normal to z, 10.5 mm thick, 2 cut(s)
        axis=2, lo=4.7500, hi=15.2500, ring=None,
        cuts=[
            ("poly", [
                (  60.0000,   40.0200,   0.000),
                (  63.1332,   40.0000,   1.000),
                ( 105.6557,   22.0711,   1.000),
                ( 121.5126,    2.3140,   0.000),
                ( 124.0000,    2.3140,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
            ("poly", [
                ( 112.1074,   -1.5033,   0.000),
                ( 116.7607,    3.0946,   1.000),
                ( 105.0836,   20.1825,  13.628),
                (  63.6947,   36.1721,   1.000),
                (  63.5983,   37.3864,   0.000),
                (  -1.5358,   36.0591,   0.000),
                (  -4.0000,   36.0591,   0.000),
                (  -4.0000,   -4.0000,   0.000),
                ( 112.1074,   -4.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 10.5 mm thick, 1 cut(s)
        axis=2, lo=4.7500, hi=15.2500, ring=None,
        cuts=[
            ("poly", [
                (   0.5033,   -1.5297,   0.000),
                (  -0.7462,    3.0052,   1.000),
                (   7.8912,    4.1158,  19.682),
                (  29.1150,    2.0547,   0.000),
                (  29.1406,    2.0994,   0.000),
                (  70.3762,    3.6843,   0.000),
                (  70.4030,    3.3714,   0.000),
                (  76.9273,    2.4006,   1.000),
                (  81.1106,    0.0088,   0.250),
                (  60.0000,   -0.0546,   0.000),
                (  60.0000,   -4.0000,   0.000),
                ( 124.0000,   -4.0000,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  -4.0000,   44.0000,   0.000),
                (  -4.0000,   -4.0000,   0.000),
                (   0.5033,   -4.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 5.5 mm thick, 3 cut(s)
        axis=2, lo=7.2500, hi=12.7500, ring=None,
        cuts=[
            ("poly", [
                ( 104.6671,   18.9580,   1.000),
                (  79.1262,    2.3124,   1.000),
                ( 116.6943,    2.6897,   0.250),
            ]),
            ("poly", [
                (  -1.0755,   -1.0709,   0.000),
                (   0.2598,    2.6615,   1.000),
                (  13.7798,    7.3657,   1.000),
                (  15.8126,    5.2921,   1.000),
                (  36.9308,    1.5499,   0.000),
                (  37.0551,    2.1514,   0.250),
                (  77.4651,    5.0182,  11.424),
                ( 101.2036,   19.9970,   0.250),
                (  64.4853,   35.6005,   1.000),
                (  60.6684,   32.4488,   1.000),
                (  56.3526,   23.3113,   0.250),
                (  53.0628,   27.0620,   0.250),
                (  58.6936,   30.7603,   1.000),
                (  61.3323,   36.6035,   1.000),
                (   9.4261,   35.1715,   0.250),
                (  18.6211,   27.5332,   0.250),
                (  17.7216,   26.2633,   0.250),
                (  -1.5157,   38.8938,   0.000),
                (  -4.0000,   38.8938,   0.000),
                (  -4.0000,   -1.0709,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0200,   0.000),
                (  64.2103,   40.0000,   1.000),
                ( 106.1670,   23.0300,   1.000),
                ( 121.5162,    2.4489,   0.000),
                ( 124.0000,    2.4489,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 3.0 mm thick, 8 cut(s)
        axis=2, lo=8.5000, hi=11.5000, ring=None,
        cuts=[
            ("poly", [
                (  80.3841,    2.9166,   0.250),
                ( 116.5211,    2.8121,   0.250),
                ( 105.1089,   18.1858,   1.000),
            ]),
            ("poly", [
                ( 100.0318,   20.2415,   1.000),
                (  64.2195,   35.5512,   1.000),
                (  59.1996,   27.0404,   1.000),
                (  76.9079,    6.1339,   1.000),
                (  81.3988,    7.2754,   0.000),
                (  81.1831,    7.7913,   0.250),
            ]),
            ("poly", [
                (  56.4840,   21.9411,   1.000),
                (  42.7937,    8.9847,   1.000),
                (  34.4366,    2.8414,   0.250),
                (  71.1406,    5.2621,   1.000),
                (  70.5131,    8.0152,   1.000),
                (  58.6817,   21.7826,   1.000),
            ]),
            ("poly", [
                (  60.2841,   35.8219,   1.000),
                (  47.5406,   36.3852,   0.250),
                (  55.5885,   30.5086,   1.000),
            ]),
            ("poly", [
                (  38.6660,   13.8925,   1.000),
                (  51.8649,   26.7529,   1.000),
                (  43.0426,   30.8080,   1.000),
                (  32.8195,   13.1029,   1.000),
            ]),
            ("poly", [
                (  21.1741,   29.0030,   0.000),
                (  21.0312,   28.7711,   0.000),
                (  31.2952,   21.1278,   1.000),
                (  33.9928,   21.6441,   1.000),
                (  40.4017,   32.4911,   1.000),
                (  13.0965,   37.2463,   1.000),
                (   8.4452,   35.9002,   0.250),
            ]),
            ("poly", [
                (  -1.5331,    0.2945,   0.000),
                (  24.9812,   16.8003,   1.000),
                (  -1.5229,   39.9054,   0.000),
                (  -4.0000,   39.9054,   0.000),
                (  -4.0000,    0.2945,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0200,   0.000),
                (  64.6306,   40.0000,   1.000),
                ( 106.5025,   23.0507,  12.819),
                ( 121.5120,    1.9084,   0.000),
                ( 124.0000,    1.9084,   0.000),
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
