"""cantilever6_d4: topology-optimized part as plates of the original part.

Original features first, then 4 plate(s); units mm.
Run:  python cantilever6_d4.py [out.step]
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
                (  60.0000,   40.0349,   0.000),
                (  71.6243,   40.0000,  41.586),
                ( 104.6089,   27.1394,  37.674),
                ( 121.5050,    3.0284,   0.000),
                ( 124.0000,    3.0284,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
            ("poly", [
                (  -1.5484,    7.5122,   0.000),
                (  21.7110,    8.9725,   0.000),
                (  21.7380,    8.2019,   0.250),
                (  89.3974,    1.6067,   0.000),
                (  89.4585,    2.2916,   0.250),
                ( 116.1463,    2.6981,   0.250),
                ( 104.1628,   22.5869,  27.376),
                (  80.2714,   34.2273,   0.000),
                (  79.7250,   33.4146,   1.250),
                (  62.3056,   35.4417,  64.107),
                (  -1.5469,   29.9699,   0.000),
                (  -4.0000,   29.9699,   0.000),
                (  -4.0000,    7.5122,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 8.0 mm thick, 3 cut(s)
        axis=2, lo=6.0000, hi=14.0000, ring=None,
        cuts=[
            ("poly", [
                (  81.8193,    2.3300,   0.250),
                ( 116.1975,    2.8638,   1.250),
                ( 105.5301,   17.3019,   1.250),
            ]),
            ("poly", [
                (  60.0000,   40.0369,   0.000),
                (  72.2873,   40.0000,  45.588),
                ( 105.6076,   27.0529,  33.799),
                ( 121.5227,    2.4116,   0.000),
                ( 124.0000,    2.4116,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
            ("poly", [
                (  -1.5108,    2.8784,   0.000),
                (   0.5349,    6.8637,   1.250),
                (  17.3398,   13.4061,  27.498),
                (  77.0294,    4.5101,  29.880),
                ( 101.4450,   17.4832,   1.250),
                ( 100.7188,   22.5235,   1.250),
                (  73.6152,   34.2326,  44.426),
                (  16.2549,   27.4297,   1.250),
                (  -1.5428,   33.5250,   0.000),
                (  -4.0000,   33.5250,   0.000),
                (  -4.0000,    2.8784,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 5.0 mm thick, 4 cut(s)
        axis=2, lo=7.5000, hi=12.5000, ring=None,
        cuts=[
            ("poly", [
                (  82.2608,    2.2334,   0.250),
                ( 116.1793,    2.9184,   1.250),
                ( 106.1434,   16.2385,   1.250),
            ]),
            ("poly", [
                (  84.2254,   29.1970,  36.843),
                (  66.2364,   32.0798,   1.250),
                (  62.3711,   29.1745,   1.250),
                (  54.2521,    9.0469,   1.250),
                (  58.4950,    8.2938,   1.250),
                (  68.2767,   10.2839,   1.250),
                (  72.4213,    7.2826,   1.250),
                (  83.9475,    9.0477,   1.250),
                ( 100.0618,   19.0601,   0.250),
                (  98.7551,   20.8068,   1.250),
            ]),
            ("poly", [
                (  -1.5405,    4.8758,   0.000),
                (  19.6662,   16.8682,   1.250),
                (  20.9846,   22.3024,   1.250),
                (  -1.5271,   34.9644,   0.000),
                (  -4.0000,   34.9644,   0.000),
                (  -4.0000,    4.8758,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0365,   0.000),
                (  72.1608,   40.0000,  43.963),
                ( 105.5916,   27.1548,  33.521),
                ( 121.5016,    2.4389,   0.000),
                ( 124.0000,    2.4389,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (  60.0000,   44.0000,   0.000),
            ]),
        ],
    ),
    dict(  # normal to z, 3.5 mm thick, 6 cut(s)
        axis=2, lo=8.2500, hi=11.7500, ring=None,
        cuts=[
            ("poly", [
                (  82.6467,    2.1846,   0.250),
                ( 116.2492,    2.9641,   1.250),
                ( 106.2364,   15.9162,   1.250),
            ]),
            ("poly", [
                (  92.8622,   20.7369,   1.250),
                (  87.0918,   12.7804,   0.250),
                (  94.7300,   16.9339,   1.250),
                (  95.9733,   20.5403,   0.250),
            ]),
            ("poly", [
                (  88.6333,   21.6093,   1.250),
                (  90.4156,   25.7269,   0.250),
                (  85.7538,   24.3354,   1.250),
                (  79.5420,   16.0202,   0.250),
                (  85.5533,   18.2769,   1.250),
            ]),
            ("poly", [
                (  81.8120,   24.4349,   0.250),
                (  82.1529,   26.2775,   0.250),
                (  79.0874,   26.4593,   0.250),
                (  79.5271,   24.8430,   1.250),
            ]),
            ("poly", [
                (  -1.5169,    4.5795,   0.000),
                (  22.1618,   19.5740,   1.250),
                (  -1.5282,   35.6641,   0.000),
                (  -4.0000,   35.6641,   0.000),
                (  -4.0000,    4.5795,   0.000),
            ]),
            ("poly", [
                (  60.0000,   40.0358,   0.000),
                (  71.9358,   40.0000,  41.927),
                ( 105.1134,   27.5834,  34.248),
                ( 121.5230,    2.7126,   0.000),
                ( 124.0000,    2.7126,   0.000),
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
