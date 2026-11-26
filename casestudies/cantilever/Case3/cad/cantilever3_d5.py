"""cantilever3_d5: topology-optimized part as plates of the original part.

Original features first, then 1 plate(s); units mm.
Run:  python cantilever3_d5.py [out.step]
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
    dict(  # normal to z, 20.0 mm thick, 7 cut(s)
        axis=2, lo=0.0000, hi=20.0000, ring=None,
        cuts=[
            ("poly", [
                (  69.5099,    1.6152,   1.152),
                ( 115.3092,    1.6876,   1.001),
                ( 116.8868,    3.1608,   1.000),
                ( 101.3658,   26.5258,   1.103),
                (  99.1135,   28.2598,   1.000),
                (  93.7398,   24.7786,   1.848),
                (  79.1480,   12.1995,  16.344),
                (  66.4581,    2.5842,   0.426),
            ]),
            ("poly", [
                (  93.7647,   29.4262,   0.622),
                (  28.2740,   32.8145,   0.362),
                (  64.4614,   13.6994,   4.324),
                (  67.9579,    7.3243,   1.947),
            ]),
            ("poly", [
                (  62.2650,    7.1802,   2.304),
                (  64.2243,    8.2713,   0.621),
                (  62.6696,    9.7078,   1.472),
                (  60.4118,    9.2978,   2.549),
                (  58.0996,    6.7537,   0.437),
            ]),
            ("poly", [
                (  -1.2967,   -0.7982,   0.000),
                (  -0.0495,    3.7099,   1.825),
                (   2.1547,    4.7560,   5.483),
                (  57.0895,    9.5201,   5.985),
                (  60.7084,   11.6444,   0.947),
                (  26.6514,   29.6570,  52.290),
                (   9.2207,   34.6547,   2.813),
                (   7.2727,   34.0925,   1.685),
                (   3.1735,   30.7562,   3.097),
                (   0.4357,   30.7682,   1.177),
                (  -1.5025,   34.3670,   0.000),
                (  -4.0000,   34.3670,   0.000),
                (  -4.0000,   -0.7982,   0.000),
            ]),
            ("poly", [
                (  -1.5180,   28.5043,   0.000),
                (   0.2472,   32.6421,   1.162),
                (   2.6663,   34.4354,   0.886),
                (  -0.2263,   37.1216,   2.319),
                (  -0.4467,   41.4341,   0.000),
                (  -0.4467,   44.0000,   0.000),
                (  -4.0000,   44.0000,   0.000),
                (  -4.0000,   28.5043,   0.000),
            ]),
            ("poly", [
                (   0.0000,   40.0578,   0.000),
                (  19.2712,   40.0000,   1.007),
                (  25.4308,   37.3145,  15.042),
                (  97.2864,   33.4897,   2.390),
                ( 100.6683,   32.0113,   5.147),
                ( 121.5037,    1.1249,   0.000),
                ( 124.0000,    1.1249,   0.000),
                ( 124.0000,   44.0000,   0.000),
                (   0.0000,   44.0000,   0.000),
            ]),
            ("arc", (55.3565, -15.1807, 17.0569, 1.096913, 2.044880), [
                (  47.5696,   -1.5000),
                (  47.5696,   -4.0000),
                (  63.1403,   -4.0000),
                (  63.1403,   -1.5000),
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
