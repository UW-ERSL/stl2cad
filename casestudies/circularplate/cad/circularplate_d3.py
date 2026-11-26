"""circularplate_d3: topology-optimized part as plates of the original part.

Original features first, then 1 plate(s); units mm.
Run:  python circularplate_d3.py [out.step]
"""

import math
import sys

from cadquery import Edge, Solid, Vector, Vertex
from OCP.BRepOffset import BRepOffset_Analyse
from OCP.ChFiDS import ChFiDS_TypeOfConcavity
from cadquery import func as cf

# ---- parameters ----------------------------------------------------------
z0 = 0.0    # level along z
z1 = 12.7    # level along z
hole_d = 25.4

DIR = (0, 0, 1)


def P(a, b):
    """Point (x, y) in the sketch plane (z = 0)."""
    return (a, b, 0)


def prism(outline, lo, hi):
    """Outline extruded from level lo to level hi."""
    body = cf.extrude(cf.face(outline), [(hi - lo) * c for c in DIR])
    return body.translate(Vector(*[lo * c for c in DIR]))


def hole(a, b, d, lo, hi):
    return prism(Edge.makeCircle(d / 2, P(a, b), DIR), lo, hi)


# ---- 1. extrusion: z 0.0 to 12.7 ------------------------------
outline_1 = cf.wire(
    Edge.makeCircle(63.5, P(63.5, 63.5), DIR),
)
part = prism(outline_1, z0, z1)

# ---- round holes: (x, y, diameter, from, to) ------------------------
for a, b, d, lo, hi in [
    (63.5, 63.5, hole_d, z0, z1),
]:
    part = part - hole(a, b, d, lo, hi)


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

BOX = ((0.0, 0.0815, 0.0), (127.0, 126.9185, 12.7))

PLATES = [
    dict(  # normal to z, 12.7 mm thick, 10 cut(s)
        axis=2, lo=0.0000, hi=12.7000, ring=None,
        cuts=[
            ("poly", [
                (  34.1303,   39.7368,   2.500),
                (   0.7133,   55.7541,   2.500),
                (  17.6539,   22.0113,  44.205),
                (  43.9717,    5.0244,   0.750),
            ]),
            ("poly", [
                (  65.2618,   98.6066,   2.500),
                (  41.9130,  118.0695,   2.500),
                (  38.1276,   88.2300,   2.500),
                (  54.4245,   87.0859,   2.500),
            ]),
            ("poly", [
                (  77.0244,   28.2245,   2.500),
                (  51.4646,    1.3926,   2.500),
                (  88.7909,    7.0770,  44.205),
                ( 113.0790,   26.8575,   0.750),
            ]),
            ("poly", [
                ( 101.2282,   65.4617,   2.500),
                ( 118.8485,   32.8614,   2.500),
                ( 124.9768,   70.1174,  44.205),
                ( 113.6699,   99.3293,   0.750),
            ]),
            ("poly", [
                (  73.2930,   99.9879,   2.500),
                ( 109.7426,  106.6717,   2.500),
                (  76.2039,  124.0128,  44.205),
                (  44.9277,  122.2863,   0.750),
            ]),
            ("poly", [
                (  31.8242,   84.0890,   2.500),
                (  36.7310,  120.8201,   2.500),
                (   9.8746,   94.2816,  44.205),
                (   1.8518,   64.0026,   0.750),
            ]),
            ("poly", [
                (  30.6561,   76.0241,   2.500),
                (   4.9306,   59.8324,   2.500),
                (  32.1399,   47.0114,   2.500),
                (  38.2640,   62.1571,   2.500),
            ]),
            ("poly", [
                (  41.4395,   36.1337,   2.500),
                (  48.8892,    6.6638,   2.500),
                (  69.4908,   28.5795,   2.500),
                (  56.9788,   39.0841,   2.500),
            ]),
            ("poly", [
                (  82.7098,   34.0626,   2.500),
                ( 113.0394,   32.0409,   2.500),
                (  98.5626,   58.4066,   2.500),
                (  84.7057,   49.7530,   2.500),
            ]),
            ("poly", [
                (  97.4328,   72.6729,   2.500),
                ( 108.7279,  100.8934,   2.500),
                (  79.1791,   95.2726,   2.500),
                (  83.1271,   79.4198,   2.500),
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
    part.exportStep(sys.argv[1] if len(sys.argv) > 1 else "CircularPlateHole_exact.step")
