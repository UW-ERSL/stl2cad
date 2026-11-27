"""cantilever7_d4: planar-truss model of a topology-optimized part (stl2cad, to/swtruss.py); units mm.

Parameters: joints, member width classes, corner radius R, bar depth t; edit and run again:
    python cantilever7_d4.py [out.step]
Members are bars between their joints; chords are strips flush with the outline (edge: a point on the
outline and the inward normal; ext: run-on at free ends). Joints are mitred, concave corners rounded with R.
The section is extruded t on both faces of the domain (depth H); webs are boxes between the faces.
"""
import math
import sys

import cadquery as cq

# ---- parameters ------------------------------------------------------------
JOINTS = {
    'J1': (0.7500, 1.5827),
    'J2': (2.0000, 37.3410),
    'J3': (5.5256, 37.1913),
    'J4': (36.5000, 37.2924),
    'J5': (37.5000, 1.6132),
    'J6': (53.4698, 2.2013),
    'J7': (79.2500, 1.5827),
    'J8': (83.4716, 36.6064),
    'J9': (119.9438, 1.7040),
}
WIDTHS = {
    'w1': 3.2286,
    'w2': 4.2525,
    'w3': 5.4243,
}
MEMBERS = [
    {'a': 'J8', 'c': 'J9', 'w': 'w2'},
    {'a': 'J1', 'c': 'J5', 'w': 'w1', 'edge': ((19.25, 0.0), (0.0, 1.0)), 'ext': (1.25, 0.0)},
    {'a': 'J5', 'c': 'J6', 'w': 'w1', 'edge': ((19.25, 0.0), (0.0, 1.0)), 'ext': (0.0, 0.0)},
    {'a': 'J2', 'c': 'J3', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (2.5, 0.0)},
    {'a': 'J3', 'c': 'J4', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (0.0, 0.0)},
    {'a': 'J4', 'c': 'J8', 'w': 'w3', 'edge': ((6.75, 40.0), (0.0, -1.0)), 'ext': (0.0, 0.0)},
    {'a': 'J3', 'c': 'J6', 'w': 'w2'},
    {'a': 'J9', 'c': 'J7', 'w': 'w1', 'edge': ((97.75, 0.0), (0.0, 1.0)), 'ext': (0.5562, 0.0)},
    {'a': 'J6', 'c': 'J8', 'w': 'w2'},
]
R = 1.3          # corner radius
t = 3.2320        # bar depth (each face)
BOX = (0.0, 0.0, 119.99999731779099, 39.99999910593033)    # domain section (x0, y0, x1, y1)
H = 20.0000        # domain depth
WEBS = [(0.0, 0.0, 91.875, 1.7886), (116.875, 0.0, 120.0, 2.5993)]    # (x0, y0, x1, y1), between the faces

# ---- geometry --------------------------------------------------------------
MITRE = 3.0          # a mitre point farther than this x the widest bar from its joint is replaced by a bevel


def _unit(v):
    n = math.hypot(v[0], v[1])
    return (v[0] / n, v[1] / n)


def _isect(p, d, q, e):
    den = d[0] * e[1] - d[1] * e[0]
    if abs(den) < 1e-9:
        return None
    s = ((q[0] - p[0]) * e[1] - (q[1] - p[1]) * e[0]) / den
    return (p[0] + s * d[0], p[1] + s * d[1])


def bars(joints, widths, members):
    """Each member as (end at a, end at c, width): chords on their centre line (foot on the outline + n w/2),
    run on by their ext at free ends."""
    out = []
    for m in members:
        A, C, w = joints[m["a"]], joints[m["c"]], widths[m["w"]]
        if m.get("edge") is not None:
            (fx, fy), (nx, ny) = m["edge"]
            t = (ny, -nx)
            sa = (A[0] - fx) * t[0] + (A[1] - fy) * t[1]
            sc = (C[0] - fx) * t[0] + (C[1] - fy) * t[1]
            u = 1.0 if sc >= sa else -1.0
            ea, ec = m.get("ext", (0.0, 0.0))
            sa, sc = sa - u * ea, sc + u * ec
            A = (fx + sa * t[0] + nx * w / 2, fy + sa * t[1] + ny * w / 2)
            C = (fx + sc * t[0] + nx * w / 2, fy + sc * t[1] + ny * w / 2)
        out.append((A, C, w))
    return out


def chord_mitres(joints, widths, members, B):
    """A chord end at a joint where another member carries on past it (no collinear chord) ends where its outer
    edge (the outline) meets that member's outer edge: the outer corner is a mitre, with no step."""
    B = [list(b) for b in B]
    for i, m in enumerate(members):
        if m.get("edge") is None:
            continue
        (fx, fy), (nx, ny) = m["edge"]
        for end in (0, 1):
            q = m["a"] if end == 0 else m["c"]
            E, O = (B[i][0], B[i][1]) if end == 0 else (B[i][1], B[i][0])
            d = _unit((O[0] - E[0], O[1] - E[1]))                       # into the chord from this end
            others = [k for k, mk in enumerate(members) if k != i and q in (mk["a"], mk["c"])]
            if not others or any(members[k].get("edge") is not None and members[k]["edge"] == m["edge"] for k in others):
                continue
            best = None
            for k in others:
                A2, C2, w2 = B[k]
                u = _unit((C2[0] - A2[0], C2[1] - A2[1])) if members[k]["a"] == q else _unit((A2[0] - C2[0], A2[1] - C2[1]))
                if u[0] * -d[0] + u[1] * -d[1] <= 0:
                    continue                                            # not carrying on past the chord's end
                nl = (-u[1], u[0])
                sgn = 1.0 if nl[0] * -nx + nl[1] * -ny > 0 else -1.0      # its edge on the outline's side
                P0 = (joints[q][0] + sgn * nl[0] * w2 / 2, joints[q][1] + sgn * nl[1] * w2 / 2)
                X = _isect((fx, fy), (ny, -nx), P0, u)
                if X is not None and (best is None or u[0] * -d[0] + u[1] * -d[1] > best[0]):
                    best = (u[0] * -d[0] + u[1] * -d[1], X)
            if best is not None:
                X = best[1]
                s = (X[0] - E[0]) * d[0] + (X[1] - E[1]) * d[1]           # < 0: past the end
                B[i][end] = (E[0] + s * d[0], E[1] + s * d[1])
    return [tuple(b) for b in B]


def outline_2d(joints, widths, members):
    """Polygons (lists of points) whose union is the truss before clipping and rounding: one rectangle per
    member, one mitred node polygon per joint of two or more members."""
    B = chord_mitres(joints, widths, members, bars(joints, widths, members))
    polys = []
    ends = {}
    for i, (m, (A, C, w)) in enumerate(zip(members, B)):
        d = _unit((C[0] - A[0], C[1] - A[1]))
        n = (-d[1] * w / 2, d[0] * w / 2)
        polys.append([(A[0] + n[0], A[1] + n[1]), (C[0] + n[0], C[1] + n[1]),
                      (C[0] - n[0], C[1] - n[1]), (A[0] - n[0], A[1] - n[1])])
        ends.setdefault(m["a"], []).append((A, d, w))
        ends.setdefault(m["c"], []).append((C, (-d[0], -d[1]), w))
    for q, inc in ends.items():
        if len(inc) < 2:
            continue
        J = joints[q]
        inc = sorted(inc, key=lambda e: math.atan2(e[1][1], e[1][0]))
        wmax = max(e[2] for e in inc)
        node = []
        for k in range(len(inc)):
            (P0, d0, w0), (P1, d1, w1) = inc[k], inc[(k + 1) % len(inc)]
            l0 = (P0[0] - d0[1] * w0 / 2, P0[1] + d0[0] * w0 / 2)      # left end corner of member k (seen from
            r1 = (P1[0] + d1[1] * w1 / 2, P1[1] - d1[0] * w1 / 2)      # the joint), right one of member k + 1
            x = _isect(l0, d0, r1, d1)                                  # mitre point of the wedge between them
            if x is None or math.hypot(x[0] - J[0], x[1] - J[1]) > MITRE * wmax:
                node += [l0, r1]                                        # bevel (or parallel edges)
            else:                                                       # an outer mitre lies behind the end
                behind0 = (x[0] - l0[0]) * d0[0] + (x[1] - l0[1]) * d0[1] < 0       # corners: they join it
                behind1 = (x[0] - r1[0]) * d1[0] + (x[1] - r1[1]) * d1[1] < 0
                node += ([l0] if behind0 else []) + [x] + ([r1] if behind1 else [])
        if len(node) >= 3:
            polys.append(node)
    return polys


def truss_face(joints, widths, members, box, R):
    """The truss section as a CadQuery face: union of outline_2d, clipped to the box (x0, y0, x1, y1), concave
    corners inside the box rounded with R."""
    import cadquery as cq
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.gp import gp_Pnt
    from OCP.TopAbs import TopAbs_IN
    area = lambda p: 0.5 * abs(sum(p[i][0] * p[i - 1][1] - p[i - 1][0] * p[i][1] for i in range(len(p))))
    # the union as unit-thick prisms (OCC's 3D fuse merges coplanar faces reliably, its 2D fuse does not)
    prisms = [cq.Solid.extrudeLinear(cq.Face.makeFromWires(cq.Wire.makePolygon(
        [cq.Vector(x, y, 0) for x, y in p], close=True)), cq.Vector(0, 0, 1))
        for p in outline_2d(joints, widths, members) if area(p) > 1e-9]
    x0, y0, x1, y1 = box
    u = prisms[0].fuse(*prisms[1:]).clean() if len(prisms) > 1 else prisms[0]
    u = u.intersect(cq.Solid.makeBox(x1 - x0, y1 - y0, 1.0, cq.Vector(x0, y0, 0))).clean()
    fs = [g for g in u.Faces() if abs(g.Center().z) < 1e-9 and g.normalAt().z < -0.5]
    if len(fs) != 1:
        raise ValueError("the truss section is %d separate regions" % len(fs))
    f = cq.Face.makeFromWires(fs[0].outerWire(), fs[0].innerWires())
    if R <= 0:
        return f
    on_box = lambda p: min(abs(p.x - x0), abs(p.x - x1), abs(p.y - y0), abs(p.y - y1)) < 1e-6
    concave = []
    for v in f.Vertices():
        p = v.Center()
        if on_box(p):
            continue
        es = [e for e in f.Edges() if any((w.Center() - p).Length < 1e-7 for w in e.Vertices())]
        if len(es) != 2:
            continue
        qs = []
        for e in es:
            t0, t1 = (0.0, 1e-3) if (e.startPoint() - p).Length < 1e-7 else (1.0, 1.0 - 1e-3)
            qs.append(e.positionAt(t1 if t0 == 0.0 else t1))
        mid = (qs[0] + qs[1]) * 0.5
        if BRepClass_FaceClassifier(f.wrapped, gp_Pnt(mid.x, mid.y, mid.z), 1e-9).State() != TopAbs_IN:
            concave.append(v)                  # the chord between the two edges leaves the face: concave
    if not concave:
        return f
    # each wire rounded on its own (Face.fillet2D drops the holes of a face)
    key = lambda v: tuple(round(c, 6) for c in v.Center().toTuple())
    ck = {key(v) for v in concave}
    ws = []
    for w in [f.outerWire()] + f.innerWires():
        vs = [v for v in w.Vertices() if key(v) in ck]
        ws.append(w.fillet2D(R, vs) if vs else w)
    return cq.Face.makeFromWires(ws[0], ws[1:])


# ---- part -----------------------------------------------------------------
face = truss_face(JOINTS, WIDTHS, MEMBERS, BOX, R)
part = cq.Solid.extrudeLinear(face, cq.Vector(0, 0, t))
part = part.fuse(part.translate(cq.Vector(0, 0, H - t)))
for x0, y0, x1, y1 in WEBS:
    part = part.fuse(cq.Solid.makeBox(x1 - x0, y1 - y0, H - 2 * t, cq.Vector(x0, y0, t)))
part = part.clean()
print("volume %.3f mm^3, valid %s" % (part.Volume(), part.isValid()))
cq.exporters.export(cq.Workplane().add(part), sys.argv[1] if len(sys.argv) > 1 else "cantilever7_d4.step")
