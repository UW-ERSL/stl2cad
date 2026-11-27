"""Member-driven CadQuery and SolidWorks models of a planar-truss TO (truss2d.fit): cantilever Case 7.

The truss is described by its parameters: joints (x, y), width classes w_k, the corner radius R, the bar depth t,
and the domain (box, depth H). Each member is a bar of its width between its joints; a chord (a member along the
outline) is a strip flush with the outline, its centre line w/2 inside it. Joints are mitred: at a joint the edges
of consecutive members (by angle) run on to their intersection (a bevel through the two end corners when that
point is farther than MITRE x the widest bar); a chord end at a joint where another member carries on past it
ends where its outer edge meets that member's outer edge. The union is clipped to the domain and every concave corner inside
it is rounded with R. The section is extruded t on both faces (z 0..t and H-t..H); webs are rectangles flush with
the outline, between the faces.

The same geometry code (GEOM) runs here and in the CadQuery script, so the script rebuilds the part from the
parameters: edit a joint, a width, R or t and run it again.

SolidWorks (one sketch, both faces):
  Truss    sketch on the Front Plane: the section's outline (lines, arcs) and one construction centerline per
           member (not chords) between its joints. Joints: the centerlines' end points, coincident at each joint,
           each joint dimensioned from the origin (horizontal, vertical). Each member edge parallel to its
           centerline at a distance linked to "w_k" / 2; each chord's inner edge parallel to the outline at a
           distance linked to "w_k"; outline edges horizontal or vertical, collinear pieces, on the origin or
           dimensioned from it; arcs tangent to their neighbours, equal, one radius linked to "R"; free chord ends
           normal to the outline, dimensioned from the origin.
  Truss face 1   Extruded Boss of Truss, depth "t" from the Front Plane;
  Truss face 2   Extruded Boss of the same sketch (shared), from an offset H, reversed, depth "t";
                 if SolidWorks refuses the shared sketch: a copy of the sketch (geometry only) on the Front Plane;
  Web n    rectangle on a plane at z = "t", Extruded Boss H - 2 "t" towards the other face.
Global variables (Equation Manager): t, R, w_1 ... w_k.

python swtruss.py ORIGINAL.stl TO.stl NAME [--units m] [--detail 4]
Writes NAME.py (CadQuery), NAME.step, NAME.stl, NAME_solidworks.py, NAME_standin.log, NAME_truss.txt.
"""
import math
import os
import re
import subprocess
import sys

import numpy as np

GEOM = r'''
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
'''

exec(GEOM)


# ---------------------------------------------------------------- parameters from truss2d.fit
def params(R, box, H, digits=4):
    """Script parameters from truss2d.fit's result: joints named J1.. (by x, then y), width classes w1.. (by
    value), members with chord edges; webs as rectangles."""
    used = sorted(R['J'], key=lambda q: (round(float(R['J'][q][0]), 3), round(float(R['J'][q][1]), 3)))
    name = {q: 'J%d' % (i + 1) for i, q in enumerate(used)}
    joints = {name[q]: tuple(round(float(c), digits) for c in R['J'][q]) for q in used}
    ws = sorted({round(m[2], digits) for m in R['members']})
    wname = {w: 'w%d' % (i + 1) for i, w in enumerate(ws)}
    members = []
    for a, c, w, edge, ext in R['members']:
        m = dict(a=name[a], c=name[c], w=wname[round(w, digits)])
        if edge is not None:
            f, n = edge
            m['edge'] = (tuple(round(float(v), 6) for v in f), tuple(round(float(v), 9) for v in n))
            m['ext'] = tuple(round(float(e), digits) for e in ext)
        members.append(m)
    webs = []
    for g, _ in R['webs']:
        x0, y0, x1, y1 = g.bounds
        webs.append(tuple(round(float(v), digits) for v in (x0, y0, x1, y1)))
    used_j = {m['a'] for m in members} | {m['c'] for m in members}
    return dict(joints={k: v for k, v in joints.items() if k in used_j}, widths={wname[w]: w for w in ws},
                members=members, box=tuple(float(v) for v in box), H=float(H), t=round(float(R['t']), digits),
                R=float(R['fillet']), webs=webs)


def solid(P):
    """Truss faces and webs (CadQuery), and the solids after each feature (for the SolidWorks checks)."""
    import cadquery as cq
    f = truss_face(P['joints'], P['widths'], P['members'], P['box'], P['R'])
    t, H = P['t'], P['H']
    f1 = cq.Solid.extrudeLinear(f, cq.Vector(0, 0, t))
    f2 = f1.translate(cq.Vector(0, 0, H - t))
    steps = [('Truss face 1', f1), ('Truss face 2', f2)]
    for i, (x0, y0, x1, y1) in enumerate(P['webs']):
        steps.append(('Web %d' % (i + 1), cq.Solid.makeBox(x1 - x0, y1 - y0, H - 2 * t, cq.Vector(x0, y0, t))))
    part, out = None, []
    for nm, s in steps:
        part = s if part is None else part.fuse(s).clean()
        out.append((nm, part))
    return f, part, out


# ---------------------------------------------------------------- CadQuery script
def write_cq(path, P, title):
    L = ['"""%s: planar-truss model of a topology-optimized part (stl2cad, to/swtruss.py); units mm.' % title, '',
         'Parameters: joints, member width classes, corner radius R, bar depth t; edit and run again:',
         '    python %s [out.step]' % os.path.basename(path),
         'Members are bars between their joints; chords are strips flush with the outline (edge: a point on the',
         'outline and the inward normal; ext: run-on at free ends). Joints are mitred, concave corners rounded with R.',
         'The section is extruded t on both faces of the domain (depth H); webs are boxes between the faces.',
         '"""', 'import math', 'import sys', '', 'import cadquery as cq', '',
         '# ---- parameters ' + '-' * 60,
         'JOINTS = {']
    for k, v in P['joints'].items():
        L.append('    %r: (%.4f, %.4f),' % (k, *v))
    L += ['}', 'WIDTHS = {']
    for k, v in P['widths'].items():
        L.append('    %r: %.4f,' % (k, v))
    L += ['}', 'MEMBERS = [']
    for m in P['members']:
        L.append('    %r,' % m)
    L += [']', 'R = %.4g          # corner radius' % P['R'], 't = %.4f        # bar depth (each face)' % P['t'],
          'BOX = %r    # domain section (x0, y0, x1, y1)' % (P['box'],), 'H = %.4f        # domain depth' % P['H'],
          'WEBS = %r    # (x0, y0, x1, y1), between the faces' % (P['webs'],), '',
          '# ---- geometry ' + '-' * 62, GEOM.strip(), '', '',
          '# ---- part ' + '-' * 65,
          'face = truss_face(JOINTS, WIDTHS, MEMBERS, BOX, R)',
          'part = cq.Solid.extrudeLinear(face, cq.Vector(0, 0, t))',
          'part = part.fuse(part.translate(cq.Vector(0, 0, H - t)))',
          'for x0, y0, x1, y1 in WEBS:',
          '    part = part.fuse(cq.Solid.makeBox(x1 - x0, y1 - y0, H - 2 * t, cq.Vector(x0, y0, t)))',
          'part = part.clean()',
          'print("volume %.3f mm^3, valid %s" % (part.Volume(), part.isValid()))',
          'cq.exporters.export(cq.Workplane().add(part), sys.argv[1] if len(sys.argv) > 1 else "%s.step")' % title, '']
    open(path, 'w').write('\n'.join(L))


# ---------------------------------------------------------------- SolidWorks sketch plan
def tag_edges(face, P):
    """The section's edges with their role: ('D', side) outline, ('M', i, s) member i's edge (s = +1 left, -1
    right of a -> c), ('C', i) chord i's inner edge, ('E', i) a chord's free end, ('A',) a fillet arc,
    ('B',) another line (a bevel)."""
    x0, y0, x1, y1 = P['box']
    B = bars(P['joints'], P['widths'], P['members'])
    out = []
    for w in [face.outerWire()] + face.innerWires():
        loop = []
        for e in w.Edges():
            p, q = e.startPoint(), e.endPoint()
            if e.geomType() == 'CIRCLE':
                loop.append(('A', (p.x, p.y), (e.positionAt(0.5).x, e.positionAt(0.5).y), (q.x, q.y), ('A',)))
                continue
            tag = ('B',)
            for side, val, k in ((0, y0, 1), (1, x1, 0), (2, y1, 1), (3, x0, 0)):
                if abs(p.toTuple()[k] - val) < 1e-6 and abs(q.toTuple()[k] - val) < 1e-6:
                    tag = ('D', side)
            if tag == ('B',):
                for i, (m, (A, C, wd)) in enumerate(zip(P['members'], B)):
                    d = _unit((C[0] - A[0], C[1] - A[1]))
                    dist = lambda r: (r.x - A[0]) * -d[1] + (r.y - A[1]) * d[0]
                    a, b = dist(p), dist(q)
                    if abs(a - b) < 1e-6 and abs(abs(a) - wd / 2) < 1e-6:
                        tag = ('C', i) if m.get('edge') is not None else ('M', i, 1 if a > 0 else -1)
                        break
                    along = lambda r: (r.x - A[0]) * d[0] + (r.y - A[1]) * d[1]
                    if m.get('edge') is not None and abs(a) <= wd / 2 + 1e-6 and abs(b) <= wd / 2 + 1e-6 \
                            and abs(along(p) - along(q)) < 1e-6:
                        tag = ('E', i)
                        break
            loop.append(('L', (p.x, p.y), (q.x, q.y), tag))
        # end points shared exactly by the edges that meet there (OCC's fillets leave differences of ~1e-7 mm)
        reps = []

        def snap(pt):
            for r in reps:
                if math.hypot(r[0] - pt[0], r[1] - pt[1]) < 1e-5:
                    return r
            reps.append(pt)
            return pt
        loop = [(e[0], snap(e[1]), e[2], e[3]) if e[0] == 'L' else (e[0], snap(e[1]), e[2], snap(e[3]), e[4])
                for e in loop]
        loop = [(e[0], e[1], snap(e[2]), e[3]) if e[0] == 'L' else e for e in loop]
        out.append(loop)
    return out


SW = r'''"""@NAME@: SolidWorks feature tree of a planar-truss TO (stl2cad, to/swtruss.py), member-driven sketch.

Feature tree:
  Truss          sketch on the Front Plane: the truss section (lines, tangent arcs) and one construction centerline
                 per member between its joints. Joints dimensioned from the origin; each member edge parallel to
                 its centerline at "w_k" / 2; each chord's inner edge at "w_k" from the outline; one corner radius
                 "R" (all fillet arcs equal); outline edges on the origin or dimensioned from it.
  Truss face 1   Extruded Boss, depth "t", from the Front Plane.
  Truss face 2   Extruded Boss of the same sketch from an offset H, reversed, depth "t" (else a copy of the sketch).
  Web n          rectangle on "Web plane" (z = "t"), Extruded Boss H - 2 "t" towards the other face.
Global variables (Equation Manager): t, R and the member widths w1 ... (edit them there).
After every feature the volume and centre of mass are compared with CadQuery's values (STEPS).

Needs: Windows, SolidWorks, Python with pywin32.   Run: python @NAME@.py [C:\\path\\part.SLDPRT]
"""
import math
import sys
import traceback

import pythoncom
import win32com.client
from win32com.client import VARIANT

sys.stdout.reconfigure(line_buffering=True)
NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
TOL = 2e-3                    # a feature's volume may differ from CadQuery's by this fraction
@DATA@


def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


def sw_constant(app, enum, member):
    import os
    exe = app.GetExecutablePath
    exe = str(exe() if callable(exe) else exe)
    for f in (exe, os.path.dirname(exe)):
        for sub in ((), ("api",)):
            path = os.path.join(f, *sub, "swconst.tlb")
            if not os.path.exists(path):
                continue
            tl = pythoncom.LoadTypeLib(path)
            for i in range(tl.GetTypeInfoCount()):
                if tl.GetDocumentation(i)[0] != enum:
                    continue
                ti = tl.GetTypeInfo(i)
                for j in range(ti.GetTypeAttr().cVars):
                    vd = ti.GetVarDesc(j)
                    if ti.GetNames(vd.memid)[0] == member:
                        return int(vd.value)
    raise LookupError(f"{enum}.{member} not found")


app = win32com.client.Dispatch("SldWorks.Application")
app.Visible = True
try:
    app.SetUserPreferenceToggle(10, False)        # no dimension input box
except Exception:
    pass
CONST = {}
for enum, member, default in (("swEndConditions_e", "swEndCondMidPlane", 6), ("swEndConditions_e", "swEndCondBlind", 0)):
    try:
        CONST[member] = sw_constant(app, enum, member)
    except Exception:
        CONST[member] = default
LOG = dict(rels=0, rels_failed=[], dims=0, dims_failed=[], eqs=0, eqs_failed=[])


def new_part():
    m = call0(app, "NewPart")
    try:
        m.Extension.SetUserPreferenceInteger(sw_constant(app, "swUserPreferenceIntegerValue_e", "swUnitSystem"),
                                             0, sw_constant(app, "swUnitSystem_e", "swUnitSystem_MMGS"))
    except Exception:
        pass
    return m


def mass(m):
    m.ForceRebuild3(False)
    mp = call0(m.Extension, "CreateMassProperty")
    if mp is None:
        return 0.0, (0.0, 0.0, 0.0)
    return mp.Volume / MM ** 3, tuple(c / MM for c in mp.CenterOfMass)


def check(m, step, f):
    ev, ec = STEPS[step]
    if f is None:
        print(f"{step}: NOT MADE")
        return False
    v, c = mass(m)
    ok = abs(v - ev) <= TOL * ev
    print(f"{step}: volume {v:.3f} mm^3 (expected {ev:.3f}); centre {math.dist(c, ec):.4f} mm off -> "
          f"{'OK' if ok else 'DIFFERS'}")
    return ok


def to_sk(T, p):
    x = [c * MM for c in p]
    return tuple((sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3))


def open_sketch(m, plane):
    sm = m.SketchManager
    for attempt in range(3):              # the first sketch call after NewPart can open nothing
        m.ClearSelection2(True)
        if isinstance(plane, str):
            m.Extension.SelectByID2(plane, "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
        else:
            plane.Select2(False, 0)
        sm.InsertSketch(True)
        sk = call0(sm, "ActiveSketch")
        if sk is not None:
            break
    else:
        raise RuntimeError("no sketch opened")
    T = list(sk.ModelToSketchTransform.ArrayData)
    sm.AddToDB = True
    helper = sm.CreateCircleByRadius(-0.5, -0.5, 0.0, 0.5 * MM)   # the first entity of a new sketch can be
    m.ClearSelection2(True)                                        # lost: a helper, deleted again
    if helper is not None and helper.Select4(False, NOTHING):
        m.EditDelete()
    return sm, T


def close_sketch(m, sm, name):
    sm.AddToDB = False
    m.ClearSelection2(True)
    sm.InsertSketch(True)
    f = m.FeatureByPositionReverse(0)
    f.Name = name
    return f


def arc_centre(s, mid, t):
    (ax, ay), (bx, by), (cx, cy) = s, mid, t
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    return ux, uy


def draw_loops(sm, T, loops):
    """Lines and arcs of the section (model mm); returns [[(entry, object)]] per loop."""
    out = []
    for loop in loops:
        objs = []
        for e in loop:
            if e[0] == "L":
                g = sm.CreateLine(*to_sk(T, (*e[1], 0.0)), *to_sk(T, (*e[2], 0.0)))
            else:
                c = arc_centre(e[1], e[2], e[3])
                cross = (e[2][0] - e[1][0]) * (e[3][1] - e[1][1]) - (e[2][1] - e[1][1]) * (e[3][0] - e[1][0])
                g = sm.CreateArc(*to_sk(T, (*c, 0.0)), *to_sk(T, (*e[1], 0.0)), *to_sk(T, (*e[3], 0.0)),
                                 1 if cross > 0 else -1)          # +1: counter-clockwise
            objs.append((e, g))
        out.append(objs)
    return out


def sel(m, objs):
    m.ClearSelection2(True)
    ok = True
    for k, o in enumerate(objs):
        if o == "ORIGIN":
            ok = ok and m.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0, k > 0, 0, NOTHING, 0)
        else:
            ok = ok and o is not None and o.Select4(k > 0, NOTHING)
    return ok


def rel(m, objs, r, what):
    try:
        if sel(m, objs):
            m.SketchAddConstraints(r); LOG["rels"] += 1
            return True
    except Exception as ex:
        what += f" ({type(ex).__name__})"
    LOG["rels_failed"].append(what)
    return False


def dim(m, objs, kind, text, what):
    try:
        d = None
        if sel(m, objs):
            tx, ty = text
            d = {"any": m.AddDimension2, "H": m.AddHorizontalDimension2, "V": m.AddVerticalDimension2}[kind](
                tx * MM, ty * MM, 0.0)
        if d is not None:
            LOG["dims"] += 1
            return call0(d.GetDimension2(0), "FullName")
    except Exception as ex:
        what += f" ({type(ex).__name__})"
    LOG["dims_failed"].append(what)
    return None


def point_at(objs, p, tol=1e-6):
    for e, g in objs:
        for nm in ("GetStartPoint2", "GetEndPoint2"):
            q = call0(g, nm)
            if q is not None and math.hypot(q.X / MM - p[0], q.Y / MM - p[1]) < tol:
                return q
    return None


def equation(m, text):
    try:
        em = call0(m, "GetEquationMgr")
        em.Add2(-1, text, True); LOG["eqs"] += 1
    except Exception as ex:
        LOG["eqs_failed"].append(f"{text} ({type(ex).__name__})")


def truss_sketch(m, dimension=True):
    sm, T = open_sketch(m, "Front Plane")
    E = draw_loops(sm, T, LOOPS)
    flat = [x for loop in E for x in loop]
    cl = {}
    for i, mb in enumerate(MEMBERS):          # construction centerlines (not chords)
        if mb.get("edge") is None:
            a, c = JOINTS[mb["a"]], JOINTS[mb["c"]]
            cl[i] = sm.CreateCenterLine(*to_sk(T, (*a, 0.0)), *to_sk(T, (*c, 0.0)))
    sm.AddToDB = False
    links = []
    if dimension:
        # joints: centerline ends coincident; one end per joint dimensioned from the origin
        for j, p in JOINTS.items():
            pts = [point_at([(None, cl[i])], p) for i in cl if j in (MEMBERS[i]["a"], MEMBERS[i]["c"])]
            pts = [q for q in pts if q is not None]
            for q in pts[1:]:
                rel(m, [pts[0], q], "sgCOINCIDENT", f"joint {j}")
            if pts:
                if abs(p[0]) > 1e-9:
                    dim(m, [pts[0], "ORIGIN"], "H", (p[0] / 2, p[1] + 4), f"joint {j} x")
                else:
                    rel(m, [pts[0], "ORIGIN"], "sgVERTICALPOINTS2D", f"joint {j} on x = 0")
                if abs(p[1]) > 1e-9:
                    dim(m, [pts[0], "ORIGIN"], "V", (p[0] + 4, p[1] / 2), f"joint {j} y")
                else:
                    rel(m, [pts[0], "ORIGIN"], "sgHORIZONTALPOINTS2D", f"joint {j} on y = 0")
        # outline edges: horizontal / vertical; the first piece of each side on the origin or dimensioned,
        # the others collinear with it
        first = {}
        for e, g in flat:
            if e[0] == "L" and e[3][0] == "D":
                side = e[3][1]
                rel(m, [g], "sgHORIZONTAL2D" if side in (0, 2) else "sgVERTICAL2D", f"outline side {side}")
                if side in first:
                    rel(m, [first[side], g], "sgCOLINEAR", f"outline side {side} piece")
                    continue
                first[side] = g
                if side == 0 and abs(BOX[1]) < 1e-9 or side == 3 and abs(BOX[0]) < 1e-9:
                    rel(m, [g, "ORIGIN"], "sgCOINCIDENT", f"outline side {side} on the origin")
                else:
                    q = call0(g, "GetStartPoint2")
                    dim(m, [q, "ORIGIN"], "V" if side in (0, 2) else "H", (e[1][0] + 4, e[1][1] + 4),
                        f"outline side {side}")
        # member edges: parallel to the centerline at w/2; chord inner edges parallel to the outline at w
        for e, g in flat:
            if e[0] != "L":
                continue
            tg = e[3]
            if tg[0] == "M":
                i = tg[1]
                rel(m, [g, cl[i]], "sgPARALLEL", f"member {i} edge")
                mid = ((e[1][0] + e[2][0]) / 2, (e[1][1] + e[2][1]) / 2)
                nm = dim(m, [g, cl[i]], "any", (mid[0] + 2, mid[1] + 2), f"member {i} edge offset")
                if nm:
                    links.append((nm, '"%s" / 2' % MEMBERS[i]["w"]))
            elif tg[0] == "C":
                i = tg[1]
                side = CHORD_SIDE[i]
                if side in first:
                    rel(m, [g, first[side]], "sgPARALLEL", f"chord {i} edge")
                    mid = ((e[1][0] + e[2][0]) / 2, (e[1][1] + e[2][1]) / 2)
                    nm = dim(m, [g, first[side]], "any", (mid[0] + 2, mid[1] + 2), f"chord {i} width")
                    if nm:
                        links.append((nm, '"%s"' % MEMBERS[i]["w"]))
            elif tg[0] == "E":
                i = tg[1]
                side = CHORD_SIDE[i]
                rel(m, [g], "sgVERTICAL2D" if side in (0, 2) else "sgHORIZONTAL2D", f"chord {i} end")
                q = call0(g, "GetStartPoint2")
                dim(m, [q, "ORIGIN"], "H" if side in (0, 2) else "V", (e[1][0] + 2, e[1][1] - 4), f"chord {i} end")
        # fillet arcs: tangent to their neighbours, equal, one radius
        arc0 = None
        for loop in E:
            for k, (e, g) in enumerate(loop):
                if e[0] != "A":
                    continue
                for nb in (loop[k - 1], loop[(k + 1) % len(loop)]):
                    if nb[0][0] == "L":
                        rel(m, [g, nb[1]], "sgTANGENT", "fillet tangency")
                if arc0 is None:
                    arc0 = g
                    nm = dim(m, [g], "any", (e[2][0] + 3, e[2][1] + 3), "corner radius")
                    if nm:
                        links.append((nm, '"R"'))
                else:
                    rel(m, [arc0, g], "sgSAMELENGTH", "equal radius")
    m.ClearSelection2(True)
    sm.InsertSketch(True)
    f = m.FeatureByPositionReverse(0)
    f.Name = "Truss" if dimension else "Truss copy"
    return f, links


def boss(m, sk, depth, offset=None, reverse=False):
    m.ClearSelection2(True)
    sk.Select2(False, 0)
    t0, off = (3, offset * MM) if offset is not None else (0, 0.0)   # 3: start offset
    return m.FeatureManager.FeatureExtrusion3(True, False, reverse, CONST["swEndCondBlind"], 0, depth * MM, 0.0,
                                              False, False, False, False, 0.0, 0.0, False, False, False, False,
                                              True, True, True, t0, off, False)


def plane_at(m, z, name):
    """Plane parallel to the Front Plane at z (mm), its side checked through a sketch transform (as the bike
    hook's plane_at_z, which SolidWorks 2025 built)."""
    for flip in (0, 256):
        m.ClearSelection2(True)
        m.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, NOTHING, 0)
        if m.FeatureManager.InsertRefPlane(8 | flip, z * MM, 0, 0, 0, 0) is None:
            continue
        f = m.FeatureByPositionReverse(0)
        sm, T = open_sketch(m, f)
        off = to_sk(T, (0.0, 0.0, z))[2] / MM
        sm.AddToDB = False
        sm.InsertSketch(True)
        sk = m.FeatureByPositionReverse(0)
        m.ClearSelection2(True); sk.Select2(False, 0); m.Extension.DeleteSelection2(0)
        if abs(off) < 1e-6:
            f.Name = name
            return f
        m.ClearSelection2(True); f.Select2(False, 0); m.Extension.DeleteSelection2(0)
    raise RuntimeError(f"no plane at z = {z}")


def web(m, plane, k, box):
    x0, y0, x1, y1 = box
    sm, T = open_sketch(m, plane)
    z = T_DEPTH
    c = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    ln = [sm.CreateLine(*to_sk(T, (*c[i], z)), *to_sk(T, (*c[(i + 1) % 4], z))) for i in range(4)]
    sm.AddToDB = False
    pts = [call0(g, "GetStartPoint2") for g in ln]
    for i, g in enumerate(ln):
        rel(m, [g], "sgHORIZONTAL2D" if i % 2 == 0 else "sgVERTICAL2D", f"web {k} side {i}")
    dim(m, [ln[0]], "any", ((x0 + x1) / 2, y0 - 3), f"web {k} length")
    dim(m, [ln[3]], "any", (x0 - 3, (y0 + y1) / 2), f"web {k} height")
    # position: a web at the domain's end is placed by that end (on the origin, or dimensioned from it)
    q, qx = (pts[1], x1) if abs(x1 - BOX[2]) < 1e-9 and abs(x0 - BOX[0]) > 1e-9 else (pts[0], x0)
    if abs(qx) > 1e-9:
        dim(m, [q, "ORIGIN"], "H", (qx / 2, y0 - 6), f"web {k} position")
    else:
        rel(m, [q, "ORIGIN"], "sgVERTICALPOINTS2D", f"web {k} on x = 0")
    rel(m, [pts[0], "ORIGIN"], "sgHORIZONTALPOINTS2D", f"web {k} on y = 0")
    m.ClearSelection2(True)
    sm.InsertSketch(True)
    sk = m.FeatureByPositionReverse(0)
    sk.Name = f"Web {k} sketch"
    f = None
    for rev in (False, True):                 # towards +z from the plane at z = t; the side checked by volume
        f = boss(m, sk, H - 2 * T_DEPTH, reverse=rev)
        if f is None:
            continue
        f.Name = f"Web {k}"
        v, c = mass(m)
        if abs(v - STEPS[f"Web {k}"][0]) <= TOL * v:
            return f
        m.ClearSelection2(True); f.Select2(False, 0); m.Extension.DeleteSelection2(0)
    return None


def status(f):
    try:
        st = call0(call0(f, "GetSpecificFeature2"), "GetConstrainedStatus")
        return {1: "unknown", 2: "under defined", 3: "fully defined", 4: "over defined", 5: "no solution",
                6: "invalid solution"}.get(st, str(st))
    except Exception as ex:
        return f"not available ({type(ex).__name__})"


def build():
    m = new_part()
    ok = True
    sk, links = truss_sketch(m)
    print(f"Truss sketch: {status(sk)}; {LOG['rels']} relations, {LOG['dims']} dimensions")
    f1 = boss(m, sk, T_DEPTH)
    if f1 is not None:
        f1.Name = "Truss face 1"
    ok &= check(m, "Truss face 1", f1)
    f2 = boss(m, sk, T_DEPTH, offset=H, reverse=True)            # the same sketch, from z = H towards the part
    if f2 is None or not check(m, "Truss face 2", f2):
        if f2 is not None:
            m.ClearSelection2(True); f2.Select2(False, 0); m.Extension.DeleteSelection2(0)
        print("Truss face 2: the shared sketch was not used; a copy of the sketch instead")
        sk2, _ = truss_sketch(m, dimension=False)
        f2 = boss(m, sk2, T_DEPTH, offset=H, reverse=True)
        ok &= check(m, "Truss face 2", f2)
    if f2 is not None:
        f2.Name = "Truss face 2"
    wp = plane_at(m, T_DEPTH, "Web plane") if WEBS else None
    for k, box in enumerate(WEBS, 1):
        ok &= check(m, f"Web {k}", web(m, wp, k, box))
    # global variables and links
    equation(m, f'"t" = {T_DEPTH:g}mm')
    equation(m, f'"R" = {R_CORNER:g}mm')
    for w, v in WIDTHS.items():
        equation(m, f'"{w}" = {v:g}mm')
    for nm, expr in links:
        equation(m, f'"{nm}" = {expr}')
    for nm in ("D1@Truss face 1", "D1@Truss face 2"):
        equation(m, f'"{nm}" = "t"')
    if WEBS:
        equation(m, '"D1@Web plane" = "t"')
    for k in range(1, len(WEBS) + 1):
        equation(m, f'"D1@Web {k}" = {H:g}mm - 2 * "t"')
    try:
        m.EditRebuild3()
    except Exception:
        pass
    v, c = mass(m)
    ev, ec = STEPS[list(STEPS)[-1]]
    print(f"after the equations: volume {v:.3f} mm^3 (expected {ev:.3f})")
    print(f"relations: {LOG['rels']} added" + (f"; failed: {LOG['rels_failed']}" if LOG['rels_failed'] else ""))
    print(f"dimensions: {LOG['dims']} added" + (f"; failed: {LOG['dims_failed']}" if LOG['dims_failed'] else ""))
    print(f"equations: {LOG['eqs']} added" + (f"; failed: {LOG['eqs_failed']}" if LOG['eqs_failed'] else ""))
    print(f"independent parameters: {N_PARAMS}")
    print("Truss sketch: " + status(sk))
    print("all features " + ("agree" if ok else "do NOT all agree"))
    return m


if __name__ == "__main__":
    try:
        m = build()
        if len(sys.argv) > 1:
            m.SaveAs3(sys.argv[1], 0, 2)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
'''


def write_sw(path, P, loops, steps, title):
    chord_side = {}
    x0, y0, x1, y1 = P['box']
    for i, m in enumerate(P['members']):
        if m.get('edge') is not None:
            (fx, fy), (nx, ny) = m['edge']
            chord_side[i] = 0 if ny > 0.5 else 2 if ny < -0.5 else 3 if nx > 0.5 else 1
    n_joint = sum((abs(p[0]) > 1e-9) + (abs(p[1]) > 1e-9) for j, p in P['joints'].items()
                  if any(j in (m['a'], m['c']) for m in P['members'] if m.get('edge') is None))
    n_end = sum(1 for lp in loops for e in lp if e[0] == 'L' and e[3][0] == 'E')
    n_web = sum(2 + (abs(b[0] - x0) > 1e-9) for b in P['webs'])
    n_params = n_joint + len(P['widths']) + 2 + n_end + n_web
    D = ['H = %.6f                 # domain depth (mm)' % P['H'],
         'BOX = %r' % (P['box'],),
         'T_DEPTH = %.6f            # global variable "t"' % P['t'],
         'R_CORNER = %.6g            # global variable "R"' % P['R'],
         'WIDTHS = %r' % (P['widths'],),
         'JOINTS = %r' % ({k: tuple(v) for k, v in P['joints'].items()},),
         'MEMBERS = %r' % (P['members'],),
         'CHORD_SIDE = %r            # outline side of each chord (0 y0, 1 x1, 2 y1, 3 x0)' % chord_side,
         'WEBS = %r' % (P['webs'],),
         'N_PARAMS = %d             # joints %d, widths %d, R and t, chord ends %d, webs %d' %
         (n_params, n_joint, len(P['widths']), n_end, n_web),
         'LOOPS = ['] + ['    %r,' % lp for lp in loops] + [']', 'STEPS = {'] + \
        ['    %r: (%.4f, %r),' % (nm, v, c) for nm, v, c in steps] + ['}']
    open(path, 'w').write(SW.replace('@NAME@', title).replace('@DATA@', '\n'.join(D)))
    return n_params


def verify(script, ref, timeout=900):
    """Run the SolidWorks script on the stand-in (swmock); compare its part with ref."""
    import tempfile
    import cadquery as cq
    import stl2cad
    d = tempfile.mkdtemp(prefix='swtruss_')
    py, out = os.path.join(d, 'part.py'), os.path.join(d, 'part.step')
    open(py, 'w').write(open(script).read())
    pkg = os.path.dirname(stl2cad.__file__)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.path.join(pkg, 'swmock'), os.path.dirname(pkg)]),
               SWMOCK_CUT_REVERSED='1', SWMOCK_BOSS_REVERSED='0', SWMOCK_EDGE_SIGHT='1')
    r = subprocess.run([sys.executable, py, out], capture_output=True, text=True, timeout=timeout, env=env)
    log = r.stdout + r.stderr
    if r.returncode != 0 or not os.path.exists(out):
        return False, log
    s = cq.importers.importStep(out).val()
    dv = s.cut(ref).Volume() + ref.cut(s).Volume()
    ok = abs(s.Volume() - ref.Volume()) < 1e-3 * ref.Volume() and dv < 1e-3 * ref.Volume()
    if not ok:              # OCC's booleans of the two solids can fail: compare volume, area, centre, inertia
        import swplates
        if swplates.invariants(s, ref, size=ref.BoundingBox().DiagonalLength):
            return True, log + '\nstand-in: volume %.3f, reference %.3f mm^3; booleans failed, volume, area, ' \
                               'centre of mass and inertia agree' % (s.Volume(), ref.Volume())
    return ok, log + '\nstand-in: volume %.3f, reference %.3f mm^3, difference volume %.3f mm^3' % (
        s.Volume(), ref.Volume(), dv)


def main():
    import argparse
    import trimesh
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import truss2d
    ap = argparse.ArgumentParser()
    ap.add_argument('original'); ap.add_argument('to'); ap.add_argument('name')
    ap.add_argument('--units', default='m'); ap.add_argument('--detail', type=int, default=4, choices=range(1, 6))
    a = ap.parse_args()
    sc = {'m': 1000.0, 'mm': 1.0, 'in': 25.4}[a.units]
    orig = trimesh.load(a.original); orig.apply_scale(sc)
    to = trimesh.load(a.to); to.apply_scale(sc)
    to = max(to.split(only_watertight=False), key=lambda b: len(b.faces))
    R = truss2d.fit(to, orig, a.detail)
    R['fillet'] = float('%.2g' % (truss2d.FILLET_SHARE * min(m[2] for m in R['members'])))
    k = R['axis']
    if k != 2:
        sys.exit('swtruss: the mirror axis must be z (the Front Plane sketch); got %s' % 'xyz'[k])
    ob = orig.bounds
    P = params(R, (ob[0][0], ob[0][1], ob[1][0], ob[1][1]), ob[1][2] - ob[0][2])
    f0 = truss_face(P['joints'], P['widths'], P['members'], P['box'], P['R'])
    P['t'] = round(R['v_truss'] / f0.Area(), 4)          # flat bars of the face's truss volume
    print('bar depth t = %.4f mm (truss volume of one face over the section area)' % P['t'])
    face, part, steps = solid(P)
    write_cq(a.name + '.py', P, a.name)
    r = subprocess.run([sys.executable, a.name + '.py', a.name + '.step'], capture_output=True, text=True)
    print('CadQuery script: ' + (r.stdout.strip() or r.stderr.strip()[-300:]))
    import cadquery as cq
    cq.exporters.export(cq.Workplane().add(part), a.name + '.stl', tolerance=0.02, angularTolerance=0.2)
    loops = tag_edges(face, P)
    tags = [e[-1][0] for lp in loops for e in lp]
    print('section: %d loops, %d edges (%s)' % (len(loops), len(tags),
          ', '.join('%s %d' % (t, tags.count(t)) for t in 'DMCEAB' if tags.count(t))))
    st = [(nm, s.Volume(), tuple(round(v, 4) for v in s.Center().toTuple())) for nm, s in steps]
    npar = write_sw(a.name + '_solidworks.py', P, loops, st, a.name + '_solidworks')
    ok, log = verify(a.name + '_solidworks.py', part)
    open(a.name + '_standin.log', 'w').write(log)
    print('SolidWorks script: %s_solidworks.py; stand-in: %s; %s' % (a.name, 'ok' if ok else 'FAILED',
                                                                      log.strip().splitlines()[-1]))
    # report against the TO
    import trimesh as tm
    m = tm.load(a.name + '.stl')
    print('volume: TO %.0f, model %.0f mm^3 (%+.1f %%); independent parameters %d'
          % (to.volume, m.volume, 100 * (m.volume / to.volume - 1), npar))
    L = ['%s: planar truss (to/swtruss.py, detail %d); units mm.' % (a.name, a.detail),
         'Bar depth t = %.4f on both faces; corner radius R = %.4g; domain depth H = %.4f.' % (P['t'], P['R'], P['H']),
         '', 'joints:'] + ['  %-4s %9.4f %9.4f' % (j, *p) for j, p in P['joints'].items()] + \
        ['', 'widths:'] + ['  %-4s %.4f' % kv for kv in P['widths'].items()] + \
        ['', 'members:'] + ['  %s - %s  %s%s' % (mb['a'], mb['c'], mb['w'], '  chord, run-on %s' % (mb['ext'],)
                                                 if mb.get('edge') is not None else '') for mb in P['members']] + \
        ['', 'webs (x0, y0, x1, y1):'] + ['  %r' % (w,) for w in P['webs']]
    open(a.name + '_truss.txt', 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
