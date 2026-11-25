"""SolidWorks script (pywin32) for a truss-like TO part built from swept members: hubs, one Loft per member
along its curved centre line, joint bosses, the original's holes.

Model (members.py): each member's centre line is lines and tangent arcs on the mid-plane; its section is an
ellipse normal to the path, w (in-plane) x d (along the hub axis), changing linearly from the start to the end.
In SolidWorks (calls checked with sw_probe_sweep.py, SolidWorks 2025):
  Hub n          Extruded Boss of a circle (Front Plane), the hub's height;
  Mid-plane      reference plane at the mid-plane level (InsertRefPlane, distance from the Front Plane);
  Member n       path sketch on the mid-plane; a plane normal to the path at each end (InsertRefPlane:
                 perpendicular to the end segment, coincident with the end point); an ellipse sketch on each;
                 Loft from start to end section with the path as centerline (InsertProtrusionBlend2;
                 sections Mark 1, centerline Mark 4). SolidWorks' loft and CadQuery's swept sections agree
                 to 0.04 % in volume (probe), so a loft is checked to 1 % of its own volume, and later
                 features against the expected values corrected by the lofts' small differences;
  Joint n        node at each joint not on a hub: ellipsoid of revolution (half-ellipse on a plane through the
                 joint parallel to the Right Plane, revolved about the vertical), 5 % over the member ends;
  Hole n         Cut-Extrude of the original's hole.
Sketch points on the new planes come from each sketch's ModelToSketchTransform (read in the probe: rotation
R applied as p R, then the translation), so no plane orientation is assumed.

python swsweep.py ORIGINAL_EXACT.py TO.stl NAME [--units in] [--detail 1-5] [--plot]
"""
import argparse, math, os, sys
import numpy as np
import cadquery as cq
import trimesh
from cadquery import Vector

from stl2cad import swgen


def mass_props(shape):
    """Volume and centre of mass (CadQuery's Gauss integration). stl2cad.core.mass_props' Gauss-Kronrod check
    (eps 1e-9) did not finish on the lofted members' B-spline faces."""
    c = shape.Center()
    return float(shape.Volume()), (c.x, c.y, c.z)
import hubweb, members, swplates, swbars

EPS = 1e-3
TAPER = 2.0      # largest ratio of a member's end widths (and end depths)
HUB_END = 0.6    # widest member end on a hub, x the hub's outer radius
HUB_FLARE = 1.25 # a member's end on a hub at most this x its other end (width and depth)
# detail level: classes of section widths and depths (0: every member end its own), significant digits
DETAIL = {1: dict(widths=1, depths=1, digits=2), 2: dict(widths=2, depths=1, digits=2),
          3: dict(widths=2, depths=2, digits=2), 4: dict(widths=3, depths=3, digits=3),
          5: dict(widths=0, depths=0, digits=3)}


def _ellipse_axes(w, d):
    """(major, minor, major along n?): semi-axes of the w x d ellipse; equal axes nudged apart (SolidWorks'
    ellipse needs a major axis)."""
    a, b = w / 2, d / 2
    if abs(a - b) < 0.02 * max(a, b):
        b = 0.98 * a
    return (a, b, True) if a >= b else (b, a, False)


def tangent(path, end):
    e = path[0] if end == 0 else path[-1]
    if e[0] == 'line':
        d = np.subtract(e[2], e[1])
    else:
        c = np.array(members.S_circle(*e[1:])); p = np.array(e[1] if end == 0 else e[3])
        r = p - c; d = np.array([-r[1], r[0]])
        nxt = np.array(e[2]) - np.array(e[1]) if end == 0 else np.array(e[3]) - np.array(e[2])
        if d @ nxt < 0:
            d = -d
    return d / np.linalg.norm(d)


def sections(m, z):
    """Start and end section of member m: (centre, tangent, major point, minor point) in model mm."""
    out = []
    for end, (w, d) in ((0, (m['w0'], m['d0'])), (1, (m['w1'], m['d1']))):
        p = m['path'][0][1] if end == 0 else m['path'][-1][-1]
        t = tangent(m['path'], end)
        n = np.array([-t[1], t[0], 0.0]); zz = np.array([0, 0, 1.0])
        c = np.array([p[0], p[1], z])
        a, b, along_n = _ellipse_axes(w, d)
        A, B = (n, zz) if along_n else (zz, n)
        out.append((c, np.r_[t, 0.0], c + a * A, c + b * B, a, b))
    return out


def fit(original, to, detail):
    hw = hubweb.detect(original, to)
    if hw is None:
        raise ValueError('no hubs: the original has no hole that the TO surrounds')
    R = members.members(hw['shadow'], hw['hubs'], to, hw['axis'], detail)
    z = R['mid']
    lo_h, hi_h = min(h['lo'] for h in hw['hubs']), max(h['hi'] for h in hw['hubs'])
    for m in R['members']:
        for k in ('d0', 'd1'):
            m[k] = min(m[k], 2 * min(z - lo_h, hi_h - z) - 0.2)
        # widest end on a hub: HUB_END x its outer radius (a wider end meeting the hub at a slant stuck out past it)
        # and no flare towards a hub: at most HUB_FLARE x the other end, in width and depth (the TO's fan at the
        # hub left wide flat end faces beside it in SolidWorks)
        for end, j in ((0, m['i']), (1, m['j'])):
            if j in R['hub_of'] and not (m['i'] in R['hub_of'] and m['j'] in R['hub_of']):
                a_, b_ = ('0', '1') if end == 0 else ('1', '0')
                m['w' + a_] = min(m['w' + a_], HUB_END * hw['hubs'][R['hub_of'][j]]['ro'], HUB_FLARE * m['w' + b_])
                m['d' + a_] = min(m['d' + a_], HUB_FLARE * m['d' + b_])
        # at most a factor TAPER between the two ends (the fan at a hub made 11 mm-wide flat flaps), mean kept
        for a_, b_ in (('w0', 'w1'), ('d0', 'd1')):
            lo_, hi_ = sorted((m[a_], m[b_]))
            if hi_ > TAPER * lo_:
                mean = 0.5 * (lo_ + hi_)
                small, big = 2 * mean / (1 + TAPER), 2 * mean * TAPER / (1 + TAPER)
                m[a_], m[b_] = (small, big) if m[a_] < m[b_] else (big, small)
        # the ellipse's major axis must stay on the same side along the member (else the loft twists)
        if (m['w0'] - m['d0']) * (m['w1'] - m['d1']) < 0:
            m['d0'] = m['d1'] = 0.5 * (m['d0'] + m['d1'])
            if (m['w0'] - m['d0']) * (m['w1'] - m['d1']) < 0:
                m['w0'] = m['w1'] = 0.5 * (m['w0'] + m['w1'])
    # section sizes in classes (the detail level sets how many; 0: every end its own): end widths and end depths
    # of all members, length-weighted; half-axes rounded to the level's digits, so that the linked global
    # variables (one per class) keep the values exactly
    D = DETAIL[detail]
    L = np.array([m_.get('length', 1.0) for m_ in R['members']])
    for key, k in (('w', D['widths']), ('d', D['depths'])):
        x = np.array([m_[key + e] for m_ in R['members'] for e in '01'])
        x = swbars.classes(x, np.repeat(L, 2), k)
        for i_, m_ in enumerate(R['members']):
            for e, v in zip('01', x[2 * i_: 2 * i_ + 2]):
                m_[key + e] = 2 * swbars._sig(v / 2, D['digits'])
    for m in R['members']:
        for end, j in ((0, m['i']), (1, m['j'])):            # the hub caps again (classes may exceed them)
            if j in R['hub_of'] and not (m['i'] in R['hub_of'] and m['j'] in R['hub_of']):
                k = 'w0' if end == 0 else 'w1'
                m[k] = min(m[k], 2 * swbars._sig(HUB_END * hw['hubs'][R['hub_of'][j]]['ro'] / 2, D['digits']))
        if (m['w0'] - m['d0']) * (m['w1'] - m['d1']) < 0:     # major axis on one side along the member
            m['d0'] = m['d1'] = max(m['d0'], m['d1'])
            if (m['w0'] - m['d0']) * (m['w1'] - m['d1']) < 0:
                m['w0'] = m['w1'] = max(m['w0'], m['w1'])
    R['hubs'], R['shadow'], R['axis'] = hw['hubs'], hw['shadow'], hw['axis']
    # members on a hub reach into it until their end section is inside the hub's disk (a wide end section on
    # the ring stuck out past the hub's round face in SolidWorks); what reaches the hole is cut with it
    for m in R['members']:
        for end, j in ((0, m['i']), (1, m['j'])):
            if j in R['hub_of']:
                hub, k = R['hubs'][R['hub_of'][j]], 'w0' if end == 0 else 'w1'
                m['path'] = into_hub(m['path'], end, hub, m[k])
    # members stop short of a joint's centre (30 % of the joint node's radius), inside the node:
    # lofts ending at one point in different planes gave invalid solids and hung OCC's Booleans
    rj = {j: r for j, r, d in joints_of(R)}
    for m in R['members']:
        for end, j in ((0, m['i']), (1, m['j'])):
            if j in rj:
                m['path'] = trim(m['path'], end, 0.3 * rj[j])
    return R


def into_hub(path, end, hub, w, margin=0.3):
    """Path extended straight along its end direction into the hub until the end section's sides (w/2 either
    side of the path) are inside the hub's disk by margin; at most to the hub's centre."""
    t = tangent(path, end) * (1 if end == 1 else -1)          # outward at this end, i.e. into the hub
    p = np.array(path[-1][-1] if end == 1 else path[0][1], float)
    c, ro = np.asarray(hub['c'], float), hub['ro']
    n = np.array([-t[1], t[0]])
    out = lambda q: max(np.linalg.norm(q + k * 0.5 * w * n - c) for k in (-1, 1)) - (ro - margin)
    # the path end stays inside the disk: s up to the chord along t; the first s with the section inside, else
    # the s where it sticks out least (a path meeting the hub at a slant went through and out the far side)
    b = (p - c) @ t
    chord = -b + math.sqrt(max(0.0, b * b - (np.linalg.norm(p - c) ** 2 - (ro - margin) ** 2)))
    # and short of the hole: the hole's cut then meets only the hub (a member reaching into the hole made the
    # cut's volume depend on the lofts' small differences)
    S = [x for x in np.arange(0.0, max(chord, 0.0), 0.1)
         if np.linalg.norm(p + x * t - c) >= hub['ri'] + margin
         and min(np.linalg.norm(p + x * t + k * 0.5 * w * n - c) for k in (-1, 1)) >= hub['ri'] + margin]
    if len(S) == 0:
        return path
    S = np.array(S)
    o = np.array([out(p + x * t) for x in S])
    s = float(S[np.argmax(o <= 0)]) if np.any(o <= 0) else float(S[np.argmin(o)])
    if s <= 0:
        return path
    q = tuple(map(float, p + s * t))
    path = list(path)
    if end == 1:
        e = path[-1]
        path[-1] = ('line', e[1], q) if e[0] == 'line' else e
        if e[0] != 'line':
            path.append(('line', tuple(map(float, p)), q))
    else:
        e = path[0]
        path[0] = ('line', q, e[2]) if e[0] == 'line' else e
        if e[0] != 'line':
            path.insert(0, ('line', q, tuple(map(float, p))))
    return path


def edge_length(e):
    if e[0] == 'line':
        return float(np.linalg.norm(np.subtract(e[2], e[1])))
    c = np.array(members.S_circle(*e[1:])); r = np.linalg.norm(np.array(e[1]) - c)
    ang = lambda v: math.atan2(v[1] - c[1], v[0] - c[0])
    a0 = ang(e[1]); sw = (ang(e[3]) - a0) % (2 * math.pi)
    if not ((ang(e[2]) - a0) % (2 * math.pi) < sw):
        sw = 2 * math.pi - sw
    return float(r * sw)


def trim(path, end, s):
    """Path with length s taken off its start (end 0) or end (end 1): end edges shorter than what is left to
    take off are dropped (an unchanged path ended at the joint's centre with the other members, and OCC's
    unions of such lofts lost material)."""
    path = list(path)
    k = 0 if end == 0 else -1
    while len(path) > 1 and edge_length(path[k]) < s + 0.5:
        s -= edge_length(path[k])
        path.pop(k)
    if s <= 0:
        return path
    return _trim_edge(path, end, s)


def _trim_edge(path, end, s):
    path = list(path)
    k = 0 if end == 0 else -1
    e = path[k]
    if e[0] == 'line':
        p, q = np.array(e[1]), np.array(e[2])
        L = np.linalg.norm(q - p)
        if L < s + 0.5:
            return path
        u = (q - p) / L
        e = ('line', tuple(p + s * u), e[2]) if end == 0 else ('line', e[1], tuple(q - s * u))
    else:
        a, m_, b = (np.array(x) for x in e[1:])
        c = np.array(members.S_circle(*e[1:])); r = np.linalg.norm(a - c)
        ang = lambda v: math.atan2(v[1] - c[1], v[0] - c[0])
        a0, a1 = ang(a), ang(b)
        sw = (a1 - a0) % (2 * math.pi)
        if not ((ang(m_) - a0) % (2 * math.pi) < sw):
            sw -= 2 * math.pi
        if abs(sw) * r < s + 0.5:
            return path
        ds = math.copysign(s / r, sw)
        if end == 0:
            a0 += ds
        else:
            a1 = a0 + sw - ds
        sw2 = (a1 - a0) if end == 1 else (a0 + (sw - ds) - a0)
        p0 = c + r * np.array([math.cos(a0), math.sin(a0)])
        p1 = c + r * np.array([math.cos(a0 + sw2), math.sin(a0 + sw2)])
        pm = c + r * np.array([math.cos(a0 + sw2 / 2), math.sin(a0 + sw2 / 2)])
        e = ('arc', tuple(p0), tuple(pm), tuple(p1))
    path[k] = tuple([e[0]] + [tuple(map(float, x)) for x in e[1:]])
    return path


def _wire3(path, z):
    es = []
    for e in path:
        P = [Vector(p[0], p[1], z) for p in e[1:]]
        es.append(cq.Edge.makeLine(*P) if e[0] == 'line' else cq.Edge.makeThreePointArc(*P))
    return cq.Wire.assembleEdges(es)


def member_solid(m, z):
    prof = []
    for c, t, pa, pb, a, b in sections(m, z):
        prof.append(cq.Wire.assembleEdges([cq.Edge.makeEllipse(a, b, Vector(*c), Vector(*t), Vector(*(pa - c)))]))
    return cq.Solid.sweep_multi(prof, _wire3(m['path'], z), True, False)


def joints_of(R):
    """Joint nodes not on a hub: (joint, horizontal semi-axis r, vertical semi-axis c) of an ellipsoid of
    revolution 5 % over the widest and deepest member end there. The members stop 0.3 r short of the centre,
    so their end sections lie inside it (0.3^2 + 1/1.05^2 < 1) and nothing touches tangentially (cylinders as
    deep as the members touched them at the top and gave invalid solids; lower ones left steps)."""
    ends = {}
    for m in R['members']:
        for j, w, d in ((m['i'], m['w0'], m['d0']), (m['j'], m['w1'], m['d1'])):
            r, dd = ends.get(j, (0, 0)); ends[j] = (max(r, w / 2), max(dd, d / 2))
    return [(j, 1.05 * r, 1.05 * c) for j, (r, c) in sorted(ends.items()) if j not in R['hub_of']]


def node_solid(p, r, c, z):
    """Ellipsoid of revolution about the vertical through p: half-ellipse (r across, c up) revolved."""
    C = Vector(p[0], p[1], z)
    arc = cq.Edge.makeEllipse(r, c, C, Vector(1, 0, 0), Vector(0, 1, 0), -90, 90) if r >= c else \
        cq.Edge.makeEllipse(c, r, C, Vector(1, 0, 0), Vector(0, 0, 1), 0, 180)
    ln = cq.Edge.makeLine(Vector(p[0], p[1], z + c), Vector(p[0], p[1], z - c))
    face = cq.Face.makeFromWires(cq.Wire.assembleEdges([arc, ln]), [])
    return cq.Solid.revolve(face, 360, Vector(p[0], p[1], z - c), Vector(p[0], p[1], z + c))


def steps(R):
    """Features in build order with the expected volume and centre of mass after each (CadQuery)."""
    z = R['mid']
    out, union = [], None

    def combine(solid, cut=False):
        """union with solid added (or cut): plain Boolean first, the fuzzy one (1e-4 mm) when the result is
        invalid (sections meeting a joint cylinder nearly tangentially); None if both are invalid."""
        if union is None:
            return solid
        n0, v0, vs = len(union.Solids()), union.Volume(), solid.Volume()
        for tol in (None, 1e-4, 1e-3):
            if cut:
                u = (union.cut(solid, tol=tol) if tol else union.cut(solid)).clean()
            else:
                u = (union.fuse(solid, tol=tol) if tol else union.fuse(solid)).clean()
            v = u.Volume()
            # a fuse keeps at least the larger part and at most both; a cut removes at most the tool
            e = 1e-4 * (v0 + vs) + 1e-2                   # integration error of lofted faces
            ok_v = (v0 - vs - e <= v <= v0 + e) if cut else (max(v0, vs) - e <= v <= v0 + vs + e)
            if u.isValid() and len(u.Solids()) <= n0 + (0 if cut else 1) and ok_v:
                return u
        return None

    def add(kind, name, solid, cut=False, **kw):
        nonlocal union
        u = combine(solid, cut)
        if u is None:
            raise RuntimeError('%s: no valid Boolean' % name)
        union = u
        v, c = mass_props(union)
        out.append(dict(kind=kind, name=name, volume=v, com=np.array(c), **kw))

    for n, h in enumerate(R['hubs'], 1):
        s = cq.Workplane().workplane(offset=h['lo']).center(*h['c']).circle(h['ro']).extrude(h['hi'] - h['lo']).val()
        add('extrude', 'Hub %d' % n, s, op='add', lo=h['lo'], hi=h['hi'], loops=[[('circle', tuple(map(float, h['c'])), float(h['ro']))]])
    for n, m in enumerate(R['members'], 1):
        add('member', 'Member %d' % n, member_solid(m, z), m=m, own=member_solid(m, z).Volume())
    J = R['joints']
    for n, (j, r, c) in enumerate(joints_of(R), 1):
        # a node that OCC cannot fuse (near-tangent contact with a member: classes make sizes equal) is made
        # slightly larger
        for f in (1.0, 1.03, 1.06, 1.1):
            sol = node_solid(J[j], f * r, f * c, z)
            if combine(sol) is not None:
                break
        add('node', 'Joint %d' % n, sol, centre=(float(J[j][0]), float(J[j][1]), z), r=f * r, c=f * c, own=sol.Volume())
    for n, h in enumerate(R['hubs'], 1):
        s = cq.Workplane().workplane(offset=h['lo'] - EPS).center(*h['c']).circle(h['ri']).extrude(h['hi'] - h['lo'] + 2 * EPS).val()
        add('extrude', 'Hole %d' % n, s, cut=True, op='cut', lo=h['lo'], hi=h['hi'], loops=[[('circle', tuple(map(float, h['c'])), float(h['ri']))]])
    return out, union


def write(path, name, R, detail=None, nd=6):
    S, union = steps(R)
    z = R['mid']
    R0, back = swgen.BUILD_FRAME[2]
    size = union.BoundingBox().DiagonalLength
    t3 = lambda q: '(%.6f, %.6f, %.6f)' % tuple(float(x) for x in q)
    P = lambda q: '(%.6f, %.6f)' % tuple(float(x) for x in q)
    allp = np.array([q for s in S if s['kind'] == 'extrude' for lp in s['loops'] for e in lp for q in e[1:2]] +
                    [s['centre'][:2] for s in S if s['kind'] == 'node'])
    L = ['"""%s: SolidWorks feature tree for a topology-optimized part, generated by stl2cad (TO swept members).' % name, '',
         'The part is built as %d hub(s), %d member(s) and %d joint node(s) (revolved half-ellipses). A member is a Loft of two ellipse sections'
         % (len(R['hubs']), len(R['members']), len(joints_of(R))),
         '(normal to its centre line at the two ends) along its centre line (lines and tangent arcs on the',
         'mid-plane at %.4g mm), so its section changes linearly along it. The original\'s holes are cut last.' % z,
         'After every feature the volume and centre of mass are checked against stl2cad\'s values.', '',
         'Needs: Windows, SolidWorks, Python with pywin32 (pip install pywin32).',
         'Run:   python %s.py [C:\\\\path\\\\part.SLDPRT]' % name, '"""', '',
         'R_FILLET = 0.0', 'ROUND_DIGITS = %d                # significant digits of the dimensions, set at the end (0: as fitted)' % DETAIL[detail or 3]['digits'],
         'DEFER_ROUNDING = True', 'DIM_DECIMALS = 2',
         'DIMENSION = True                # relations and dimensions in the member sketches (False: geometry only)',
         '# detail level %s (1: fewest dimensions ... 5: smallest model error)' % detail, '',
         'AXIS_I, AXIS_J = 0, 1', 'PLANE_NAME = "Front Plane"', 'PLANE_INDEX = 0',
         'ROT = %s' % [list(map(float, r)) for r in R0], 'ROT_BACK = %s' % (tuple(back),),
         'FINAL_CENTROID = %s' % t3(S[-1]['com']), 'SIZE = %s' % swgen._f(size, 3),
         'HELPER_AB = (%s, %s)' % (swgen._f(allp[:, 0].min() - 10, 3), swgen._f(allp[:, 1].min() - 10, 3)),
         'MID = %.6f                # mm, mid-plane of the members (along z)' % z, '',
         '# ---- features in build order. extrude: as stl2cad\'s FEATURES; member: path (lines and arcs, mm, on the',
         '# mid-plane), sections (centre, major-axis point, minor-axis point, path direction; model mm), expected volume ----',
         'STEPS = [']
    for s in S:
        cb = np.asarray(R0, float) @ s['com']
        if s['kind'] == 'node':
            L += ['    dict(kind="node", name="%s", volume=%.4f, centroid=%s, own=%.4f, centre=%s, r=%.6f, c=%.6f),'
                  % (s['name'], s['volume'], t3(cb), s['own'], t3(s['centre']), s['r'], s['c'])]
            continue
        if s['kind'] == 'extrude':
            L += ['    dict(kind="extrude", name="%s", op="%s", lo=%s, hi=%s, volume=%.4f, centroid=%s, holes=False,'
                  % (s['name'], s['op'], swgen._f(s['lo'], nd), swgen._f(s['hi'], nd), s['volume'], t3(cb)),
                  '         expect="fully defined", plans=[],',
                  '         loops=[[(\'circle\', %s, %s)]]),' % (P(s['loops'][0][0][1]), swgen._f(s['loops'][0][0][2], nd))]
        else:
            m = s['m']
            secs = sections(m, z)
            L += ['    dict(kind="member", name="%s", volume=%.4f, centroid=%s, own=%.4f,' % (s['name'], s['volume'], t3(cb), s['own']),
                  '         w=(%.4f, %.4f), d=(%.4f, %.4f),' % (m['w0'], m['w1'], m['d0'], m['d1']),
                  '         path=[']
            for e in m['path']:
                L.append('             (%r, %s),' % (e[0], ', '.join(P(q) for q in e[1:])))
            L += ['         ],', '         sections=[']
            for c, t, pa, pb, a, b in secs:
                L.append('             (%s, %s, %s, %s),' % (t3(c), t3(pa), t3(pb), t3(t)))
            L += ['         ]),']
    L += [']', 'FEATURES = []', 'CYLINDERS = []', 'REVOLVES = []', 'REV_TURN = None', 'REV_DONE = None', 'CUTS = []',
          'XCUTS = []', 'BLENDS = []', 'JUNCTIONS = []', '']
    head, _ = swgen._RUNTIME.split('if __name__ == "__main__":')
    L.append(head + SWEEP_CODE)
    open(path, 'w').write('\n'.join(L) + '\n')
    return dict(volume=S[-1]['volume'], union=union, steps=len(S))


def verify(script_path, ref, timeout=900):
    """Run the script on the API stand-in (swmock) and compare the result with ref by volume, centre of mass
    and bounding box (swplates.verify's Booleans and Gauss-Kronrod volumes do not finish on lofted faces)."""
    import re, subprocess, tempfile
    import stl2cad
    d = tempfile.mkdtemp(prefix='swsweep_')
    py, out = os.path.join(d, 'part.py'), os.path.join(d, 'part.step')
    open(py, 'w').write(re.sub(r'^R_FILLET = [0-9.]+', 'R_FILLET = 0.0', open(script_path).read(), flags=re.M))
    pkg = os.path.dirname(stl2cad.__file__)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.path.join(pkg, 'swmock'), os.path.dirname(pkg)]),
               SWMOCK_CUT_REVERSED='1', SWMOCK_BOSS_REVERSED='0', SWMOCK_EDGE_SIGHT='1')
    r = subprocess.run([sys.executable, py, out], capture_output=True, text=True, timeout=timeout, env=env)
    log = r.stdout + r.stderr
    if r.returncode != 0 or not os.path.exists(out):
        return False, log
    s = cq.Shape.importBrep(out + '.brep') if os.path.exists(out + '.brep') else cq.importers.importStep(out).val()
    (v, c), (vr, cr) = mass_props(s), mass_props(ref)
    b, br = s.BoundingBox(), ref.BoundingBox()
    db = max(abs(x - y) for x, y in zip((b.xmin, b.ymin, b.zmin, b.xmax, b.ymax, b.zmax),
                                         (br.xmin, br.ymin, br.zmin, br.xmax, br.ymax, br.zmax)))
    dc = math.dist(c, cr)
    # 0.1 %: Booleans on lofted faces differ by a few mm^3 between runs of OCC; boxes of B-spline faces are loose
    ok = abs(v - vr) <= 1e-3 * vr and dc <= 0.05 and db <= 0.5
    return ok, log + '\nstand-in: volume %.3f, reference %.3f mm^3; centres of mass %.5f mm apart; bounding boxes %.5f mm' % (v, vr, dc, db)


SWEEP_CODE = r'''
# ---- swept members ----------------------------------------------------------------------------------------
def to_sk(T, p):
    """Model point (mm) -> sketch coordinates (m), by the sketch's ModelToSketchTransform (ArrayData: rotation
    R row-major applied as p R, translation, scale)."""
    x = [c * MM for c in p]
    return tuple((sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3))


def open_sketch(sw, plane):
    """Sketch on a plane feature; its ModelToSketchTransform as a list."""
    sw.model.ClearSelection2(True)
    plane.Select2(False, 0)
    sm = sw.model.SketchManager
    sm.InsertSketch(True)
    settle()
    sk = call0(sm, "ActiveSketch")
    if sk is None:
        sk = call0(sw.model, "GetActiveSketch2")
    T = list(sk.ModelToSketchTransform.ArrayData)
    sm.AddToDB = True
    helper = sm.CreateCircleByRadius(-0.5, -0.5, 0.0, HELPER_R * MM)     # lines need another entity first
    return sm, T, helper


def close_sketch(sw, sm, helper):
    sm.AddToDB = False
    if helper is not None:
        sw.remove_helper(helper)
    sm.InsertSketch(True)
    sw.report_status()
    return sw.model.FeatureByPositionReverse(0)


def sketch_points(sm):
    """The open sketch's points as (x, y, object), sketch coordinates in m (GetSketchPoints2: a line's, arc's or
    ellipse's own point methods are not reachable over late binding, sw_probe_dims)."""
    sk = call0(sm, "ActiveSketch")
    return [(p.X, p.Y, p) for p in (call0(sk, "GetSketchPoints2") or [])]


def at(pts, q, tol=1e-7):
    return [o for x, y, o in pts if abs(x - q[0]) < tol and abs(y - q[1]) < tol]


def rel(sw, name, objs, what):
    """Relation between selected objects (SketchAddConstraints returns nothing: counted as added)."""
    try:
        sw.model.ClearSelection2(True)
        ok = all(o.Select4(k > 0, NOTHING) for k, o in enumerate(objs))
        if not ok:
            sw.rels_failed.append(what + " (selection)"); return
        sw.model.SketchAddConstraints(name)
        sw.rels_ok += 1
    except Exception as ex:
        sw.rels_failed.append(f"{what}: {type(ex).__name__}: {ex}")


def dim(sw, objs, kind, text, what, origin=False):
    """Dimension of selected objects (and the origin): kind H, V or any; text position in sketch m."""
    try:
        sw.model.ClearSelection2(True)
        ok = all(o.Select4(k > 0, NOTHING) for k, o in enumerate(objs))
        if origin:
            ok = ok and sw.model.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0, True, 0, NOTHING, 0)
        d = None
        if ok:
            d = {"any": sw.model.AddDimension2, "H": sw.model.AddHorizontalDimension2,
                 "V": sw.model.AddVerticalDimension2}[kind](text[0], text[1], 0.0)
        sw._count(d, what)
    except Exception as ex:
        sw.dims_failed.append(f"{what}: {type(ex).__name__}: {ex}")


def last_feature(sw):
    return sw.model.FeatureByPositionReverse(0)


def offset_plane(sw, base, dist, point, name):
    """Reference plane at distance dist (mm) from the base plane (by name), through the model point (mm):
    the side is checked through a sketch's transform (the other side tried, flag 256, if the first is wrong)."""
    for flip in (0, 256):                     # swRefPlaneReferenceConstraint_OptionFlip
        sw.model.ClearSelection2(True)
        if not sw.model.Extension.SelectByID2(base, "PLANE", 0, 0, 0, False, 0, NOTHING, 0):
            f = sw.model.FeatureByName(base)
            if f is None or not f.Select2(False, 0):
                raise RuntimeError(f"could not select {base}")
        p = sw.model.FeatureManager.InsertRefPlane(8 | flip, abs(dist) * MM, 0, 0, 0, 0)
        if p is None:
            continue
        f = last_feature(sw)
        sm, T, h = open_sketch(sw, f)
        z = to_sk(T, point)[2]
        sm.AddToDB = False
        sm.InsertSketch(True)
        sk = last_feature(sw)
        sw.model.ClearSelection2(True); sk.Select2(False, 0); sw.model.Extension.DeleteSelection2(0)
        if abs(z) < 1e-6:
            f.Name = name
            return f
        print(f"  {name}, try {flip}: plane on the wrong side ({z / MM:.4f} mm off)")
        sw.model.ClearSelection2(True); f.Select2(False, 0); sw.model.Extension.DeleteSelection2(0)
    raise RuntimeError(f"could not make {name}")


def mid_plane(sw):
    f = offset_plane(sw, PLANE_NAME, MID, (0.0, 0.0, MID), "Mid-plane")
    print(f"mid-plane at z = {MID:.4f} mm")
    return f


def node_step(sw, st, drift, made):
    """Joint node: a half-ellipse (r across, c up) on a plane through the joint parallel to the Right Plane,
    revolved about the vertical through the joint (FeatureRevolve2, as stl2cad's revolves)."""
    x, y, z = st["centre"]
    r, c = st["r"], st["c"]
    if abs(x) > 1e-6:
        plane = offset_plane(sw, "Right Plane", x, (x, y, z), st["name"] + " plane")
    else:
        plane = sw.model.FeatureByName("Right Plane")
    made["planes"].append(plane)
    sm, T, h = open_sketch(sw, plane)
    Q = lambda p: to_sk(T, p)
    C, top, bot = Q((x, y, z)), Q((x, y, z + c)), Q((x, y, z - c))
    side = Q((x, y + r, z))
    if r >= c:
        A, B = side, top
    else:
        A, B = top, side
    turn = arc_turn(top[:2], side[:2], bot[:2])
    make(sm.CreateEllipticalArc, C[0], C[1], 0.0, A[0], A[1], 0.0, B[0], B[1], 0.0,
         top[0], top[1], 0.0, bot[0], bot[1], 0.0, turn)
    make(sm.CreateLine, bot[0], bot[1], 0.0, top[0], top[1], 0.0)
    e0, e1 = Q((x, y, z - c - 1.0)), Q((x, y, z + c + 1.0))
    cl = make(sm.CreateCenterLine, e0[0], e0[1], 0.0, e1[0], e1[1], 0.0)
    sk = close_sketch(sw, sm, h)
    made["sketches"].append(sk)
    fm = sw.model.FeatureManager
    v, cc = expected(st, drift)
    for mark in (None, 4, 16):
        sw.model.ClearSelection2(True)
        sk.Select2(False, 0)
        if mark is not None:
            sd = call0(sw.model.SelectionManager, "CreateSelectData"); sd.Mark = mark
            cl.Select4(True, sd)
        f = fm.FeatureRevolve2(True, True, False, False, False, False, 0, 0, 2 * math.pi, 0.0,
                               False, False, 0.0, 0.0, 0, 0.0, 0.0, True, True, True)
        if f is None:
            print(f"  {st['name']}: revolve{'' if mark is None else f' (centerline Mark {mark})'} returned no feature")
            continue
        vv, c2 = sw.mass()
        # 5 % of the node's own volume: it overlaps the lofts, whose small differences from CadQuery change
        # the overlap (SolidWorks: 0.15 mm^3 on a 30 mm^3 node; OCC unions on lofts: up to 5 mm^3 on 240 mm^3)
        if abs(vv - v) <= 0.05 * st["own"] + 1e-3 and math.dist(c2, cc) <= 1e-4 * SIZE + 1e-2:
            f.Name = st["name"]
            print(f"{st['name']}: node, volume {vv:.3f} mm^3 (expected {v:.3f})")
            carry(drift, vv, c2, v, cc)
            return
        print(f"  {st['name']}: revolve gave volume {vv:.3f} mm^3, expected {v:.3f}")
        sw.model.ClearSelection2(True); f.Select2(False, 0); sw.model.Extension.DeleteSelection2(0)
    raise RuntimeError(f"{st['name']}: no revolve gives the expected node")


def hide(sw, made):
    """Planes and member sketches hidden (they stay in the tree)."""
    n = 0
    for kind, call in (("planes", "BlankRefGeom"), ("sketches", "BlankSketch")):
        for f in made[kind]:
            try:
                sw.model.ClearSelection2(True)
                f.Select2(False, 0)
                getattr(sw.model, call)()
                n += 1
            except Exception:
                pass
    sw.model.ClearSelection2(True)
    print(f"hid {n} plane(s) and sketch(es)")


def arc_turn(s, m, t):
    return 1 if (m[0] - s[0]) * (t[1] - s[1]) - (m[1] - s[1]) * (t[0] - s[0]) > 0 else -1


def path_sketch(sw, plane, path, name=""):
    """Centre line on the mid-plane; also a sketch point at each end. Returns the sketch feature, the first and
    last segment and the two end points (COM objects, selected later by object: a pick by coordinates took a
    nearby point of another member's sketch, 1.4 mm away, in SolidWorks)."""
    sm, T, h = open_sketch(sw, plane)
    Q = lambda q: to_sk(T, (q[0], q[1], MID))[:2]
    segs = []
    for e in path:
        if e[0] == "line":
            segs.append(make(sm.CreateLine, *Q(e[1]), 0.0, *Q(e[2]), 0.0))
        else:
            s, m, t = (Q(q) for q in e[1:])
            xc, yc = circle_centre(s, m, t)
            segs.append(make(sm.CreateArc, xc, yc, 0.0, s[0], s[1], 0.0, t[0], t[1], 0.0, arc_turn(s, m, t)))
    ends = [Q(path[0][1]), Q(path[-1][-1])]
    pts = [sm.CreatePoint(*ends[0], 0.0), sm.CreatePoint(*ends[1], 0.0)]
    sm.AddToDB = False
    if DIMENSION:
        dimension_path(sw, sm, path, segs, ends, Q, name)
    f = close_sketch(sw, sm, h)
    return f, (segs[0], segs[-1]), pts


def dimension_path(sw, sm, path, segs, ends, Q, name):
    """Path sketch fully defined: line-arc tangencies; the extra end points coincident with the path's ends;
    the ends and each arc's centre by horizontal and vertical distance from the origin; each arc's radius."""
    P = sketch_points(sm)
    for k in range(len(segs) - 1):
        if (path[k][0] == "arc") != (path[k + 1][0] == "arc"):
            rel(sw, "sgTANGENT", [segs[k], segs[k + 1]], f"{name}: tangency {k + 1}")
    for e, q in zip(("start", "end"), ends):
        two = at(P, q)
        if len(two) >= 2:
            rel(sw, "sgCOINCIDENT", two[:2], f"{name}: {e} point")
    def xy(q, what):
        o = at(P, q)
        if not o:
            sw.dims_failed.append(f"{what}: no point"); return
        if abs(q[0]) > 1e-7:
            dim(sw, [o[0]], "H", (q[0] / 2, q[1] - 0.004), f"{what} x", origin=True)
        if abs(q[1]) > 1e-7:
            dim(sw, [o[0]], "V", (q[0] - 0.004, q[1] / 2), f"{what} y", origin=True)
    xy(ends[0], f"{name}: start")
    xy(ends[1], f"{name}: end")
    for k, e in enumerate(path):
        if e[0] != "arc":
            continue
        s_, m_, t_ = (Q(q) for q in e[1:])
        c = circle_centre(s_, m_, t_)
        xy(c, f"{name}: arc {k + 1} centre")
        txt = (m_[0] + (m_[0] - c[0]) * 0.5, m_[1] + (m_[1] - c[1]) * 0.5)
        dim(sw, [segs[k]], "any", txt, f"{name}: radius of path arc {k + 1}")


def plane_through(sw, plane, p, t):
    """Whether the plane feature holds the model point p (mm) and is normal to the direction t: a sketch on it
    is opened, the transform read, the sketch deleted."""
    sm, T, h = open_sketch(sw, plane)
    z0 = to_sk(T, p)[2]
    z1 = to_sk(T, tuple(p[k] + t[k] for k in range(3)))[2]
    sm.AddToDB = False
    sm.InsertSketch(True)
    sk = last_feature(sw)
    sw.model.ClearSelection2(True); sk.Select2(False, 0); sw.model.Extension.DeleteSelection2(0)
    return abs(z0) < 1e-6 and abs(abs(z1 - z0) - MM) < 1e-6, z0 / MM


def end_plane(sw, path, end, seg, pt, tan):
    """Plane normal to the path at its start (end 0) or end (end 1): perpendicular to the end segment (Mark 0),
    coincident with the end point (Mark 1). The segment and point are selected as objects; if the plane is not
    through the end point, by coordinates; the plane is checked either way."""
    e = path[0] if end == 0 else path[-1]
    p = e[1] if end == 0 else e[-1]
    P = (p[0], p[1], MID)
    mid = e[2] if e[0] == "arc" else tuple(0.5 * (e[1][k] + e[2][k]) for k in range(2))
    ext, fm, model = sw.model.Extension, sw.model.FeatureManager, sw.model

    def by_object():
        sd0 = call0(model.SelectionManager, "CreateSelectData"); sd0.Mark = 0
        sd1 = call0(model.SelectionManager, "CreateSelectData"); sd1.Mark = 1
        return seg.Select4(False, sd0) and pt.Select4(True, sd1)

    def by_coordinates():
        return (ext.SelectByID2("", "EXTSKETCHSEGMENT", mid[0] * MM, mid[1] * MM, MID * MM, False, 0, NOTHING, 0)
                and ext.SelectByID2("", "EXTSKETCHPOINT", p[0] * MM, p[1] * MM, MID * MM, True, 1, NOTHING, 0))

    tried = []
    for how, sel in (("objects", by_object), ("coordinates", by_coordinates)):
        model.ClearSelection2(True)
        try:
            ok = sel()
        except Exception as ex:
            tried.append(f"{how}: {type(ex).__name__}"); continue
        if not ok or fm.InsertRefPlane(2, 0, 4, 0, 0, 0) is None:      # perpendicular, coincident
            tried.append(f"{how}: {'no selection' if not ok else 'no plane'}"); continue
        f = last_feature(sw)
        good, off = plane_through(sw, f, P, tan)
        if good:
            if tried:
                print(f"  path end {end}: plane by {how} ({'; '.join(tried)})")
            return f
        tried.append(f"{how}: plane {off:+.4f} mm off the end point")
        model.ClearSelection2(True); f.Select2(False, 0); ext.DeleteSelection2(0)
    raise RuntimeError(f"path end {end}: no plane normal to the path through its end ({'; '.join(tried)})")


def section_sketch(sw, plane, sec, seg=None, name=""):
    c, pa, pb = sec[:3]
    sm, T, h = open_sketch(sw, plane)
    C, A, B = (to_sk(T, q) for q in (c, pa, pb))
    if max(abs(C[2]), abs(A[2]), abs(B[2])) > 1e-6:
        raise RuntimeError(f"section not on its plane (sketch z {C[2] / MM:.4f}, {A[2] / MM:.4f}, {B[2] / MM:.4f} mm)")
    make(sm.CreateEllipse, C[0], C[1], 0.0, A[0], A[1], 0.0, B[0], B[1], 0.0)
    sm.AddToDB = False
    if DIMENSION and seg is not None:
        # fully defined (sw_probe_dims): centre pierced by the path's end segment, major axis along a sketch axis,
        # the two half-axes (labelled radius: equal ones are linked across sketches)
        P = sketch_points(sm)
        oc, oa, ob = at(P, C), at(P, A), at(P, B)
        if oc and oa and ob:
            rel(sw, "sgATPIERCE", [oc[0], seg], f"{name}: centre on the path")
            if abs(A[1] - C[1]) < 1e-9:
                rel(sw, "sgHORIZONTALPOINTS2D", [oa[0], oc[0]], f"{name}: major axis")
            elif abs(A[0] - C[0]) < 1e-9:
                rel(sw, "sgVERTICALPOINTS2D", [oa[0], oc[0]], f"{name}: major axis")
            for o, q, what in ((oa, A, "major"), (ob, B, "minor")):
                txt = (q[0] + 0.5 * (q[0] - C[0]) + 0.001, q[1] + 0.5 * (q[1] - C[1]) + 0.001)
                dim(sw, [oc[0], o[0]], "any", txt, f"{name}: radius of section ({what} half-axis)")
        else:
            sw.dims_failed.append(f"{name}: ellipse points not found")
    return close_sketch(sw, sm, h)


def member_step(sw, plane, st, drift, made):
    v0, c0 = sw.mass()
    path, segs, pts = path_sketch(sw, plane, st["path"], st["name"] + " path")
    tans = [tuple(st["sections"][k][3]) for k in (0, 1)]
    pl = [end_plane(sw, st["path"], k, segs[k], pts[k], tans[k]) for k in (0, 1)]
    secs = [section_sketch(sw, pl[k], st["sections"][k], segs[k], f"{st['name']} {('start', 'end')[k]} section")
            for k in (0, 1)]
    made["planes"] += pl
    made["sketches"] += [path] + secs
    m = sw.model
    m.ClearSelection2(True)
    secs[0].Select2(False, 1); secs[1].Select2(True, 1); path.Select2(True, 4)
    f = m.FeatureManager.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False,
                                                0.0, 0.0, 0, True, True, True, 0)
    if f is None:
        raise RuntimeError(f"{st['name']}: no loft")
    f.Name = st["name"]
    for k, p in enumerate(pl):
        p.Name = f"{st['name']} plane {'start' if k == 0 else 'end'}"
    v, c = sw.mass()
    ev, ec = expected(st, drift)
    added = v - v0
    if abs(v - ev) > 0.01 * st["own"] + 1e-3:
        raise RuntimeError(f"{st['name']}: volume {v:.3f} mm^3, expected {ev:.3f} (member alone {st['own']:.3f})")
    # the loft's small difference from CadQuery's sweep is carried into later checks
    carry(drift, v, c, ev, ec)
    print(f"{st['name']}: loft, volume {v:.3f} mm^3 (expected {ev:.3f}; difference {v - ev:+.3f})")


def carry(drift, v, c, ev, ec):
    """Add a step's difference from the expected volume and first moment to the drift."""
    drift[0] += v - ev
    drift[1] = [drift[1][k] + v * c[k] - ev * ec[k] for k in range(3)]


def expected(st, drift):
    v = st["volume"] + drift[0]
    c = tuple((st["volume"] * st["centroid"][k] + drift[1][k]) / v for k in range(3))
    return v, c


def build_sweep():
    sw = SW()
    print(f"connected to SolidWorks; new part created; {len(STEPS)} step(s)")
    set_units(sw)
    drift = [0.0, [0.0, 0.0, 0.0]]
    plane = None
    made = dict(planes=[], sketches=[])
    for st in STEPS:
        if st["kind"] == "member":
            if plane is None:
                plane = mid_plane(sw)
                made["planes"].append(plane)
            member_step(sw, plane, st, drift, made)
        elif st["kind"] == "node":
            node_step(sw, st, drift, made)
        else:
            v, c = expected(st, drift)
            build_step(sw, dict(st, volume=v, centroid=c))
            carry(drift, *sw.mass(), v, c)
    v0 = sw.mass()[0]
    print(sw.link_dimensions())
    v1 = sw.mass()[0]
    if abs(v1 - v0) > 1e-4 * v0 + 1e-3:
        print(f"NOTE: linking changed the volume from {v0:.3f} to {v1:.3f} mm^3; the equations are removed")
        sw.unlink_dimensions()
    print(sw.round_dimensions())
    hide(sw, made)
    linked = getattr(sw, "linked", [])
    nvar = len({val for _, val in linked})
    print(f"independent dimensions: {sw.dims_ok - len(linked) + nvar} "
          f"({sw.dims_ok} dimensions, {len(linked)} of them driven by {nvar} global variable(s))")
    print(f"dimensions: {sw.dims_ok} added" + (f", {len(sw.dims_failed)} failed: {sw.dims_failed}" if sw.dims_failed else ""))
    print(f"relations: {sw.rels_ok} added" + (f", {len(sw.rels_failed)} failed: {sw.rels_failed}" if sw.rels_failed else ""))
    st_ = sw.status
    print("sketch status: " + ", ".join(f"{k} x {st_.count(k)}" for k in sorted(set(st_))))
    sw.restore_preferences()
    return sw


if __name__ == "__main__":
    sw = build_sweep()
    if len(sys.argv) > 1:
        sw.model.SaveAs3(sys.argv[1], 0, 0)
'''


if __name__ == '__main__':
    import time
    ap = argparse.ArgumentParser()
    ap.add_argument('base'); ap.add_argument('to'); ap.add_argument('name')
    ap.add_argument('--units', default='in'); ap.add_argument('--detail', type=int, default=3, choices=range(1, 6))
    ap.add_argument('--plot', action='store_true'); ap.add_argument('--no-verify', action='store_true')
    a = ap.parse_args()
    t0 = time.time()
    ns = swplates._namespace(a.base)
    to = trimesh.load(a.to); to.apply_scale({'m': 1000.0, 'mm': 1.0, 'in': 25.4}[a.units])
    to = max(to.split(only_watertight=False), key=lambda b: len(b.faces))
    R = fit(ns['original'], to, a.detail)
    info = write(a.name + '_solidworks.py', a.name + '_solidworks', R, a.detail)
    s = info['union']
    cq.exporters.export(cq.Workplane().add(s), a.name + '.step')
    cq.exporters.export(cq.Workplane().add(s), a.name + '.stl', tolerance=0.02, angularTolerance=0.2)
    m = trimesh.load(a.name + '.stl')
    import sweepref
    narc = sum(e[0] == 'arc' for mm in R['members'] for e in mm['path'])
    print('%d members (%d arcs in the paths), %d joints; model valid %s, solids %d' % (
        len(R['members']), narc, len(joints_of(R)), s.isValid(), len(s.Solids())))
    print('volume: TO %.0f, model %.0f mm^3 (%+.1f %%)' % (to.volume, m.volume, 100 * (m.volume / to.volume - 1)))
    print('voxel overlap: %.0f %% of the TO in the model, %.0f %% of the model in the TO' % sweepref.compare(to, m))
    print('SolidWorks script: %s_solidworks.py (%d steps)' % (a.name, info['steps']))
    if not a.no_verify:
        ok, log = verify(a.name + '_solidworks.py', s)
        open(a.name + '_standin.log', 'w').write(log)
        print('API stand-in: %s; %s' % ('ok' if ok else 'FAILED', log.strip().splitlines()[-1]))
    print('time %.0f s' % (time.time() - t0))
