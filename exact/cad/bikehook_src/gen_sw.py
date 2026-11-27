"""Writes ../bikehook_solidworks.py from ../bikehook_exact.py (CadQuery): the feature data and the volume and centre
of mass after each feature, in SolidWorks' build order (sw_template.py).

python gen_sw.py      (from exact/cad/bikehook_src)
"""
import importlib.util, math, os, sys
import numpy as np, cadquery as cq
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, '..')
spec = importlib.util.spec_from_file_location('bh', os.path.join(OUT, 'bikehook_exact.py'))
bh = importlib.util.module_from_spec(spec); spec.loader.exec_module(bh)


def r3(p): return tuple(round(float(v), 6) for v in p)


def ring(V, z):
    tg = bh.fillet_polygon(V, bh.radii(V, bh.R_F)); ents = []
    for i in range(len(V)):
        t1, m, t2 = tg[i]
        ents.append(('A', r3((*t1, z)), r3((*m, z)), r3((*t2, z))))
        nx = tg[(i + 1) % len(V)][0]
        if np.linalg.norm(nx - t2) > 1e-6:
            ents.append(('L', r3((*t2, z)), r3((*nx, z))))
    return ents


def j_ents(top):
    a = math.radians(45); r = bh.R_J
    P = [(bh.C_R[0] + r, 0, top), (bh.C_R[0] + r, 0, bh.C_R[1])]
    e = [('L', P[0], P[1]),
         ('A', P[1], (bh.C_R[0] + r * math.cos(a), 0, bh.C_R[1] - r * math.sin(a)), (bh.C_R[0], 0, bh.Z_FLAT)),
         ('L', (bh.C_R[0], 0, bh.Z_FLAT), (bh.C_L[0], 0, bh.Z_FLAT)),
         ('A', (bh.C_L[0], 0, bh.Z_FLAT), (bh.C_L[0] - r * math.cos(a), 0, bh.C_L[1] - r * math.sin(a)),
          (bh.C_L[0] - r, 0, bh.C_L[1]))]
    return e


def data():
    D = {}
    path = j_ents(bh.Z_LEG) + [('L', (bh.C_L[0] - bh.R_J, 0, bh.C_L[1]), (bh.C_L[0] - bh.R_J, 0, bh.Z_LIP))]
    D['PATH'] = [tuple([e[0]] + [r3(p) for p in e[1:]]) for e in path]
    u0 = bh.W - bh.E / bh.SLOPE
    loc = [(-u0, 0), (u0, 0), (bh.W, bh.E), (bh.W, bh.T - bh.E), (u0, bh.T), (-u0, bh.T), (-bh.W, bh.T - bh.E), (-bh.W, bh.E)]
    rad = [bh.R_BLEND, bh.R_BLEND, bh.R_SIDE, bh.R_SIDE, bh.R_BLEND, bh.R_BLEND, bh.R_SIDE, bh.R_SIDE]
    tg = bh.fillet_polygon(loc, rad); X0 = bh.C_R[0] + bh.R_J
    to3 = lambda q: r3((X0 + q[1], q[0], bh.Z_LEG))
    prof = []
    for i in range(len(loc)):
        t1, m, t2 = tg[i]
        prof.append(('A', to3(t1), to3(m), to3(t2)))
        nx = tg[(i + 1) % len(loc)][0]
        if np.linalg.norm(nx - t2) > 1e-9:
            prof.append(('L', to3(t2), to3(nx)))
    D['PROFILE'] = prof
    rib = bh.rib().val()
    f = [f for f in rib.Faces() if f.geomType() == 'PLANE' and abs(f.Center().y - bh.RIB_W) < 1e-6][0]
    q = lambda v: r3((v.x, 0.0, v.z))
    D['RIB'] = [('L', q(e.startPoint()), q(e.endPoint())) if e.geomType() == 'LINE'
                else ('A', q(e.startPoint()), q(e.positionAt(0.5)), q(e.endPoint())) for e in f.outerWire().Edges()]
    D['RIB_PTS'] = [r3(p) for s in (1, -1) for p in ((0.0, s * bh.RIB_W, bh.RIB_Z),
                    (-4 + 20 * math.cos(math.radians(-120)), s * bh.RIB_W, bh.RIB_Z + 20 + 20 * math.sin(math.radians(-120))),
                    (4 + 27 * math.cos(math.radians(-70)), s * bh.RIB_W, bh.RIB_Z + 27 + 27 * math.sin(math.radians(-70))))]
    D['SECS'] = [(float(z), ring(np.array(V), z)) for z, V in bh.SECTIONS]
    arms = []
    for a in bh.ARMS:
        P = np.array(a['points'], float)
        f, d = bh.spline(P)
        secs = []
        for i, ra, rb, mj in a['sec']:
            t = np.asarray(f(d[i], 1), float).reshape(-1)[:3]; t /= np.linalg.norm(t)
            mj = np.array(mj, float); mj = mj - (mj @ t) * t; mj /= np.linalg.norm(mj)
            if ra < rb:
                ra, rb, mj = rb, ra, np.cross(t, mj)
            n = np.cross(t, mj)
            c = P[i]
            secs.append((i, r3(c), r3(c + ra * mj), r3(c + rb * n)))
        arms.append(dict(points=[r3(p) for p in P], sec=secs))
    D['ARMS'] = arms
    tyre = j_ents(bh.Z_LEG + 20) + [('L', (bh.C_L[0] - bh.R_J, 0, bh.C_L[1]), (bh.C_L[0] - bh.R_J, 0, bh.Z_LEG + 20)),
                                    ('L', (bh.C_L[0] - bh.R_J, 0, bh.Z_LEG + 20), (bh.C_R[0] + bh.R_J, 0, bh.Z_LEG + 20))]
    D['TYRE'] = [tuple([e[0]] + [r3(p) for p in e[1:]]) for e in tyre]
    D['TOP_BOX'] = [(-150.0, -150.0), (150.0, -150.0), (150.0, 150.0), (-150.0, 150.0)]
    zb = bh.H_BOSS + 25.0
    box = [('L', (a[0], a[1], zb), (b[0], b[1], zb)) for a, b in zip(D['TOP_BOX'], D['TOP_BOX'][1:] + D['TOP_BOX'][:1])]
    P = {'path': path_plan(), 'profile': plan_of(D['PROFILE'], 'xy', bh.Z_LEG, 'cradle section '),
         'rib': plan_of(D['RIB'], 'xz', 0.0, 'rib '), 'tyre': plan_of(D['TYRE'], 'xz', 0.0, 'tyre space '),
         'top box': plan_of(box, 'xy', zb, 'top cut ')}
    for k, (z, ents) in enumerate(D['SECS']):
        P['section %d' % (k + 1)] = plan_of(ents, 'xy', z, 'trunk section %d ' % (k + 1))
    D['PLANS'] = P
    D['END_Z'] = (bh.Z_LEG, bh.Z_LIP)
    for k in ('BOSSES', 'R_BOSS', 'H_BOSS', 'R_HOLE', 'RIB_W', 'RIB_F', 'END_R'):
        D[k] = getattr(bh, k)
    return D


def chain(E):
    """Edges ordered into a closed chain, each starting where the previous one ends."""
    E = list(E); out = [E.pop(0)]
    while E:
        end = out[-1][-1]
        for k, e in enumerate(E):
            if math.dist(e[1], end) < 1e-5:
                out.append(E.pop(k)); break
            if math.dist(e[-1], end) < 1e-5:
                e = E.pop(k); out.append((e[0],) + tuple(reversed(e[1:]))); break
        else:
            raise RuntimeError('open outline')
    return out


def plan_of(ents, frame, z0=0.0, tag=''):
    """(frame, z0, ops) for a closed outline of ('L', p, q) / ('A', s, mid, t) entities in model mm."""
    sys.path.insert(0, os.path.join(HERE, '..', '..', '..'))
    from stl2cad import sketchplan
    f = (lambda p: (round(p[0], 6), round(p[1], 6))) if frame == 'xy' else (lambda p: (round(p[0], 6), round(p[2], 6)))
    E = [('line', f(e[1]), f(e[2])) if e[0] == 'L' else ('arc', f(e[1]), f(e[2]), f(e[3])) for e in ents]
    ops, full = sketchplan.plan_outline(chain(E), tag=tag)
    def plain(x):                                 # numpy scalars written as plain numbers
        if isinstance(x, (list, tuple)):
            return type(x)(plain(v) for v in x)
        return round(float(x), 6) if isinstance(x, (float, np.floating, np.integer)) and not isinstance(x, bool) else x
    ops = [plain(op) for op in ops]
    if not full:
        print('WARNING: %splan does not fully define the outline' % tag)
    return (frame, z0, ops)


def path_plan():
    """The J path (open): its lines vertical / horizontal, tangent to the arcs, the arcs' radii, the leg's top
    end from the origin, the bottom line's height and length, the lip's top end height."""
    r = bh.R_J; a = math.radians(45)
    L1 = (bh.C_R[0] + r, (bh.Z_LEG + bh.C_R[1]) / 2)
    A1 = (bh.C_R[0] + r * math.cos(a), bh.C_R[1] - r * math.sin(a))
    L2 = (0.0, bh.Z_FLAT)
    A2 = (bh.C_L[0] - r * math.cos(a), bh.C_L[1] - r * math.sin(a))
    L3 = (bh.C_L[0] - r, (bh.C_L[1] + bh.Z_LIP) / 2)
    S, P_ = 'SKETCHSEGMENT', 'SKETCHPOINT'
    top, lip = (bh.C_R[0] + r, bh.Z_LEG), (bh.C_L[0] - r, bh.Z_LIP)
    e2 = ((bh.C_R[0], bh.Z_FLAT), (bh.C_L[0], bh.Z_FLAT))
    ops = [('rel', 'sgVERTICAL2D', [(S, L1)], 'path leg vertical'),
           ('rel', 'sgHORIZONTAL2D', [(S, L2)], 'path bottom horizontal'),
           ('rel', 'sgVERTICAL2D', [(S, L3)], 'path lip vertical')]
    for x, y, w in ((L1, A1, 'leg'), (A1, L2, 'leg arc'), (L2, A2, 'bottom'), (A2, L3, 'lip arc')):
        ops.append(('rel', 'sgTANGENT', [(S, x), (S, y)], 'path %s tangent' % w))
    ops += [('dim', 'any', [(S, A1)], (A1[0] + 3, A1[1] - 3), 'path leg arc radius'),
            ('dim', 'any', [(S, A2)], (A2[0] - 3, A2[1] - 3), 'path lip arc radius'),
            ('dim', 'H', [(P_, top), ('ORIGIN', (0.0, 0.0))], (top[0] / 2, top[1] + 4), 'path leg x'),
            ('dim', 'V', [(P_, top), ('ORIGIN', (0.0, 0.0))], (top[0] + 4, top[1] / 2), 'path leg top'),
            ('dim', 'V', [(S, L2), ('ORIGIN', (0.0, 0.0))], (4.0, bh.Z_FLAT / 2), 'path bottom height'),
            ('dim', 'H', [(P_, e2[0]), (P_, e2[1])], (0.0, bh.Z_FLAT - 4), 'path bottom length'),
            ('dim', 'V', [(P_, lip), ('ORIGIN', (0.0, 0.0))], (lip[0] - 4, lip[1] / 2), 'path lip top')]
    return ('xz', 0.0, ops)


def steps():
    """(name, volume, centre, own volume) after each feature, in SolidWorks' order."""
    S = []
    bos = None
    for x, y in bh.BOSSES:
        c = cq.Solid.makeCylinder(bh.R_BOSS, bh.H_BOSS, cq.Vector(x, y, 0))
        bos = c if bos is None else bos.fuse(c)
    vb, cb = bos.Volume(), np.array(bos.Center().toTuple())
    S.append(('Bosses', vb, cb, vb))
    prev = vb
    for name, s in (('Cradle', bh.cradle().val()), ('Rib', bh.cradle().union(bh.rib()).val()),
                    ('Rib fillet', bh.cradle_with_rib().val()), ('End fillet', bh.end_fillets(bh.cradle_with_rib()).val())):
        v = vb + s.Volume(); c = (vb * cb + s.Volume() * np.array(s.Center().toTuple())) / v
        S.append((name, v, c, max(abs(v - prev), 1.0))); prev = v
    # later features: each adds its exact volume (OCC) minus its overlap with the part so far, from a fine mesh
    # intersection (manifold3d), so the mesh error is that of the overlap only (as frame3d.steps_mesh)
    sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'to'))
    import frame3d
    tol = (0.005, 0.05)
    cur = bh.end_fillets(bh.cradle_with_rib()).val()
    M = frame3d.to_manifold(cur, tol) + frame3d.to_manifold(bos, tol)
    v, mom = S[-1][1], S[-1][1] * np.array(S[-1][2])
    w = [bh.section_wire(np.array(V), z) for z, V in bh.SECTIONS]
    def add(name, solid, cut=False):
        nonlocal M, v, mom
        P = frame3d.to_manifold(solid, tol)
        I = P ^ M
        vi, ci = frame3d.manifold_props(I) if (not I.is_empty() and I.volume() > 0) else (0.0, (0.0, 0.0, 0.0))
        if cut:
            v, mom = v - vi, mom - vi * np.array(ci)
            M = M - P
            S.append((name, v, mom / v, vi))
        else:
            ve = solid.Volume()
            v, mom = v + ve - vi, mom + ve * np.array(solid.Center().toTuple()) - vi * np.array(ci)
            M = M + P
            S.append((name, v, mom / v, ve))
    for k in range(len(w) - 1):
        add('Trunk %d' % (k + 1), cq.Solid.makeLoft(w[k:k + 2], ruled=True))
    for n, a in enumerate(bh.ARMS):
        add('Arm %d' % (n + 1), bh.arm(a))
    add('Top cut', cq.Solid.makeBox(300, 300, 50, cq.Vector(-150, -150, bh.H_BOSS)), cut=True)
    add('Tyre cut', bh.keepout(), cut=True)
    holes = None
    for x, y in bh.BOSSES:
        c = cq.Solid.makeCylinder(bh.R_HOLE, 60, cq.Vector(x, y, -30))
        holes = c if holes is None else holes.fuse(c)
    add('Holes', holes, cut=True)
    part = M
    out = []
    for i, (n, v, c, own) in enumerate(S):
        own = own if own > 0 else abs(v - S[i - 1][1])
        out.append((n, round(v, 3), r3(c), round(max(own, 1.0), 3)))
    return out, part


if __name__ == '__main__':
    D = data()
    S, part = steps()
    print('final %.3f' % S[-1][1])
    lines = []
    for k, v in D.items():
        if k == 'ARMS':
            lines.append('ARMS = [  # spline points; sections (point index, centre, end of major axis, end of minor axis)\n'
                         + ''.join('    dict(points=%r,\n         sec=%r),\n' % (a['points'], a['sec']) for a in v) + ']')
        elif k == 'PLANS':
            lines.append('PLANS = {  # sketch: (plane frame, plane level, relations and dimensions from stl2cad.sketchplan)\n'
                         + ''.join('    %r: %r,\n' % (n, p) for n, p in v.items()) + '}')
        elif k in ('PATH', 'PROFILE', 'RIB', 'TYRE', 'RIB_PTS'):
            lines.append('%s = [\n' % k + ''.join('    %r,\n' % (e,) for e in v) + ']')
        elif k == 'SECS':
            lines.append('SECS = [\n' + ''.join('    (%r, [\n' % z + ''.join('        %r,\n' % (e,) for e in ents) + '    ]),\n'
                                               for z, ents in v) + ']')
        else:
            lines.append('%s = %r' % (k, v))
    lines.append('STEPS = {  # feature: (volume, centre of mass, own volume) after it, mm^3 and mm\n' +
                 ''.join('    %r: (%r, %r, %r),\n' % (n, v, c, o) for n, v, c, o in S) + '}')
    lines.append('FINAL = %r' % ((S[-1][1], S[-1][2]),))
    code = open(os.path.join(HERE, 'sw_template.py')).read()
    code = code.replace('@NAME@', 'bikehook_solidworks').replace('@Z_LEG@', repr(bh.Z_LEG)).replace('@DATA@', '\n'.join(lines))
    open(os.path.join(OUT, 'bikehook_solidworks.py'), 'w').write(code)
    for n, v, c, o in S: print('%-12s %10.3f  own %9.3f' % (n, v, o))
