"""SolidWorks script (pywin32) for a truss-like TO part built from its skeleton: hubs + straight bars.

The model (all along the hub axis, which is the sketch normal):
  Hub n      Extruded Boss of a disk of the TO ring's outer radius, the hub's height;
  Bar n      Extruded Boss of a straight slot (two lines, two end arcs) between two joints, of the
             bar's width, symmetric about one mid-plane, of the bar's depth; slots meeting at a joint
             overlap in their round ends, which forms the joint;
  Hole n     Cut-Extrude of the original's hole;
then the re-entrant edges where bars meet are filleted (swplates' junction fillets).
Each bar is its own feature, so it can be suppressed. Bar widths and depths come in a few classes (the
detail level sets how many); equal end radii in different sketches are linked to one global variable by
swgen's link_dimensions, so a width class is one parameter.

python swbars.py ORIGINAL_EXACT.py TO.stl NAME [--units in] [--detail 1-5] [--plot]
"""
import argparse, math, os, sys
import numpy as np
import cadquery as cq
import trimesh

from stl2cad import swgen
from stl2cad.core import mass_props
from stl2cad.sketchplan import plan_outline
import hubweb, skeleton, swplates

# detail level: number of width and depth classes (0: every bar its own), skeleton tolerances (mm)
DETAIL = {1: dict(widths=1, depths=1, straight=2.0, short=6.0, digits=2),
          2: dict(widths=2, depths=1, straight=1.5, short=5.0, digits=2),
          3: dict(widths=2, depths=2, straight=1.0, short=4.0, digits=2),
          4: dict(widths=3, depths=3, straight=1.0, short=4.0, digits=3),
          5: dict(widths=0, depths=0, straight=0.5, short=3.0, digits=3)}
EPS = 1e-3      # mm: cut prisms reach this far past the faces (OCC leaves coplanar faces whole)


def classes(x, w, k):
    """Values x (weights w) in k classes: optimal 1D partition of the sorted values (least weighted
    squared error, dynamic programming); each value replaced by its class's weighted mean. k = 0 or
    k >= the number of distinct values: x unchanged."""
    x, w = np.asarray(x, float), np.asarray(w, float)
    n = len(x)
    if k <= 0 or k >= n:
        return x.copy()
    o = np.argsort(x); xs, ws = x[o], w[o]
    cw, cx, cxx = (np.concatenate([[0], np.cumsum(a)]) for a in (ws, ws * xs, ws * xs * xs))
    sse = lambda i, j: cxx[j] - cxx[i] - (cx[j] - cx[i]) ** 2 / (cw[j] - cw[i])    # xs[i:j]
    INF = float('inf')
    C = np.full((k + 1, n + 1), INF); B = np.zeros((k + 1, n + 1), int); C[0, 0] = 0
    for m in range(1, k + 1):
        for j in range(m, n + 1):
            for i in range(m - 1, j):
                c = C[m - 1, i] + sse(i, j)
                if c < C[m, j]:
                    C[m, j], B[m, j] = c, i
    out = np.empty(n); j = n
    for m in range(k, 0, -1):
        i = B[m, j]; out[o[i:j]] = (cx[j] - cx[i]) / (cw[j] - cw[i]); j = i
    return out


def _sig(x, n):
    """x to n significant digits, never coarser than 0.1 (as the scripts' rounding)."""
    return x if not n or x == 0 else round(x, max(1, n - 1 - int(math.floor(math.log10(abs(x))))))


def fit(original, to, detail=3):
    """Hubs (hubweb) and skeleton bars with sizes in classes. Returns dict(hubs, joints, bars, mid, axis)
    with bars as dict(i, j, w, d)."""
    D = DETAIL[detail]
    skeleton.STRAIGHT, skeleton.SHORT = D['straight'], D['short']
    hw = hubweb.detect(original, to)
    if hw is None:
        raise ValueError('no hubs: the original has no hole that the TO surrounds')
    S = skeleton.fit(hw['shadow'], hw['hubs'], to, hw['axis'])
    B = S['bars']
    L = np.array([b['length'] for b in B])
    sig = lambda v: _sig(v, D['digits'])
    mid = sig(float(np.average([b['centre'] for b in B], weights=L)))
    lo_h = min(h['lo'] for h in hw['hubs']); hi_h = max(h['hi'] for h in hw['hubs'])
    W = classes([b['width_area'] for b in B], L, D['widths'])
    Dp = classes([b['depth'] for b in B], L * np.array([b['width'] for b in B]), D['depths'])
    bars = []
    for b, w, d in zip(B, W, Dp):
        r = sig(w / 2)                                       # end radius: the sketch dimension
        d = min(sig(d), 2 * min(mid - lo_h, hi_h - mid))     # inside the hubs' height range
        bars.append(dict(i=b['i'], j=b['j'], w=2 * r, d=d))
    return dict(hubs=hw['hubs'], joints=[tuple(map(float, p)) for p in S['joints']], bars=bars, mid=mid,
                axis=hw['axis'], shadow=hw['shadow'], hub_of=S['hub_of'])


def slot_loop(p, q, w):
    """Straight slot between centres p and q, width w: line, end arc at q, line, end arc at p (head to tail)."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    u = (q - p) / np.linalg.norm(q - p); n = np.array([-u[1], u[0]]); r = w / 2
    t = lambda a: tuple(float(c) for c in a)
    return [('line', t(p + r * n), t(q + r * n)),
            ('arc', t(q + r * n), t(q + r * u), t(q - r * n)),
            ('line', t(q - r * n), t(p - r * n)),
            ('arc', t(p - r * n), t(p - r * u), t(p + r * n))]


def slot_plan(loop, i, j):
    """Relations and dimensions of a slot sketch (slot_loop's edge order): the four line-arc tangencies
    and equal end radii as planned by stl2cad, then x and y of each end centre from the origin (labelled
    'joint k x' / 'joint k y': swgen links the same joint's coordinates in all sketches to one global
    variable) and the radius of one end."""
    ops, _ = plan_outline([tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in loop])
    rels = [op for op in ops if op[0] == 'rel']
    a1, a0 = loop[1], loop[3]                                   # end arcs at q (joint j) and p (joint i)
    cen = lambda e: tuple(float(c) for c in (np.array(e[1]) + np.array(e[3])) / 2)
    dims = []
    for k, e in ((i, a0), (j, a1)):
        c = cen(e)
        dims += [('dim', 'H', [('SKETCHPOINT', c), ('ORIGIN', None)], (c[0] / 2, c[1] - 8.0), 'joint %d x' % k),
                 ('dim', 'V', [('SKETCHPOINT', c), ('ORIGIN', None)], (c[0] - 8.0, c[1] / 2), 'joint %d y' % k)]
    dims.append(('dim', 'any', [('SKETCHSEGMENT', a1[2])], tuple(np.array(a1[2]) + (np.array(a1[2]) - cen(a1))),
                 'radius of the bar end at joint %d' % j))
    return rels + dims


def _prism(loop, lo, hi):
    """Solid of a slot or circle loop between levels lo and hi (z)."""
    wp = cq.Workplane('XY').workplane(offset=lo)
    if loop[0][0] == 'circle':
        return wp.center(*loop[0][1]).circle(loop[0][2]).extrude(hi - lo).val()
    e0 = loop[0]
    wp = wp.moveTo(*e0[1])
    for e in loop:
        wp = wp.lineTo(*e[2]) if e[0] == 'line' else wp.threePointArc(e[2], e[3])
    return wp.close().extrude(hi - lo).val()


def features(M):
    """Feature list (name, op, loops, lo, hi, volume, com) and the final solid, built in order."""
    feats, union = [], None

    def add(name, op, loop, lo, hi):
        nonlocal union
        if op == 'add':
            s = _prism(loop, lo, hi)
            union = s if union is None else union.fuse(s).clean()
        else:
            union = union.cut(_prism(loop, lo - EPS, hi + EPS)).clean()
        v, c = mass_props(union)
        feats.append(dict(name=name, op=op, loops=[loop], lo=lo, hi=hi, volume=v, com=np.array(c)))

    for n, h in enumerate(M['hubs'], 1):
        add('Hub %d' % n, 'add', [('circle', tuple(map(float, h['c'])), float(h['ro']))], h['lo'], h['hi'])
    J = M['joints']
    for n, b in enumerate(M['bars'], 1):
        add('Bar %d w%.3g d%.3g' % (n, b['w'], b['d']), 'add', slot_loop(J[b['i']], J[b['j']], b['w']),
            M['mid'] - b['d'] / 2, M['mid'] + b['d'] / 2)
        feats[-1]['joints'] = (b['i'], b['j'])
    for n, h in enumerate(M['hubs'], 1):
        add('Hole %d' % n, 'cut', [('circle', tuple(map(float, h['c'])), float(h['ri']))], h['lo'], h['hi'])
    return feats, union


def write(path, name, M, r_fillet=1.0, nd=6, digits=3, detail=None):
    feats, union = features(M)
    if M['axis'] != 2:
        raise NotImplementedError('hub axis along %s: only z for now' % 'xyz'[M['axis']])
    k = M['axis']
    R0, back = swgen.BUILD_FRAME[k]
    i0, j0 = [m for m in range(3) if m != k]
    bb = union.BoundingBox()
    size = bb.DiagonalLength
    t3 = lambda q: '(%.5f, %.5f, %.5f)' % tuple(float(x) for x in q)

    def plans_of(loops, joints=None):
        if joints is not None:          # a slot: stl2cad's relations, then the joints' coordinates and one radius
            return [slot_plan(loops[0], *joints)], True
        plans, full = [], True
        for lp in loops:
            if len(lp) == 1:
                plans.append([]); continue
            try:
                ops, ok = plan_outline([tuple([e[0]] + [tuple(map(float, q)) for q in e[1:]]) for e in lp])
            except Exception:
                ops, ok = [], False
            plans.append(ops); full = full and ok
        return plans, full

    nb = len(M['bars'])
    ws = sorted({round(b['w'], 6) for b in M['bars']}); ds = sorted({round(b['d'], 6) for b in M['bars']})
    L = ['"""%s: SolidWorks feature tree for a topology-optimized part, generated by stl2cad (TO skeleton).' % name, '',
         'The part is built as %d hub(s) and %d straight bar(s) along the skeleton of the TO result. A bar is' % (len(M['hubs']), nb),
         'an Extruded Boss of a slot sketch (two lines and two end arcs) between two joints, symmetric about',
         'the mid-plane at %.4g mm; slots meeting at a joint overlap in their round ends. Bar widths: %s mm;' % (M['mid'], ', '.join('%.3g' % w for w in ws)),
         'depths: %s mm. The original\'s holes are cut last. After every feature the volume and centre of' % ', '.join('%.3g' % d for d in ds),
         'mass are checked against stl2cad\'s values. Re-entrant edges are filleted last (radius R_FILLET).', '',
         'Needs: Windows, SolidWorks, Python with pywin32 (pip install pywin32).',
         'Run:   python %s.py [C:\\\\path\\\\part.SLDPRT]' % name, '"""', '',
         'R_FILLET = %.2f                # mm, fillet radius on re-entrant edges (0: none)' % r_fillet,
         'ROUND_DIGITS = %d                # significant digits of the dimensions, set after the checks (0: as fitted)' % digits,
         'DEFER_ROUNDING = True           # rounding after the junction fillets',
         'DIM_DECIMALS = 2                # decimal places shown on length dimensions (document units: mm)',
         '# detail level %s (1: fewest dimensions ... 5: smallest model error)' % detail, '',
         'AXIS_I, AXIS_J = %d, %d           # model axes of the sketch coordinates (a, b)' % (i0, j0),
         'PLANE_NAME = "Front Plane"', 'PLANE_INDEX = 0',
         'ROT = %s' % [list(map(float, r)) for r in R0],
         'ROT_BACK = %s' % (tuple(back),),
         'FINAL_CENTROID = %s   # after the first feature, STL frame, mm' % t3(feats[0]['com']),
         'SIZE = %s' % swgen._f(size, 3)]
    allp = np.array([q for f in feats for lp in f['loops'] for e in lp for q in e[1:] if np.ndim(q) == 1])
    L.append('HELPER_AB = (%s, %s)' % (swgen._f(allp[:, 0].min() - 10, 3), swgen._f(allp[:, 1].min() - 10, 3)))
    L += ['', '# ---- features: sketch loops (a, b in mm), levels along the hub axis, expected volume (mm^3) ----',
          'FEATURES = [']
    for ft in feats:
        plans, full = plans_of(ft['loops'], ft.get('joints'))
        cb = np.asarray(R0, float) @ ft['com']
        L += ['    dict(',
              '        name="%s", op="%s", lo=%s, hi=%s, volume=%.4f, centroid=%s, holes=False,'
              % (ft['name'], ft['op'], swgen._f(ft['lo'], nd), swgen._f(ft['hi'], nd), ft['volume'], t3(cb)),
              '        expect="%s",' % ('fully defined' if full else 'under defined'),
              '        plans=' + repr(swgen._round_plans(plans)) + ',', '        loops=[']
        for lp in ft['loops']:
            L.append('        ' + swgen._loop_code(swplates._Loop(lp), max(nd, 9)) + ',')
        L += ['        ],', '    ),']
    L += [']', '', 'CYLINDERS = []', 'REVOLVES = []', 'REV_TURN = None', 'REV_DONE = None', 'CUTS = []',
          'XCUTS = []', 'BLENDS = []']
    E = swplates.concave_edges(union)
    L += ['', '# ---- re-entrant edges (STL frame, mm): midpoint, length, direction into the empty corner,',
          '# fillet area / r^2 per mm (angle of the empty wedge) ----', 'JUNCTIONS = [']
    nj = 0
    for e in E:
        P = e.positionAt(0.5)
        b = swplates._corner_direction(union, e, P)
        kk, phi = swplates._wedge(union, e, P, b)
        if phi > swplates.MAX_WEDGE:
            continue
        nj += 1
        L.append('    (%s, %.4f, %s, %.5f),   # %.0f deg' % (t3(P.toTuple()), e.Length(), t3(b), kk, phi))
    L += [']', '']
    head, _ = swgen._RUNTIME.split('if __name__ == "__main__":')
    L.append(head + swplates.FILLET_CODE)
    open(path, 'w').write('\n'.join(L) + '\n')
    return dict(volume=feats[-1]['volume'], concave=nj, union=union, features=len(feats))


if __name__ == '__main__':
    import time
    ap = argparse.ArgumentParser()
    ap.add_argument('base', help="stl2cad's exact script of the original (out_NAME/NAME_exact.py)")
    ap.add_argument('to'); ap.add_argument('name')
    ap.add_argument('--units', default='in'); ap.add_argument('--detail', type=int, default=3, choices=range(1, 6))
    ap.add_argument('--fillet', type=float, default=1.0)
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--no-verify', action='store_true')
    a = ap.parse_args()
    t0 = time.time()
    ns = swplates._namespace(a.base)
    to = trimesh.load(a.to); to.apply_scale({'m': 1000.0, 'mm': 1.0, 'in': 25.4}[a.units])
    to = max(to.split(only_watertight=False), key=lambda b: len(b.faces))
    M = fit(ns['original'], to, a.detail)
    D = DETAIL[a.detail]
    print('skeleton: %d joints (%d on hubs), %d bars; mid-plane %.3g mm; widths %s; depths %s'
          % (len(M['joints']), len(M['hub_of']), len(M['bars']), M['mid'],
             sorted({'%.3g' % b['w'] for b in M['bars']}), sorted({'%.3g' % b['d'] for b in M['bars']})))
    info = write(a.name + '_solidworks.py', a.name + '_solidworks', M, a.fillet, digits=D['digits'], detail=a.detail)
    s = info['union']
    cq.exporters.export(cq.Workplane().add(s), a.name + '.step')
    cq.exporters.export(cq.Workplane().add(s), a.name + '.stl', tolerance=0.02, angularTolerance=0.2)
    m = trimesh.load(a.name + '.stl')
    print('SolidWorks script: %s_solidworks.py (%d features, %d re-entrant edge(s) to fillet)'
          % (a.name, info['features'], info['concave']))
    print('model: valid %s, solids %d; volume: TO %.0f, model %.0f mm^3 (%+.1f %%)'
          % (s.isValid(), len(s.Solids()), to.volume, m.volume, 100 * (m.volume / to.volume - 1)))
    lo = np.minimum(to.bounds[0], m.bounds[0]) - 1; pitch = 1.0
    shape = tuple(np.ceil((np.maximum(to.bounds[1], m.bounds[1]) + 1 - lo) / pitch).astype(int) + 1)

    def vox(mesh):
        g = mesh.voxelized(pitch).fill(); idx = np.argwhere(g.matrix)
        pts = trimesh.transformations.transform_points(idx, g.transform)
        G = np.zeros(shape, bool); k = np.round((pts - lo) / pitch).astype(int); G[k[:, 0], k[:, 1], k[:, 2]] = True
        return G
    A, B = vox(to), vox(m); both = (A & B).sum()
    print('voxel overlap (1 mm voxels): %.0f %% of the TO is in the model, %.0f %% of the model is in the TO'
          % (100 * both / A.sum(), 100 * both / B.sum()))
    if a.plot:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        from shapely.geometry import LineString
        fig, ax = plt.subplots(figsize=(14, 4))
        g = M['shadow']; ax.fill(*g.exterior.xy, c='0.85'); [ax.fill(*h.xy, c='w') for h in g.interiors]
        J = M['joints']
        for b in M['bars']:
            pg = LineString([J[b['i']], J[b['j']]]).buffer(b['w'] / 2)
            ax.plot(*pg.exterior.xy, 'r-', lw=0.8)
        for h in M['hubs']:
            t = np.linspace(0, 2 * np.pi, 60)
            for r in (h['ri'], h['ro']):
                ax.plot(h['c'][0] + r * np.cos(t), h['c'][1] + r * np.sin(t), 'b', lw=0.8)
        ax.set_aspect('equal'); ax.set_title('%s: %d bars (red: slot outlines), TO shadow grey' % (a.name, len(M['bars'])))
        plt.tight_layout(); plt.savefig(a.name + '_skeleton.png', dpi=90)
    if not a.no_verify:
        ok, log = swplates.verify(a.name + '_solidworks.py', s)
        open(a.name + '_standin.log', 'w').write(log)
        print('API stand-in: %s; %s' % ('ok' if ok else 'FAILED', log.strip().splitlines()[-1]))
    print('time %.0f s' % (time.time() - t0))
