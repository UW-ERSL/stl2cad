"""SolidWorks script (pywin32) for a TO result built as plates of the original part.

Each plate becomes one Extruded Boss from a sketch of its exact region (original section minus
the plate's cuts, as lines, arcs and circles). The first plate is built along its own axis; the
others with the body turned so their axis is z (stl2cad's Move/Copy Body sequence), each checked
against the volume and centre of mass of the CadQuery plate model. Re-entrant edges of the result
(face normals: BRepOffset_Analyse) are filleted last, picked at their midpoints; the radius is a
parameter of the script.
"""
import math, os
import numpy as np
import cadquery as cq
from cadquery import Vector
from OCP.BRepTools import BRepTools_WireExplorer
from OCP.BRepOffset import BRepOffset_Analyse
from OCP.ChFiDS import ChFiDS_TypeOfConcavity

from stl2cad import swgen
from stl2cad.core import mass_props
from stl2cad.sketchplan import plan_outline
import codegen3d


def _namespace(base_script):
    """Original part (exec of stl2cad's exact script) and the plate helpers of codegen3d."""
    ns = {'__name__': 'base'}
    exec(compile(open(base_script).read(), base_script, 'exec'), ns)
    exec(compile('from cadquery import Edge, Solid, Vertex\n'
                 'from OCP.BRepOffset import BRepOffset_Analyse\n'
                 'from OCP.ChFiDS import ChFiDS_TypeOfConcavity\n' + codegen3d.HELPERS, 'helpers', 'exec'), ns)
    return ns


def _ordered_edges(wire):
    ex = BRepTools_WireExplorer(wire.wrapped)
    out = []
    while ex.More():
        out.append(cq.Edge(ex.Current()))
        ex.Next()
    return out


def region_loops(solid, axis, level, tol=1e-6):
    """Loops of the plate's faces lying in the plane at level along axis, in (u, v) coordinates:
    ('line', p, q), ('arc', start, mid, end), ('circle', centre, r)."""
    u, v = [i for i in range(3) if i != axis]
    uv = lambda p: (float(p[u]), float(p[v]))
    loops = []
    for f in solid.Faces():
        if f.geomType() != 'PLANE':
            continue
        n = f.normalAt().toTuple(); c = f.Center().toTuple()
        if abs(abs(n[axis]) - 1) > 1e-9 or abs(c[axis] - level) > tol:
            continue
        for w in [f.outerWire()] + f.innerWires():
            es = []
            for e in _ordered_edges(w):
                t = e.geomType()
                p0, p1 = e.startPoint().toTuple(), e.endPoint().toTuple()
                if t == 'LINE':
                    es.append(('line', uv(p0), uv(p1)))
                elif t == 'CIRCLE':
                    if e.IsClosed():
                        es.append(('circle', uv(e.arcCenter().toTuple()), float(e.radius())))
                    else:
                        es.append(('arc', uv(p0), uv(e.positionAt(0.5).toTuple()), uv(p1)))
                else:
                    raise ValueError('plate outline edge of type %s' % t)
            loops.append(_chain(es))
    return loops


def _rev(e):
    return (e[0], e[-1]) + tuple(reversed(e[2:-1])) + (e[1],) if e[0] == 'arc' else (e[0], e[2], e[1])


def _chain(es):
    """Edges of a closed loop head to tail (the explorer's order, each edge turned as needed)."""
    if len(es) == 1 or es[0][0] == 'circle':
        return es
    nxt = es[1]
    if min(math.dist(es[0][1], nxt[1]), math.dist(es[0][1], nxt[-1])) < min(math.dist(es[0][-1], nxt[1]), math.dist(es[0][-1], nxt[-1])):
        es[0] = _rev(es[0])
    for i in range(1, len(es)):
        if math.dist(es[i][1], es[i - 1][-1]) > math.dist(es[i][-1], es[i - 1][-1]):
            es[i] = _rev(es[i])
    return es


class _Loop:                                   # the loop object swgen._loop_code expects
    def __init__(self, edges):
        self.edges, self.pattern, self.linear = edges, None, None


EPS = 1e-3                 # mm
MAX_WEDGE = 150.0          # deg: re-entrant edges with a flatter empty wedge are not filleted


def _axis_pick(solid, axis, centre, d=0.5):
    """A ray (start, direction; mm) that hits a cylindrical face of solid whose axis is the plate axis
    through centre (in-plane coordinates), fired from the empty side; None if there is none."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_Cylinder
    u, v = [i for i in range(3) if i != axis]
    for f in sorted(solid.Faces(), key=lambda f: -f.Area()):
        ad = BRepAdaptor_Surface(f.wrapped)
        if ad.GetType() != GeomAbs_Cylinder:
            continue
        ax = ad.Cylinder().Axis(); L = ax.Location(); D = ax.Direction()
        Lp, Dp = np.array([L.X(), L.Y(), L.Z()]), np.array([D.X(), D.Y(), D.Z()])
        if abs(abs(Dp[axis]) - 1) > 1e-6 or math.hypot(Lp[u] - centre[0], Lp[v] - centre[1]) > 1e-3:
            continue
        R = ad.Cylinder().Radius()
        bb = f.BoundingBox()
        mid = 0.5 * ((bb.xmin, bb.ymin, bb.zmin)[axis] + (bb.xmax, bb.ymax, bb.zmax)[axis])
        for t in np.linspace(0, 2 * math.pi, 24, endpoint=False):    # a point of the face half way up
            r = np.zeros(3); r[u], r[v] = math.cos(t), math.sin(t)
            P = Lp.copy(); P[axis] = mid; P = P + R * r
            if f.distance(cq.Vertex.makeVertex(*P)) > 1e-6:
                continue
            for s_ in (1, -1):
                if not solid.isInside(Vector(*(P + s_ * d * r)), 1e-5):
                    return (P + s_ * d * r, -s_ * r)
    return None


def _wedge(solid, e, P, b):
    """Fillet area per unit length and radius^2 at a re-entrant edge: cot(phi/2) - (pi - phi)/2, phi the
    angle of the empty wedge between the two faces (pi - the angle between their outward normals,
    oriented here to point into the corner, direction b); 1 - pi/4 for perpendicular faces, 0 when flat."""
    n = [v if v @ b >= 0 else -v for v in _face_normals(solid, e, P)]
    phi = math.pi - math.acos(float(np.clip(n[0] @ n[1], -1, 1)))
    return 1.0 / math.tan(phi / 2) - (math.pi - phi) / 2, math.degrees(phi)


def _corner_direction(solid, e, P, d=0.05):
    """Unit vector at P (a point of re-entrant edge e) into the empty corner: the one of +-n1 +- n2
    (n1, n2 the two face normals, whatever their orientation) that leads outside the solid."""
    n = _face_normals(solid, e, P)
    for s1 in (1, -1):
        for s2 in (1, -1):
            v = s1 * n[0] + s2 * n[1]
            v = v / (np.linalg.norm(v) or 1.0)
            if not solid.isInside(P + Vector(*(d * v)), 1e-5):
                return v
    return v


def _face_normals(solid, e, P):
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
    m = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(solid.wrapped, TopAbs_EDGE, TopAbs_FACE, m)
    return [np.array(cq.Face(f).normalAt(P).toTuple()) for f in m.FindFromKey(e.wrapped)]


def concave_edges(solid, d=0.05):
    """Re-entrant edges: BRepOffset_Analyse's concave edges, confirmed by the material around the
    edge midpoint: of the four points P + d(+-n1 +- n2) (n1, n2 the two face normals) a re-entrant
    edge has three inside the solid, a convex edge one, a seam between coplanar faces two. The
    analyser alone also returned convex edges of the plate union (face orientation)."""
    an = BRepOffset_Analyse(solid.wrapped, 0.01)
    out = []
    for e in solid.Edges():
        if not any(i.Type() == ChFiDS_TypeOfConcavity.ChFiDS_Concave for i in an.Type(e.wrapped)):
            continue
        P = e.positionAt(0.5)
        n = _face_normals(solid, e, P)
        if len(n) != 2:
            continue
        inside = sum(solid.isInside(P + Vector(*(d * (s1 * n[0] + s2 * n[1]))), 1e-5)
                     for s1 in (1, -1) for s2 in (1, -1))
        if inside == 3:
            out.append(e)
    return out


def write(path, base_script, name, plates, box, r_fillet=1.0, nd=6, digits=3, detail=None, per_cut=True):
    ns = _namespace(base_script)
    t3 = lambda q: '(%.5f, %.5f, %.5f)' % tuple(float(x) for x in q)
    solids = [ns['plate'](p['axis'], p['lo'], p['hi'], p['cuts'], box, p.get('ring')) for p in plates]
    # largest plate first; plates with an 'order' key (nested height-field plates) in that order
    order = sorted(range(len(plates)), key=lambda i: (plates[i].get('order', 0), -mass_props(solids[i])[0]))
    k = plates[order[0]]['axis']
    R0, back = swgen.BUILD_FRAME[k]
    i0, j0 = [m for m in range(3) if m != k]
    size = float(np.linalg.norm(np.array(box[1]) - np.array(box[0])))

    def prism_of(cut, a, lo, hi):
        """The cut's prism between the plate's two levels (the Cut-Extrude's extent)."""
        if cut[0] == 'poly':
            w = ns['rounded_polygon'](cut[1], a, lo - 1)
        else:
            (cx, cy, r, a0, a1), closing = cut[1], cut[2]
            w = ns['arc_notch'](cx, cy, r, a0, a1, closing, a, lo - 1)
        # 1e-3 mm past the plate's faces: a prism ending in them leaves those faces whole (OCC)
        return ns['through'](w, a, lo, hi).intersect(ns['slab'](a, lo - EPS, hi + EPS, box))

    # Each plate: the original between its two levels (Extruded Boss), then, if per_cut, one Cut-Extrude per
    # cut through the plate (each can be suppressed). The part is built plate by plate in one body, so a
    # plate's cut also trims other plates where they lie inside it, and a plate's boss refills their cuts
    # there (about 0.3 % of the volume for the I-beam cantilevers); the expected values follow that order.
    union, steps = None, []
    for n, idx in enumerate(order, 1):
        p = plates[idx]; a, lo, hi = p['axis'], p['lo'], p['hi']
        before = None if union is None else np.array(mass_props(union)[1])
        tag = '%s %d ' % (p.get('kind', 'Plate'), n) if len(plates) > 1 else ''
        feats = []
        if per_cut:
            boss = ns['plate'](a, lo, hi, [], box, p.get('ring'))
            union = boss if union is None else union.fuse(boss).clean()
            v, c = mass_props(union)
            kind = p.get('kind', 'Plate')
            feats.append(dict(name='%s %d %s %.4g to %.4g' % (kind, n, swgen.AX[a], lo, hi) if len(plates) > 1
                              else '%s %s %.4g to %.4g' % (kind, swgen.AX[a], lo, hi), op='add',
                              loops=region_loops(boss, a, lo), volume=v, com=np.array(c)))
            nv = nn = 0
            sec = p.get('sec')
            pat = p.get('pattern') if n == 1 else None      # a Circular Pattern only on the first plate
            cuts = p['cuts']
            if pat:        # seeds, then the cuts outside the pattern, then the copies (seed by seed)
                ns_, nc = pat['seeds'], pat['copies']
                seq = list(range(ns_)) + list(range(ns_ + nc, len(cuts))) + list(range(ns_, ns_ + nc))
            else:
                seq = list(range(len(cuts)))
            names, copies = {}, []
            for i in seq:
                cut = cuts[i]
                prism = prism_of(cut, a, lo, hi)
                nxt = union - prism
                v1, c1 = mass_props(nxt)
                is_copy = pat is not None and ns_ <= i < ns_ + nc
                if v - v1 < 1e-3 and not is_copy:
                    continue                          # removes nothing
                if is_copy:
                    k = i - ns_; seed = k // (pat['n'] - 1)
                    nm = '%s copy %d' % (names.get(seed, 'cut %d' % seed), k % (pat['n'] - 1) + 2)
                else:
                    notch = sec is not None and not sec.cut_polygon(cut).within(sec.Oout.buffer(1e-6))
                    nn, nv = (nn + 1, nv) if notch else (nn, nv + 1)
                    nm = tag + ('Notch %d' % nn if notch else 'Void %d' % nv)
                    names[i] = nm
                ft = dict(name=nm, op='cut', loops=region_loops(prism, a, lo - EPS), volume=v1, com=np.array(c1))
                (copies if is_copy else feats).append(ft)
                union, v = nxt, v1
            if pat:
                feats.append(dict(name=tag + 'Pattern of cuts', op='pattern', n=pat['n'],
                                  seeds=[names[i] for i in range(ns_) if i in names], copies=copies,
                                  axis_pick=_axis_pick(union, a, pat['centre']), volume=v, com=np.array(c1),
                                  loops=[]))
        else:
            s_ = solids[idx]
            union = s_ if union is None else union.fuse(s_).clean()
            v, c = mass_props(union)
            feats.append(dict(name='Plate %d %s %.4g to %.4g' % (n, swgen.AX[a], lo, hi), op='add',
                              loops=region_loops(s_, a, lo), volume=v, com=np.array(c)))
        steps.append(dict(p=p, idx=idx, feats=feats, loops=[lp for f in feats for lp in f['loops']],
                          volume=feats[-1]['volume'], com=feats[-1]['com'], before=before))

    def plans_of(loops):
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

    L = ['"""%s: SolidWorks feature tree for a topology-optimized part, generated by stl2cad (TO plates).' % name, '',
         'The part is built as %d plate(s) of the original part: each plate is one Extruded Boss from a' % len(plates),
         'sketch of its region (original section minus the material the optimization removed), between',
         'two levels along its axis. The first plate is built along its own axis on the Front Plane; the',
         'others with the body turned (Move/Copy Body) so that their axis is z, then turned back. After',
         'every feature the volume and centre of mass are checked against stl2cad\'s values. Re-entrant',
         'edges where plates meet are filleted last (radius R_FILLET; 0 = no fillets).', '',
         'Needs: Windows, SolidWorks, Python with pywin32 (pip install pywin32).',
         'Run:   python %s.py [C:\\\\path\\\\part.SLDPRT]' % name, '"""', '',
         'R_FILLET = %.2f                # mm, fillet radius on re-entrant edges (0: none)' % r_fillet,
         'ROUND_DIGITS = %d                # significant digits of the dimensions, set after the checks (0: as fitted)' % digits,
         'DEFER_ROUNDING = True           # rounding after the junction fillets',
         'DIM_DECIMALS = 2                # decimal places shown on length dimensions (document units: mm)',
         '# detail level %s (1: fewest dimensions ... 5: smallest model error)' % detail, '',
         'AXIS_I, AXIS_J = %d, %d           # model axes of the sketch coordinates (a, b) of the first plate' % (i0, j0),
         'PLANE_NAME = "Front Plane"', 'PLANE_INDEX = 0',
         'ROT = %s' % [list(map(float, r)) for r in R0],
         'ROT_BACK = %s' % (tuple(back),),
         'FINAL_CENTROID = %s   # after the first plate, STL frame, mm' % t3(steps[0]['com']),
         'SIZE = %s' % swgen._f(size, 3)]
    allp = np.array([q for lp in steps[0]['loops'] for e in lp for q in e[1:] if np.ndim(q) == 1])
    L.append('HELPER_AB = (%s, %s)' % (swgen._f(allp[:, 0].min() - 10, 3), swgen._f(allp[:, 1].min() - 10, 3)))
    # ---- first plate: FEATURES -----------------------------------------------------------------
    st = steps[0]; p = st['p']; feats = st['feats']
    L += ['', '# ---- first plate: sketch loops (a, b in mm), levels along its axis, expected volume (mm^3) ----',
          'FEATURES = [']
    def feature_code(ft, ind):
        plans, full = plans_of(ft['loops'])
        cb = np.asarray(R0, float) @ ft['com']
        out = [ind + 'dict(',
               ind + '    name="%s", op="%s", lo=%s, hi=%s, volume=%.4f, centroid=%s, holes=False,'
               % (ft['name'], ft['op'], swgen._f(p['lo'], nd), swgen._f(p['hi'], nd), ft['volume'], t3(cb)),
               ind + '    expect="%s",' % ('fully defined' if full else 'under defined'),
               ind + '    plans=' + repr(swgen._round_plans(plans)) + ',', ind + '    loops=[']
        for lp in ft['loops']:
            out.append(ind + '    ' + swgen._loop_code(_Loop(lp), max(nd, 9)) + ',')
        return out + [ind + '    ],', ind + '),']

    for ft in feats:
        if ft['op'] != 'pattern':
            L += feature_code(ft, '    ')
            continue
        ap = ft['axis_pick']
        L += ['    dict(',
              '        name="%s", op="pattern", n=%d, seeds=%r, volume=%.4f,' % (ft['name'], ft['n'], ft['seeds'], ft['volume']),
              '        axis_pick=%s,   # ray at a cylindrical face on the pattern axis (STL frame)'
              % ('None' if ap is None else '(%s, %s)' % (t3(ap[0]), t3(ap[1]))),
              '        copies=[   # built one by one when the pattern fails']
        for c in ft['copies']:
            L += feature_code(c, '            ')
        L += ['        ],', '    ),']
    L += [']', '', 'CYLINDERS = []', 'REVOLVES = []', 'REV_TURN = None', 'REV_DONE = None', 'CUTS = []']
    # ---- other plates: XCUTS entries (body turned so the plate axis is z: its boss, then its cuts) ----
    L += ['', '# ---- other plates: the body is turned so the plate axis is z, the plate added and cut, the body',
          '# turned back (turn_in on the first entry of a plate, turn_out on its last) ----',
          'XCUTS = [']
    for st in steps[1:]:
        p = st['p']; j = p['axis']
        m, ang = swgen.XCUT_TURN[j]
        R = swgen._rot(m, ang)
        u, v = [i for i in range(3) if i != j]

        def to_sk(q):
            P3 = np.zeros(3); P3[u], P3[v] = q
            return tuple(float(x) for x in (R @ P3)[:2])
        sgn = float((R @ np.eye(3)[j])[2])
        lo, hi = sorted((sgn * p['lo'], sgn * p['hi']))
        allq = []
        fl = []
        for ft in st['feats']:
            loops = []
            for lp in ft['loops']:
                es = []
                for e in lp:
                    if e[0] == 'circle':
                        es.append(('circle', to_sk(e[1]), e[2]))
                    else:
                        es.append((e[0],) + tuple(to_sk(q) for q in e[1:]))
                loops.append(es)
            fl.append(loops)
            allq += [q for lp in loops for e in lp for q in e[1:] if np.ndim(q) == 1]
        P = np.array(allq)
        nf = len(st['feats'])
        for i, (ft, loops) in enumerate(zip(st['feats'], fl)):
            plans, full = plans_of(loops)
            L += ['    dict(',
                  '        name="%s", op="%s", axis=%s, angle=%.12f, lo=%s, hi=%s, volume=%.4f, turn_in=%s, turn_out=%s,'
                  % (ft['name'], ft['op'], m, ang, swgen._f(lo, nd), swgen._f(hi, nd), ft['volume'], i == 0, i == nf - 1),
                  '        turned=%s, centroid=%s, back=%s, holes=False,'
                  % (t3(R @ st['before']), t3(R @ ft['com']), t3(ft['com'])),
                  '        helper=(%s, %s),' % (swgen._f(P[:, 0].min() - 10, 3), swgen._f(P[:, 1].min() - 10, 3)),
                  '        expect="%s",' % ('fully defined' if full else 'under defined'),
                  '        plans=' + repr(swgen._round_plans(plans)) + ',', '        loops=[']
            for lp in loops:
                L.append('        ' + swgen._loop_code(_Loop(lp), max(nd, 9)) + ',')
            L += ['        ],', '    ),']
    L += [']', 'BLENDS = []']
    # ---- re-entrant edges: picked at their midpoints, with the volume each fillet adds ----------
    E = concave_edges(union)
    if os.environ.get('SWPLATES_DUMP'):
        cq.exporters.export(cq.Workplane().add(union), os.environ['SWPLATES_DUMP'])
    L += ['', '# ---- re-entrant edges (STL frame, mm): midpoint, length, direction into the empty corner (the',
          '# edge is picked with a ray fired back along it), fillet area / r^2 per mm (angle of the empty wedge) ----',
          'JUNCTIONS = [']
    nj = 0
    for e in E:
        P = e.positionAt(0.5)
        b = _corner_direction(union, e, P)
        k, phi = _wedge(union, e, P, b)
        if phi > MAX_WEDGE:                   # nearly flat: a fillet adds next to nothing
            continue
        nj += 1
        L.append('    (%s, %.4f, %s, %.5f),   # %.0f deg' % (t3(P.toTuple()), e.Length(), t3(b), k, phi))
    L += [']', '']
    rt = swgen._RUNTIME
    head, main = rt.split('if __name__ == "__main__":')
    L.append(head + FILLET_CODE)
    open(path, 'w').write('\n'.join(L) + '\n')
    return dict(volume=steps[-1]['volume'], concave=nj, union=union)


FILLET_CODE = r'''
def junction_fillets(sw, r):
    """Fillet the re-entrant edges with radius r: all in one Fillet feature, else edge by edge.
    A fillet is kept when the volume grows by at least 0.25 times its area (r^2 (cot(phi/2) - (pi - phi)/2)
    per mm of edge, phi the angle of the empty wedge) and at most 1 % of the part; SolidWorks' failures are
    reported, not fatal."""
    if r <= 0 or not JUNCTIONS:
        return
    if not sw.turned_back:
        print("NOTE: part not in the STL's orientation; no fillets")
        return
    fm, ext = sw.model.FeatureManager, sw.model.Extension
    E = VARIANT(pythoncom.VT_EMPTY, None)
    rm = r * MM
    makes = [("FeatureFillet3 (no propagation)", lambda: fm.FeatureFillet3(194, rm, rm, 0, 0, 0, 0, E, E, E, E, E, E, E)),
             ("FeatureFillet3", lambda: fm.FeatureFillet3(195, rm, rm, 0, 0, 0, 0, E, E, E, E, E, E, E)),
             ("FeatureFillet3 (Nothing)", lambda: fm.FeatureFillet3(195, rm, rm, 0, 0, 0, 0, NOTHING, NOTHING,
                                                                    NOTHING, NOTHING, NOTHING, NOTHING, NOTHING)),
             ("FeatureFillet", lambda: fm.FeatureFillet(195, rm, 0, 0, E, E, E))]

    def pick(p, b, mark):
        """Select the edge at p. SelectByID2 with a point picks the first entity along the line of sight
        (seen in SolidWorks: an outer edge in front of the junction edge), so the edge is picked with a ray
        fired at p from 0.3 mm inside the empty corner (direction b); SelectByID2 is the fall-back."""
        q, d = sw.in_model_frame(p), sw.in_model_frame(b)
        start = [(q[i] + 0.3 * d[i]) * MM for i in range(3)]
        n0 = sw.model.SelectionManager.GetSelectedObjectCount2(-1)
        if ext.SelectByRay(*start, -d[0], -d[1], -d[2], 0.05 * MM, 1, True, mark, 0):     # 1: swSelEDGES
            if sw.model.SelectionManager.GetSelectedObjectCount2(-1) > n0:
                return True
        return ext.SelectByID2("", "EDGE", *(c * MM for c in q), True, mark, NOTHING, 0)

    def attempt(edges, name):
        v0 = sw.mass()[0]
        est = sum(L * k for _, L, _, k in edges) * r * r
        for mark in (1, 0):
            for how, mk in makes:
                sw.model.ClearSelection2(True)
                if not all(pick(p, b, mark) for p, _, b, _ in edges):
                    continue
                f = None
                try:
                    f = mk()
                except Exception as ex:
                    print(f"  {name}: {how} failed: {type(ex).__name__}: {ex}")
                if f is None:
                    continue
                dv = sw.mass()[0] - v0
                # kept if it adds material (a fillet on a convex edge removes it), at least about the
                # amount for the picked edges; more is allowed: SolidWorks runs the fillet on along
                # tangent-continuous re-entrant edges (a notch side meeting a flange)
                if 0.25 * est <= dv <= 0.01 * v0:
                    f.Name = name
                    return True
                print(f"  {name}: {how} added {dv:.3f} mm^3, expected at least {0.25 * est:.3f}; removed")
                sw.model.ClearSelection2(True)
                f.Select2(False, 0)
                ext.DeleteSelection2(0)
        return False

    if attempt(JUNCTIONS, "Fillet re-entrant edges"):
        print(f"re-entrant edges: {len(JUNCTIONS)} filleted (R {r} mm) in one feature")
        return
    done, failed = 0, []
    for n, e in enumerate(JUNCTIONS, 1):
        if attempt([e], f"Fillet re-entrant edge {n}"):
            done += 1
        else:
            failed.append(e[0])
    print(f"re-entrant edges: {done} of {len(JUNCTIONS)} filleted (R {r} mm)"
          + (f"; not filleted at {failed}" if failed else ""))


if __name__ == "__main__":
    sw = build()
    junction_fillets(sw, R_FILLET)              # before rounding: the edges are picked at their computed midpoints
    print(sw.round_dimensions())
    if len(sys.argv) > 1:
        sw.model.SaveAs3(sys.argv[1], 0, 0)
'''


def invariants(a, b, size, rtol=1e-4):
    """Volume, surface area, centre of mass and inertia tensor (about the centre of mass) of two
    solids agree."""
    from OCP.GProp import GProp_GProps
    from OCP.BRepGProp import BRepGProp
    def props(x):
        g = GProp_GProps(); BRepGProp.VolumeProperties_s(x.wrapped, g)
        c = g.CentreOfMass(); M = g.MatrixOfInertia()
        I = np.array([[M.Value(i, j) for j in (1, 2, 3)] for i in (1, 2, 3)])
        return g.Mass(), x.Area(), np.array([c.X(), c.Y(), c.Z()]), I
    va, aa, ca, Ia = props(a); vb, ab, cb, Ib = props(b)
    return (abs(va - vb) <= rtol * vb and abs(aa - ab) <= 1e-3 * ab and np.linalg.norm(ca - cb) <= 1e-5 * size + 1e-4
            and np.abs(Ia - Ib).max() <= rtol * np.abs(Ib).max())


def verify(script_path, ref, fillets=False, timeout=900):
    """Run the script on stl2cad's CadQuery stand-in for the SolidWorks API (swmock) and compare the
    result with ref (the plate model without fillets). With fillets=False, R_FILLET is set to 0."""
    import os, re, subprocess, sys, tempfile
    import stl2cad
    d = tempfile.mkdtemp(prefix='swplates_')
    py, out = os.path.join(d, 'part.py'), os.path.join(d, 'part.step')
    src = open(script_path).read()
    if not fillets:
        src = re.sub(r'^R_FILLET = [0-9.]+', 'R_FILLET = 0.0', src, flags=re.M)
    open(py, 'w').write(src)
    pkg = os.path.dirname(stl2cad.__file__)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.path.join(pkg, 'swmock'), os.path.dirname(pkg)]),
               SWMOCK_CUT_REVERSED='1', SWMOCK_BOSS_REVERSED='0', SWMOCK_EDGE_SIGHT='1')
    r = subprocess.run([sys.executable, py, out], capture_output=True, text=True, timeout=timeout, env=env)
    log = r.stdout + r.stderr
    if r.returncode != 0 or not os.path.exists(out):
        return False, log
    s = cq.importers.importStep(out).val()
    v, vr = mass_props(s)[0], mass_props(ref)[0]

    def diff(a, b):
        p, q = a.cut(b), b.cut(a)
        return abs(mass_props(p)[0]) + abs(mass_props(q)[0]), p.isValid() and q.isValid()
    dv, ok = diff(s, ref)
    if not (ok and dv < 1e-3 * vr):
        # OCC booleans of the two solids can fail (cut returns the whole solid, common returns nothing);
        # then the solids are compared by volume, surface area, centre of mass and inertia tensor
        bb = ref.BoundingBox()
        if invariants(s, ref, size=bb.DiagonalLength):
            return True, log + ('\nstand-in: volume %.3f, reference %.3f mm^3; booleans failed, volume, area, '
                                'centre of mass and inertia agree' % (v, vr))
    return (abs(v - vr) < 1e-3 * vr and dv < 1e-3 * vr), log + '\nstand-in: volume %.3f, reference %.3f mm^3, difference volume %.3f mm^3' % (v, vr, dv)
