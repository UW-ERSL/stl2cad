"""SolidWorks script (pywin32) for a TO result built as plates of the original part.

Each plate becomes one Extruded Boss from a sketch of its exact region (original section minus
the plate's cuts, as lines, arcs and circles). The first plate is built along its own axis; the
others with the body turned so their axis is z (stl2cad's Move/Copy Body sequence), each checked
against the volume and centre of mass of the CadQuery plate model. Re-entrant edges of the result
(face normals: BRepOffset_Analyse) are filleted last, picked at their midpoints; the radius is a
parameter of the script.
"""
import math
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


def concave_edges(solid):
    an = BRepOffset_Analyse(solid.wrapped, 0.01)
    return [e for e in solid.Edges()
            if any(i.Type() == ChFiDS_TypeOfConcavity.ChFiDS_Concave for i in an.Type(e.wrapped))]


def write(path, base_script, name, plates, box, r_fillet=1.0, nd=6, digits=3, detail=None):
    ns = _namespace(base_script)
    t3 = lambda q: '(%.5f, %.5f, %.5f)' % tuple(float(x) for x in q)
    solids = [ns['plate'](p['axis'], p['lo'], p['hi'], p['cuts'], box) for p in plates]
    order = sorted(range(len(plates)), key=lambda i: -mass_props(solids[i])[0])   # largest plate first
    k = plates[order[0]]['axis']
    R0, back = swgen.BUILD_FRAME[k]
    i0, j0 = [m for m in range(3) if m != k]
    union, steps = None, []
    for n, idx in enumerate(order):
        p, s = plates[idx], solids[idx]
        before = None if union is None else np.array(mass_props(union)[1])
        union = s if union is None else union.fuse(s).clean()
        vol, com = mass_props(union)
        steps.append(dict(p=p, idx=idx, loops=region_loops(s, p['axis'], p['lo']),
                          volume=vol, com=np.array(com), before=before))
    size = float(np.linalg.norm(np.array(box[1]) - np.array(box[0])))

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
         '# detail level %s (1: fewest dimensions ... 5: smallest model error)' % detail, '',
         'AXIS_I, AXIS_J = %d, %d           # model axes of the sketch coordinates (a, b) of the first plate' % (i0, j0),
         'PLANE_NAME = "Front Plane"', 'PLANE_INDEX = 0',
         'ROT = %s' % [list(map(float, r)) for r in R0],
         'ROT_BACK = %s' % (tuple(back),),
         'FINAL_CENTROID = %s   # after the first plate, STL frame, mm' % t3(steps[0]['com']),
         'SIZE = %s' % swgen._f(size, 3)]
    allp = np.array([q for st in steps for lp in st['loops'] for e in lp for q in e[1:] if np.ndim(q) == 1])
    L.append('HELPER_AB = (%s, %s)' % (swgen._f(allp[:, 0].min() - 10, 3), swgen._f(allp[:, 1].min() - 10, 3)))
    # ---- first plate: FEATURES -----------------------------------------------------------------
    st = steps[0]; p = st['p']
    plans, full = plans_of(st['loops'])
    cb = np.asarray(R0, float) @ st['com']
    L += ['', '# ---- first plate: sketch loops (a, b in mm), levels along its axis, expected volume (mm^3) ----',
          'FEATURES = [', '    dict(',
          '        name="Plate 1 %s %.4g to %.4g", op="add", lo=%s, hi=%s, volume=%.4f, centroid=%s, holes=False,'
          % (swgen.AX[p['axis']], p['lo'], p['hi'], swgen._f(p['lo'], nd), swgen._f(p['hi'], nd), st['volume'], t3(cb)),
          '        expect="%s",' % ('fully defined' if full else 'under defined'),
          '        plans=' + repr(swgen._round_plans(plans)) + ',', '        loops=[']
    for lp in st['loops']:
        L.append('        ' + swgen._loop_code(_Loop(lp), max(nd, 9)) + ',')
    L += ['        ],', '    ),', ']', '', 'CYLINDERS = []', 'REVOLVES = []', 'REV_TURN = None', 'REV_DONE = None', 'CUTS = []']
    # ---- other plates: XCUTS entries with op "add" -----------------------------------------------
    L += ['', '# ---- other plates: the body is turned so the plate axis is z, the plate added, the body turned back ----',
          'XCUTS = [']
    for n, st in enumerate(steps[1:], 2):
        p = st['p']; j = p['axis']
        m, ang = swgen.XCUT_TURN[j]
        R = swgen._rot(m, ang)
        u, v = [i for i in range(3) if i != j]
        def to_sk(q):
            P3 = np.zeros(3); P3[u], P3[v] = q
            return tuple(float(x) for x in (R @ P3)[:2])
        loops = []
        for lp in st['loops']:
            es = []
            for e in lp:
                if e[0] == 'circle':
                    es.append(('circle', to_sk(e[1]), e[2]))
                else:
                    es.append((e[0],) + tuple(to_sk(q) for q in e[1:]))
            loops.append(es)
        sgn = float((R @ np.eye(3)[j])[2])
        lo, hi = sorted((sgn * p['lo'], sgn * p['hi']))
        plans, full = plans_of(loops)
        P = np.array([q for lp in loops for e in lp for q in e[1:] if np.ndim(q) == 1])
        L += ['    dict(',
              '        name="Plate %d %s %.4g to %.4g", op="add", axis=%s, angle=%.12f, lo=%s, hi=%s, volume=%.4f,'
              % (n, swgen.AX[j], p['lo'], p['hi'], m, ang, swgen._f(lo, nd), swgen._f(hi, nd), st['volume']),
              '        turned=%s, centroid=%s, back=%s, holes=False,' % (t3(R @ st['before']), t3(R @ st['com']), t3(st['com'])),
              '        helper=(%s, %s),' % (swgen._f(P[:, 0].min() - 10, 3), swgen._f(P[:, 1].min() - 10, 3)),
              '        expect="%s",' % ('fully defined' if full else 'under defined'),
              '        plans=' + repr(swgen._round_plans(plans)) + ',', '        loops=[']
        for lp in loops:
            L.append('        ' + swgen._loop_code(_Loop(lp), max(nd, 9)) + ',')
        L += ['        ],', '    ),']
    L += [']', 'BLENDS = []']
    # ---- re-entrant edges: picked at their midpoints, with the volume each fillet adds ----------
    E = concave_edges(union)
    L += ['', '# ---- re-entrant edges (STL frame, mm): midpoint, length; filleted last ----', 'JUNCTIONS = [']
    for e in E:
        L.append('    (%s, %.4f),' % (t3(e.positionAt(0.5).toTuple()), e.Length()))
    L += [']', '']
    rt = swgen._RUNTIME
    head, main = rt.split('if __name__ == "__main__":')
    L.append(head + FILLET_CODE)
    open(path, 'w').write('\n'.join(L) + '\n')
    return dict(volume=steps[-1]['volume'], concave=len(E), union=union)


FILLET_CODE = r'''
def junction_fillets(sw, r):
    """Fillet the re-entrant edges with radius r: all in one Fillet feature, else edge by edge.
    A fillet is kept when the volume grows by about r^2 (1 - pi/4) per mm of edge (the exact
    amount for two perpendicular plane faces); SolidWorks' failures are reported, not fatal."""
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

    def attempt(edges, name):
        v0 = sw.mass()[0]
        est = sum(L for _, L in edges) * r * r * (1 - math.pi / 4)
        for mark in (1, 0):
            for how, mk in makes:
                sw.model.ClearSelection2(True)
                if not all(ext.SelectByID2("", "EDGE", *(c * MM for c in sw.in_model_frame(p)), True, mark, NOTHING, 0)
                           for p, _ in edges):
                    continue
                f = None
                try:
                    f = mk()
                except Exception as ex:
                    print(f"  {name}: {how} failed: {type(ex).__name__}: {ex}")
                if f is None:
                    continue
                dv = sw.mass()[0] - v0
                if dv > 0 and abs(dv - est) <= 0.6 * est + 1e-3:
                    f.Name = name
                    return True
                print(f"  {name}: {how} added {dv:.3f} mm^3, expected about {est:.3f}; removed")
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
    junction_fillets(sw, R_FILLET)
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
               SWMOCK_CUT_REVERSED='1', SWMOCK_BOSS_REVERSED='0')
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
