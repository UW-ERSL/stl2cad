"""TO reconstruction as plates of the original part (case 2: original + TO).

python run_plates.py ORIGINAL.stl TO.stl NAME [--units m|mm|in] [--detail 1-5] [--sw] [--plot]
Writes NAME.py (CadQuery reference script), NAME.step, with --sw NAME_solidworks.py (checked on the
SolidWorks API stand-in), with --plot NAME_sections.png. The original's exact model is built once, in
out_<original>/ under the current folder.
"""
import argparse, os, subprocess, sys, time
import numpy as np, trimesh
import plates, codegen3d
from tofit import Section

ap = argparse.ArgumentParser()
ap.add_argument('original'); ap.add_argument('to'); ap.add_argument('name')
ap.add_argument('--units', default='m'); ap.add_argument('--pitch', type=float, default=0.5)
ap.add_argument('--plot', action='store_true')
ap.add_argument('--fillet', type=float, default=1.0, help='radius on re-entrant edges (mm)')
ap.add_argument('--sw', action='store_true', help='SolidWorks script NAME_solidworks.py, checked on the API stand-in')
ap.add_argument('--detail', type=int, default=4, choices=range(1, 6),
                help='1: fewest dimensions (larger model error) ... 5: smallest model error')
ap.add_argument('--radii', type=int, default=None, help='global corner radii (0: radius per corner); default from --detail')
ap.add_argument('--fit-only', action='store_true')
a = ap.parse_args()
import tofit
D = tofit.set_detail(a.detail)
if a.radii is None:
    a.radii = D['radii']
scale = {'m': 1000.0, 'mm': 1.0, 'in': 25.4}[a.units]
T = {}; t0 = time.time()

# ---- 1. original part: exact features from stl2cad (once per original) ----------
stem = os.path.splitext(os.path.basename(a.original))[0]
base = os.path.join('out_' + stem, stem + '_exact.py')
if not os.path.exists(base):
    subprocess.run([sys.executable, '-m', 'stl2cad.cli', a.original, '--units', a.units, '--no-solidworks',
                    '-o', 'out_' + stem], check=True, stdout=subprocess.DEVNULL)
T['original'] = time.time() - t0

# ---- 2. meshes, plates -----------------------------------------------------------
t = time.time()
orig = trimesh.load(a.original); orig.apply_scale(scale)
to = trimesh.load(a.to); to.apply_scale(scale)
bodies = to.split(only_watertight=False)
to = max(bodies, key=lambda b: len(b.faces))
print('TO: %d bodies, kept the largest (%d facets, %.0f mm^3)' % (len(bodies), len(to.faces), to.volume))
P = plates.detect(to, a.pitch, bounds=orig.bounds)
T['plates'] = time.time() - t

# ---- 3. one section fit per plate ------------------------------------------------
t = time.time()
for p in P:
    lvl = 0.5 * (p['lo'] + p['hi'])
    sec = Section(plates.section(orig, p['axis'], lvl), plates.section(to, p['axis'], lvl))
    p['cuts'] = sec.fit(); p['stats'] = sec.stats; p['sec'] = sec; p['refs'] = sec.refs
    s = sec.stats
    print('  plate normal to %s at %.2f: %d cuts, %d segments, mean %.3f max %.3f mm%s%s'
          % ('xyz'[p['axis']], lvl, len(p['cuts']), s['segments'], s['mean'], s['max'],
             ', refit finely: %s' % s['fine'] if s['fine'] else '',
             ', %d separate island(s) of TO material ignored' % s['islands'] if s['islands'] else ''))
if a.radii:
    import radii
    RQ = radii.quantize(P, a.radii, accept=D['accept'])
T['fit'] = time.time() - t
if a.fit_only:
    print('fit %.1f s' % T['fit']); sys.exit()

# ---- 4. script, model ------------------------------------------------------------
t = time.time()
codegen3d.write(a.name + '.py', base, a.name, P, orig.bounds, 0.0)   # internal reference: no 3D fillets (OCC)
subprocess.run([sys.executable, a.name + '.py', a.name + '.step'], check=True)
T['model'] = time.time() - t
if a.sw:
    t = time.time()
    import swplates
    info = swplates.write(a.name + '_solidworks.py', base, a.name + '_solidworks', P, orig.bounds, a.fillet,
                          digits=D['digits'], detail=a.detail)
    print('SolidWorks script: %s_solidworks.py (%d plate(s), %d re-entrant edge(s) to fillet)' % (a.name, len(P), info['concave']))
    ok, log = swplates.verify(a.name + '_solidworks.py', info['union'])
    open(a.name + '_standin.log', 'w').write(log)
    print('  API stand-in: %s; %s' % ('ok' if ok else 'FAILED', log.strip().splitlines()[-1]))
    T['solidworks'] = time.time() - t

# ---- 5. checks -------------------------------------------------------------------
t = time.time()
import cadquery as cq
s = cq.importers.importStep(a.name + '.step').val()
cq.exporters.export(cq.Workplane().add(s), a.name + '.stl', tolerance=0.02, angularTolerance=0.2)
m = trimesh.load(a.name + '.stl')
print('model: valid %s, solids %d, faces %d' % (s.isValid(), len(s.Solids()), len(s.Faces())))
ns = {}; exec(compile('from OCP.BRepOffset import BRepOffset_Analyse\nfrom OCP.ChFiDS import ChFiDS_TypeOfConcavity\n'
     + codegen3d.HELPERS[codegen3d.HELPERS.index('def concave_edges'):codegen3d.HELPERS.index('def fillet_reentrant')], 'chk', 'exec'), ns)
left = ns['concave_edges'](s)
print('re-entrant edges (filleted in SolidWorks): %d%s' % (len(left), ''.join('\n  at (%.1f, %.1f, %.1f)' % e.positionAt(0.5).toTuple() for e in left[:5])))
print('volume: TO %.0f, model %.0f mm^3 (%+.1f %%)' % (to.volume, m.volume, 100 * (m.volume / to.volume - 1)))
# voxel overlap on a common grid
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
for p in P:
    st = Section.thin_strips(p['sec'], p['sec'].material(p['cuts']))
    inTO = Section.thin_strips(p['sec'], p['sec'].M)
    bad = [g for g in st if not any(g.distance(h) < 0.5 for h in inTO)]
    if bad:
        print('  plate normal to %s: %d material strip(s) thinner than 1 mm' % ('xyz'[p['axis']], len(bad)))
T['checks'] = time.time() - t
print('time: ' + ', '.join('%s %.1f s' % kv for kv in T.items()) + ', total %.1f s' % (time.time() - t0))

if a.plot:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(P), figsize=(6 * len(P), 4))
    axs = np.atleast_1d(axs)
    for ax, p in zip(axs, P):
        sec = p['sec']; lvl = 0.5 * (p['lo'] + p['hi'])
        def draw(g, **k):
            for q in (g.geoms if hasattr(g, 'geoms') else [g]):
                if q.is_empty: continue
                ax.plot(*q.exterior.xy, **k)
                for h in q.interiors: ax.plot(*h.xy, **k)
        draw(sec.O, color='k', lw=0.8, ls='--'); draw(sec.M, color='0.55', lw=2.2)
        draw(plates.section(m, p['axis'], lvl), color='r', lw=0.9)
        u, v = plates.frame(p['axis'])
        ax.set_aspect('equal'); ax.set_xlabel('xyz'[u]); ax.set_ylabel('xyz'[v])
        ax.set_title('plate normal to %s, %s = %.1f to %.1f' % ('xyz'[p['axis']], 'xyz'[p['axis']], p['lo'], p['hi']))
    plt.tight_layout(); plt.savefig(a.name + '_sections.png', dpi=90)
