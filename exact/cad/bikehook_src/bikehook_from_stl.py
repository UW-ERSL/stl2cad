"""Bike hook (exact/bikehook.stl): CadQuery script and STEP of the generative-design hook.

Preserve regions (preserve.py): bosses and the J cradle with its rib, exact; their dimensions were measured from
the mesh's CAD faces (planar patches, ray casts and circle fits), not detected by the program.
Connector (automatic): the mesh minus the cradle grown by 0.6 mm (organic.py); the six arms (arms.py) from the
bars of a 3D frame (to/frame3d.py, detail 5, hubs = the bosses, hub reach 5 mm, members rising above z = -30), each
a natural cubic spline fitted to the centroids of the mesh sections normal to it, with elliptical sections of the
same area and principal directions; the trunk as a ruled loft through 11 horizontal sections z = -28 .. -85 of the mesh at x >= 16 (loft_bh.py;
not -22, where SolidWorks refuses the loft), each reduced to 12 corners (to/tofit.reduce_sides) with R2 tangent fillets, corners matched between sections.

python bikehook_from_stl.py [detail]      (from exact/cad/bikehook_src; writes ../bikehook_exact.py, .step, .stl)
"""
import os, sys, time, subprocess
import numpy as np, trimesh
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'to'))
import preserve as P, organic, loft_bh, gen_cq, frame3d, arms
from evalmodel import compare
from organic import to_trimesh

if __name__ == '__main__':
    detail = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t0 = time.time()
    mesh = trimesh.load(organic.STL)
    W = organic.web_mesh(mesh, P.inflated(0.6).val())
    frame3d.REACH = 5.0
    R = frame3d.fit(organic.original(), W, detail=detail)
    print('frame %.0f s' % (time.time() - t0))
    A = arms.fit_from_frame(mesh, R, P.BOSSES, P.R_BOSS, detail=detail)
    print('arms %.0f s' % (time.time() - t0))
    secs = [(z, np.round(V, 1)) for z, V in loft_bh.sections(mesh)]
    out = os.path.join(HERE, '..', 'bikehook_exact')
    gen_cq.write(out + '.py', A, secs)
    subprocess.run([sys.executable, out + '.py', out + '.step'], check=True)
    import cadquery as cq
    s = cq.importers.importStep(out + '.step').val()
    model = to_trimesh(s, 0.02); model.export(out + '.stl')
    c = compare(mesh, model)
    print('volume: mesh %.0f, model %.0f mm^3 (%+.1f %%)' % (c['vol_mesh'], c['vol_model'], 100 * (c['vol_model'] / c['vol_mesh'] - 1)))
    print('overlap: %.1f %% of the mesh in the model, %.1f %% of the model in the mesh' % (100 * c['mesh_in_model'], 100 * c['model_in_mesh']))
    print('distance mesh -> model: median %.2f, p95 %.2f, max %.2f mm' % c['mesh_to_model'])
    print('distance model -> mesh: median %.2f, p95 %.2f, max %.2f mm' % c['model_to_mesh'])
    print('total %.0f s' % (time.time() - t0))
