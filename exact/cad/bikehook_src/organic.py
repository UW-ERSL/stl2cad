"""Organic connector of the bike hook: the mesh minus the cradle, fitted as a 3D frame (to/frame3d.py) whose hubs
are the six bosses."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_TO = os.path.join(HERE, '..', '..', '..', 'to')
STL = os.path.join(HERE, '..', '..', 'bikehook.stl')

import sys, time
import numpy as np, trimesh, manifold3d as mf, cadquery as cq
sys.path.insert(0, REPO_TO)
import frame3d
import preserve as P


def to_trimesh(shape, tol=0.02):
    vs, fs = shape.tessellate(tol, 0.1)
    return trimesh.Trimesh(np.array([v.toTuple() for v in vs]), np.array(fs))


def to_manifold(t):
    return mf.Manifold(mf.Mesh(vert_properties=np.asarray(t.vertices, np.float32), tri_verts=np.asarray(t.faces, np.uint32)))


def from_manifold(M):
    g = M.to_mesh(); return trimesh.Trimesh(np.asarray(g.vert_properties)[:, :3], np.asarray(g.tri_verts))


def original():
    """Stand-in for frame3d's original part: a block round the whole hook, the bosses' holes through its top
    5 mm (so each hole is a hub of the boss height)."""
    lo = cq.Workplane('XY').box(140, 120, 100, centered=(True, True, False)).translate((30, 0, -100))
    top = cq.Workplane('XY').box(140, 120, P.H_BOSS, centered=(True, True, False)).translate((30, 0, 0))
    for x, y in P.BOSSES:
        top = top.cut(cq.Workplane('XY').center(x, y).circle(P.R_HOLE).extrude(P.H_BOSS))
    return lo.union(top).val()


def web_mesh(mesh, cradle):
    return from_manifold(to_manifold(mesh) - to_manifold(to_trimesh(cradle)))


if __name__ == '__main__':
    t0 = time.time()
    m = trimesh.load(STL)
    cr = P.inflated(float(sys.argv[2]) if len(sys.argv) > 2 else 0.6).val()
    W = web_mesh(m, cr)
    print('web: %d faces, %.0f mm^3 (mesh %.0f, cradle %.0f) %.1fs' % (len(W.faces), W.volume, m.volume, cr.Volume(), time.time() - t0))
    W.export('web.stl')
    detail = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    frame3d.REACH = float(sys.argv[3]) if len(sys.argv) > 3 else frame3d.REACH
    R = frame3d.fit(original(), W, detail=detail)
    import pickle; pickle.dump(R, open('frame_d%d.pkl' % detail, 'wb'))
    J = np.array(R['joints'])
    for n, mm in enumerate(R['members']):
        a, b = J[mm['i']], J[mm['j']]
        print('%2d  %s -> %s  L %.1f  r %.2f-%.2f' % (n, np.round(a, 1), np.round(b, 1), np.linalg.norm(b - a), mm['r0'], mm['r1']))
    print('%.1fs' % (time.time() - t0))
