"""Mesh-level assembly of preserve + members and its comparison with the input mesh."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_TO = os.path.join(HERE, '..', '..', '..', 'to')
STL = os.path.join(HERE, '..', '..', 'bikehook.stl')

import sys, pickle, numpy as np, trimesh, manifold3d as mf, cadquery as cq
sys.path.insert(0, REPO_TO); sys.path.insert(0, HERE)
import frame3d, preserve as P
from organic import to_trimesh, to_manifold, from_manifold


def keepout():
    """The tyre's space inside the J (no material there): the J's inside, x between the legs, above the bottom,
    across the full width."""
    w = cq.Workplane('XZ').moveTo(P.C_R[0] + P.R_J, P.Z_LEG + 20).lineTo(P.C_R[0] + P.R_J, P.C_R[1]) \
        .threePointArc((P.C_R[0] + P.R_J * 0.7071, P.C_R[1] - P.R_J * 0.7071), (P.C_R[0], P.Z_FLAT)) \
        .lineTo(P.C_L[0], P.Z_FLAT) \
        .threePointArc((P.C_L[0] - P.R_J * 0.7071, P.C_L[1] - P.R_J * 0.7071), (P.C_L[0] - P.R_J, P.C_L[1])) \
        .lineTo(P.C_L[0] - P.R_J, P.Z_LEG + 20).close()
    return w.extrude(40, both=True).val()


def assemble(R, members=True):
    pres = P.end_fillets(P.cradle_with_rib()).val()
    M = to_manifold(to_trimesh(pres, 0.01))
    for x, y in P.BOSSES:
        M = M + to_manifold(to_trimesh(cq.Solid.makeCylinder(P.R_BOSS, P.H_BOSS, cq.Vector(x, y, 0)), 0.01))
    if members:
        for m in R['members']:
            M = M + to_manifold(to_trimesh(frame3d.member_solid(R, m), 0.01))
    # keep-outs: above the boss tops, the tyre space, the holes
    M = M - to_manifold(to_trimesh(cq.Solid.makeBox(300, 300, 50, cq.Vector(-150, -150, P.H_BOSS)), 0.01))
    M = M - to_manifold(to_trimesh(keepout(), 0.01))
    for x, y in P.BOSSES:
        M = M - to_manifold(to_trimesh(cq.Solid.makeCylinder(P.R_HOLE, 50, cq.Vector(x, y, -20)), 0.005))
    return from_manifold(M)


def compare(mesh, model, pitch=0.5):
    out = dict(vol_mesh=mesh.volume, vol_model=model.volume)
    A = to_manifold(mesh); B = to_manifold(model)
    i = (A ^ B).volume()
    out['mesh_in_model'] = i / A.volume(); out['model_in_mesh'] = i / B.volume()
    p = mesh.sample(30000); q = model.sample(30000)
    sd = lambda M, X: np.concatenate([np.abs(trimesh.proximity.signed_distance(M, X[i:i + 1000])) for i in range(0, len(X), 1000)])
    d1 = sd(model, p); d2 = sd(mesh, q)
    out['mesh_to_model'] = (np.median(d1), np.percentile(d1, 95), d1.max())
    out['model_to_mesh'] = (np.median(d2), np.percentile(d2, 95), d2.max())
    return out


if __name__ == '__main__':
    R = pickle.load(open(sys.argv[1], 'rb'))
    mesh = trimesh.load(STL)
    model = assemble(R)
    model.export(sys.argv[1].replace('.pkl', '_model.stl'))
    for k, v in compare(mesh, model).items():
        print(k, np.round(v, 3) if not isinstance(v, float) else round(v, 3))
