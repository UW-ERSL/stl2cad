"""Writes bikehook_exact.py: a standalone CadQuery script of the bike hook (exact preserve regions, trunk loft,
arm members, keep-out cuts)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_TO = os.path.join(HERE, '..', '..', '..', 'to')
STL = os.path.join(HERE, '..', '..', 'bikehook.stl')

import sys, pickle, inspect, numpy as np, trimesh
import preserve as P, loft_bh, frame3d
sys.path.insert(0, REPO_TO)

def members(R, ztop=-30):
    J = np.array(R['joints']); out = []
    for m in R['members']:
        if max(J[m['i']][2], J[m['j']][2]) <= ztop:
            continue
        if R['kind'][m['i']] == 'hub' and R['kind'][m['j']] == 'hub':
            continue
        p0, p1, t, r0, r1 = frame3d.member_geom(R, m)
        out.append((tuple(round(float(v), 2) for v in p0), tuple(round(float(v), 2) for v in p1), round(float(r0), 2), round(float(r1), 2)))
    return out

def write(path, ARMS, secs):
    src = inspect.getsource(P)
    head, body = src.split('def fillet_polygon', 1)
    pres = 'def fillet_polygon' + body.split("if __name__ == '__main__':")[0]
    pres = pres.split('def inflated(')[0]
    consts = head.split('import cadquery as cq', 1)[1]
    L = []
    L.append('"""Bike hook (bikehook.stl), written by stl2cad.\n\nPreserve regions (exact): six bosses D%.0f x %.0f with D%.1f holes; the J cradle, one sweep of a\ntapered section (%.0f mm at the centre, %.0f mm at the sides) along lines and R%.0f arcs, with a rib under the bottom.\nConnector: a ruled loft through %d horizontal sections (%d corners each, R%.0f tangent fillets) and %d curved\narms (natural cubic splines; elliptical sections at spline points at most 7 mm apart; here ruled lofts between\nthem, in SolidWorks one loft along the spline). Cut last: above the boss tops, the tyre space inside the J, the holes.\n\npython bikehook_exact.py [out.step]\n"""' % (
        2 * P.R_BOSS, P.H_BOSS, 2 * P.R_HOLE, P.T, P.T - 2 * P.E, P.R_J, len(secs), loft_bh.K, loft_bh.R_F, len(ARMS)))
    L.append('import math, sys\nimport numpy as np\nimport cadquery as cq\nfrom cadquery import Vector\n' + consts)
    L.append(pres)
    L.append('\n# ---- connector ----\nR_F = %r\nSECTIONS = [  # (z, corners (x, y))\n' % loft_bh.R_F + ''.join(
        '    (%.1f, [%s]),\n' % (z, ', '.join('(%.1f, %.1f)' % tuple(v) for v in V)) for z, V in secs) + ']\n')
    L.append('ARMS = [  # spline points (natural cubic, chord-length parameters); sections (point index, half-axes a, b,\n'
             '          # major-axis direction)\n' + ''.join('    dict(points=%r,\n         sec=%r),\n' % (a['points'], a['sec']) for a in ARMS) + ']\n')
    import arms as _arms, frame3d as _f3
    L.append(inspect.getsource(_f3.natural_spline) + '\n\n' + inspect.getsource(_f3.arc_lengths) + '\n\n'
             + 'def spline(P):\n    return natural_spline(np.asarray(P, float))\n\n\n'
             + inspect.getsource(_arms.sample) + '\n\n' + inspect.getsource(_arms.ellipse_wire) + '\n\n'
             + inspect.getsource(_arms.arm_solid) + '\n\n'
             + inspect.getsource(_arms.arm_loft) + '\n')
    L.append(inspect.getsource(loft_bh.radii) + '\n\n' + inspect.getsource(loft_bh.section_wire).replace('r=R_F', 'r=R_F') + '\n')
    L.append('''
def trunk():
    return cq.Solid.makeLoft([section_wire(np.array(V), z) for z, V in SECTIONS], ruled=True)


def arm(a):
    """The arm: ruled lofts between the ellipses at its spline points (SolidWorks lofts the same sections with the
    spline as centerline, smooth)."""
    sec = [(i, a_, b_, np.array(m)) for i, a_, b_, m in a['sec']]
    return arm_loft(np.array(a['points'], float), sec)


def keepout():
    """The tyre's space inside the J."""
    a = math.radians(45)
    w = (cq.Workplane('XZ').moveTo(C_R[0] + R_J, Z_LEG + 20).lineTo(C_R[0] + R_J, C_R[1])
         .threePointArc((C_R[0] + R_J * math.cos(a), C_R[1] - R_J * math.sin(a)), (C_R[0], Z_FLAT))
         .lineTo(C_L[0], Z_FLAT)
         .threePointArc((C_L[0] - R_J * math.cos(a), C_L[1] - R_J * math.sin(a)), (C_L[0] - R_J, C_L[1]))
         .lineTo(C_L[0] - R_J, Z_LEG + 20).close())
    return w.extrude(40, both=True).val()


def build():
    part = end_fillets(cradle_with_rib()).val()
    for x, y in BOSSES:
        part = part.fuse(cq.Solid.makeCylinder(R_BOSS, H_BOSS, cq.Vector(x, y, 0)))
    core = trunk()                      # arms fused to the trunk first, then the cradle and bosses
    for a in ARMS:
        core = core.fuse(arm(a))
    part = part.fuse(core)
    part = part.cut(cq.Solid.makeBox(300, 300, 50, cq.Vector(-150, -150, H_BOSS)))
    part = part.cut(keepout())
    for x, y in BOSSES:
        part = part.cut(cq.Solid.makeCylinder(R_HOLE, 60, cq.Vector(x, y, -30)))
    return part.clean()


if __name__ == '__main__':
    part = build()
    print('valid %s, solids %d, faces %d, volume %.1f mm^3' % (part.isValid(), len(part.Solids()), len(part.Faces()), part.Volume()))
    cq.exporters.export(cq.Workplane().add(part), sys.argv[1] if len(sys.argv) > 1 else 'bikehook_exact.step')
''')
    open(path, 'w').write('\n'.join(L))

if __name__ == '__main__':
    ARMS = pickle.load(open(sys.argv[1], 'rb'))
    W = trimesh.load(STL)
    secs = [(z, np.round(V, 1)) for z, V in loft_bh.sections(W)]
    write(sys.argv[2], ARMS, secs)
