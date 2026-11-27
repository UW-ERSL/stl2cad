"""Writes a SolidWorks probe for the bike hook from ../bikehook_exact.py.

python gen_probe.py [TEMPLATE] [OUT]   (from exact/cad/bikehook_src; default probe_template.py -> ../sw_probe_hook.py;
                                       probe_template2.py ../sw_probe_hook2.py for probe 2)
"""
import importlib.util, math, os, numpy as np, cadquery as cq
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, '..')
spec = importlib.util.spec_from_file_location('bh', os.path.join(OUT, 'bikehook_exact.py')); bh = importlib.util.module_from_spec(spec); spec.loader.exec_module(bh)

def r3(p): return tuple(round(float(v), 6) for v in p)

def ring(V, z, radii):
    tg = bh.fillet_polygon(V, radii); ents = []
    n = len(V)
    for i in range(n):
        t1, m, t2 = tg[i]
        ents.append(('A', r3((*t1, z)), r3((*m, z)), r3((*t2, z))))
        nx = tg[(i + 1) % n][0]
        if np.linalg.norm(nx - t2) > 1e-6:
            ents.append(('L', r3((*t2, z)), r3((*nx, z))))
    return ents

# J path (model xz plane, y = 0)
a = math.radians(45); r = bh.R_J
P0 = (bh.C_R[0] + r, 0, bh.Z_LEG); P1 = (bh.C_R[0] + r, 0, bh.C_R[1])
PATH = [('L', P0, P1),
        ('A', P1, (bh.C_R[0] + r * math.cos(a), 0, bh.C_R[1] - r * math.sin(a)), (bh.C_R[0], 0, bh.Z_FLAT)),
        ('L', (bh.C_R[0], 0, bh.Z_FLAT), (bh.C_L[0], 0, bh.Z_FLAT)),
        ('A', (bh.C_L[0], 0, bh.Z_FLAT), (bh.C_L[0] - r * math.cos(a), 0, bh.C_L[1] - r * math.sin(a)), (bh.C_L[0] - r, 0, bh.C_L[1])),
        ('L', (bh.C_L[0] - r, 0, bh.C_L[1]), (bh.C_L[0] - r, 0, bh.Z_LIP))]
PATH = [tuple([e[0]] + [r3(p) for p in e[1:]]) for e in PATH]
# profile at the leg's top end (plane z = Z_LEG), local u = y, v = +x outward from x = 20
u0 = bh.W - bh.E / bh.SLOPE
loc = [(-u0, 0), (u0, 0), (bh.W, bh.E), (bh.W, bh.T - bh.E), (u0, bh.T), (-u0, bh.T), (-bh.W, bh.T - bh.E), (-bh.W, bh.E)]
rad = [bh.R_BLEND, bh.R_BLEND, bh.R_SIDE, bh.R_SIDE, bh.R_BLEND, bh.R_BLEND, bh.R_SIDE, bh.R_SIDE]
tg = bh.fillet_polygon(loc, rad)
X0 = bh.C_R[0] + bh.R_J
to3 = lambda q: r3((X0 + q[1], q[0], bh.Z_LEG))
PROFILE = []
for i in range(len(loc)):
    t1, m, t2 = tg[i]
    PROFILE.append(('A', to3(t1), to3(m), to3(t2)))
    nx = tg[(i + 1) % len(loc)][0]
    if np.linalg.norm(nx - t2) > 1e-9:
        PROFILE.append(('L', to3(t2), to3(nx)))
# rib outline (xz plane): read back from the rib solid's face at y = +RIB_W
rib = bh.rib().val()
f = [f for f in rib.Faces() if f.geomType() == 'PLANE' and abs(f.Center().y - bh.RIB_W) < 1e-6][0]
RIB = []
for e in f.outerWire().Edges():
    p0, p1 = e.startPoint(), e.endPoint()
    q = lambda v: r3((v.x, 0.0, v.z))
    if e.geomType() == 'LINE': RIB.append(('L', q(p0), q(p1)))
    else: RIB.append(('A', q(p0), q(e.positionAt(0.5)), q(p1)))
# loft sections
SECS = [(z, ring(np.array(V), z, bh.radii(np.array(V), bh.R_F))) for z, V in bh.SECTIONS]
w = [bh.section_wire(np.array(V), z) for z, V in bh.SECTIONS]
pair = cq.Solid.makeLoft(w[:2], ruled=True)
lofts = [cq.Solid.makeLoft(w[k:k + 2], ruled=True) for k in range(len(w) - 1)]
allr = lofts[0]
for L in lofts[1:]: allr = allr.fuse(L)
cr = bh.cradle().val(); cr2 = bh.cradle().union(bh.rib()).val()
cr3 = bh.cradle_with_rib().val(); cr4 = bh.end_fillets(bh.cradle_with_rib()).val()
E = dict(sweep=(cr.Volume(), cr.Center().toTuple()), rib=(cr2.Volume(), cr2.Center().toTuple()),
         ribfil=(cr3.Volume(), cr3.Center().toTuple()), ends=(cr4.Volume(), cr4.Center().toTuple()),
         pair=(pair.Volume(), pair.Center().toTuple()), lofts=(allr.Volume(), allr.Center().toTuple()))
crr = bh.cradle().union(bh.rib()).val()
es = [e for f in crr.Faces() if f.geomType() == 'PLANE' and abs(abs(f.Center().y) - bh.RIB_W) < 1e-6 and f.Center().z < -80
      for e in f.Edges() if (e.geomType() == 'LINE' and abs(e.Center().z - bh.RIB_Z) < 1e-3) or (e.geomType() == 'CIRCLE' and e.Center().z < -85.5)]
r1 = cq.Workplane().add(crr).newObject(es).fillet(1.0).val()
E['ribfil1'] = (r1.Volume(), r1.Center().toTuple())
E = {k: (round(v, 3), r3(c)) for k, (v, c) in E.items()}
E['pairs'] = [(round(L.Volume(), 3), r3(L.Center().toTuple())) for L in lofts]
fmt = lambda L: '[\n' + ''.join('    %r,\n' % (e,) for e in L) + ']'
import sys
TPL = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'probe_template.py')
OUTF = sys.argv[2] if len(sys.argv) > 2 else os.path.join(OUT, 'sw_probe_hook.py')
code = open(TPL).read()
code = code.replace('@PATH@', fmt(PATH)).replace('@PROFILE@', fmt(PROFILE)).replace('@RIB@', fmt(RIB))
code = code.replace('@SECS@', '[\n' + ''.join('    (%r, %s),\n' % (z, fmt(ents).replace('\n', '\n    ')) for z, ents in SECS) + ']')
code = code.replace('@RIB_PTS@', '[(0.0, 13.0, -90.0), (-14.0, 13.0, -87.320508), (13.234544, 13.0, -88.371701), (0.0, -13.0, -90.0), (-14.0, -13.0, -87.320508), (13.234544, -13.0, -88.371701)]')
code = code.replace('@EXPECT@', repr(E)).replace('@Z_LEG@', repr(bh.Z_LEG)).replace('@Z_LIP@', repr(bh.Z_LIP))
code = code.replace('@RIB_W@', repr(bh.RIB_W)).replace('@RIB_Z@', repr(bh.RIB_Z)).replace('@RIB_F@', repr(bh.RIB_F)).replace('@END_R@', repr(bh.END_R))
open(OUTF, 'w').write(code)
print(E)
