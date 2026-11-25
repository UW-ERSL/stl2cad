"""CadQuery reference of the swept-member model: hubs, one swept member per skeleton member (rectangle
w x d normal to the path, w and d linear from start to end), a vertical cylinder at each joint that is not
on a hub, the original's holes cut last.

python sweepref.py ORIGINAL_EXACT.py TO.stl NAME [--units in] [--detail 1-5]
"""
import argparse, math, sys
import numpy as np
import cadquery as cq
import trimesh
from cadquery import Vector

import hubweb, members, swplates
from stl2cad.core import mass_props


def wire3(path, z):
    es = []
    for e in path:
        P = [Vector(p[0], p[1], z) for p in e[1:]]
        es.append(cq.Edge.makeLine(*P) if e[0] == 'line' else cq.Edge.makeThreePointArc(*P))
    return cq.Wire.assembleEdges(es)


def tangent(path, end):
    e = path[0] if end == 0 else path[-1]
    if e[0] == 'line':
        d = np.subtract(e[2], e[1])
    else:
        c = np.array(members.S_circle(*e[1:])); p = np.array(e[1] if end == 0 else e[3])
        r = p - c; d = np.array([-r[1], r[0]])
        m = np.array(e[2]) - np.array(e[1])                  # orient along the arc's direction
        if end == 0 and d @ m < 0 or end == 1 and d @ (np.array(e[3]) - np.array(e[2])) < 0:
            d = -d
    return d / np.linalg.norm(d)


SHAPE = 'ellipse'        # section: 'rect', 'ellipse' or 'round' (rectangle with corner radius ROUND x the smaller side)
ROUND = 0.35


def section(p, t, w, d, z):
    """Section w (in-plane, across the path) x d (along z) centred at (p, z), normal t."""
    n = np.array([-t[1], t[0]])
    c = np.array([p[0], p[1], z])
    if SHAPE != 'rect':
        pl = cq.Plane(origin=Vector(*c), xDir=Vector(n[0], n[1], 0), normal=Vector(t[0], t[1], 0))
        wp = cq.Workplane(pl)
        if SHAPE == 'ellipse':
            return wp.ellipse(w / 2, d / 2).wires().val() if False else cq.Wire.makeEllipse(
                w / 2, d / 2, Vector(*c), Vector(t[0], t[1], 0), Vector(n[0], n[1], 0))
        return wp.rect(w, d).vertices().fillet2D(ROUND * min(w, d)).wires().val() if False else \
            wp.sketch().rect(w, d).vertices().fillet(ROUND * min(w, d)).finalize().val().outerWire() if False else \
            cq.Wire.makePolygon([Vector(*x) for x in _round_rect(c, n, w, d, ROUND * min(w, d))], close=True)
    a, b = np.r_[n * w / 2, 0], np.array([0, 0, d / 2])
    pts = [c + a + b, c - a + b, c - a - b, c + a - b]
    return cq.Wire.makePolygon([Vector(*q) for q in pts], close=True)


def _round_rect(c, n, w, d, r, k=6):
    """Polygon of a rectangle with rounded corners (k segments per corner), centre c, across n, up z."""
    zz = np.array([0, 0, 1.0]); nn = np.r_[n, 0]
    pts = []
    for sx, sy, a0 in ((1, 1, 0), (-1, 1, 90), (-1, -1, 180), (1, -1, 270)):
        cc = c + nn * sx * (w / 2 - r) + zz * sy * (d / 2 - r)
        for a in np.radians(a0 + np.linspace(0, 90, k + 1)):
            pts.append(cc + r * (math.cos(a) * nn + math.sin(a) * zz))
    return pts


def member_solid(m, z):
    path = m['path']
    p0, p1 = path[0][1], path[-1][-1]
    prof = [section(p0, tangent(path, 0), m['w0'], m['d0'], z), section(p1, tangent(path, 1), m['w1'], m['d1'], z)]
    return cq.Solid.sweep_multi(prof, wire3(path, z), makeSolid=True, isFrenet=False)


def build(R, hubs):
    z = R['mid']
    parts = []
    for h in hubs:
        parts.append(cq.Workplane().workplane(offset=h['lo']).center(*h['c']).circle(h['ro']).extrude(h['hi'] - h['lo']).val())
    lo_h, hi_h = min(h['lo'] for h in hubs), max(h['hi'] for h in hubs)
    ends = {}
    for m in R['members']:
        for k in ('d0', 'd1'):
            m[k] = min(m[k], 2 * min(z - lo_h, hi_h - z))
        parts.append(member_solid(m, z))
        for j, w, d in ((m['i'], m['w0'], m['d0']), (m['j'], m['w1'], m['d1'])):
            r, dd = ends.get(j, (0, 0)); ends[j] = (max(r, w / 2), max(dd, d))
    for j, (r, d) in ends.items():
        if j in R['hub_of']:
            continue
        # radius 5 % over the widest end: the end corners of that member would lie on the cylinder (a
        # degenerate contact that crashed OCC's Boolean)
        parts.append(cq.Workplane().workplane(offset=z - d / 2).center(*R['joints'][j]).circle(1.05 * r).extrude(d).val())
    # fuzzy Booleans (1e-4 mm): sweeps and cylinders meeting at a joint touch along nearly coincident faces
    est = sum(p.Volume() for p in parts)
    for tol in (1e-4, 1e-3, None, 1e-2):
        u = parts[0]
        for s in parts[1:]:
            u = u.fuse(s, tol=tol) if tol else u.fuse(s)
        u = u.clean()
        if len(u.Solids()) == 1 and u.isValid() and u.Volume() > 0.5 * est:
            break
        print('  fuse with tol %s: %d solid(s), volume %.0f of %.0f; trying again' % (tol, len(u.Solids()), u.Volume(), est))
    for h in hubs:
        u = u.cut(cq.Workplane().workplane(offset=h['lo'] - 1).center(*h['c']).circle(h['ri']).extrude(h['hi'] - h['lo'] + 2).val())
    return u


def compare(to, m, pitch=1.0):
    lo = np.minimum(to.bounds[0], m.bounds[0]) - 1
    shape = tuple(np.ceil((np.maximum(to.bounds[1], m.bounds[1]) + 1 - lo) / pitch).astype(int) + 1)

    def vox(mesh):
        g = mesh.voxelized(pitch).fill(); idx = np.argwhere(g.matrix)
        pts = trimesh.transformations.transform_points(idx, g.transform)
        G = np.zeros(shape, bool); k = np.round((pts - lo) / pitch).astype(int); G[k[:, 0], k[:, 1], k[:, 2]] = True
        return G
    A, B = vox(to), vox(m); both = (A & B).sum()
    return 100 * both / A.sum(), 100 * both / B.sum()


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('base'); ap.add_argument('to'); ap.add_argument('name')
    ap.add_argument('--units', default='in'); ap.add_argument('--detail', type=int, default=3)
    a = ap.parse_args()
    ns = swplates._namespace(a.base)
    to = trimesh.load(a.to); to.apply_scale({'m': 1000.0, 'mm': 1.0, 'in': 25.4}[a.units])
    hw = hubweb.detect(ns['original'], to)
    R = members.members(hw['shadow'], hw['hubs'], to, hw['axis'], a.detail)
    s = build(R, hw['hubs'])
    cq.exporters.export(cq.Workplane().add(s), a.name + '.stl', tolerance=0.02, angularTolerance=0.2)
    cq.exporters.export(cq.Workplane().add(s), a.name + '.step')
    m = trimesh.load(a.name + '.stl')
    narc = sum(e[0] == 'arc' for mm in R['members'] for e in mm['path'])
    print('%d members (%d arcs in the paths), valid %s, solids %d' % (len(R['members']), narc, s.isValid(), len(s.Solids())))
    print('volume: TO %.0f, model %.0f mm^3 (%+.1f %%)' % (to.volume, m.volume, 100 * (m.volume / to.volume - 1)))
    print('voxel overlap: %.0f %% of the TO in the model, %.0f %% of the model in the TO' % compare(to, m))
