"""Trunk of the bike hook's connector as a loft through horizontal sections: each section is the main region of
the connector there, reduced to K corners (to/tofit.reduce_sides) with a tangent fillet at every corner."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_TO = os.path.join(HERE, '..', '..', '..', 'to')
STL = os.path.join(HERE, '..', '..', 'bikehook.stl')

import sys, math, numpy as np, trimesh, cadquery as cq
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
sys.path.insert(0, REPO_TO); sys.path.insert(0, HERE)
import plates, tofit
from preserve import fillet_polygon

ZS = [-28.0, -34.0, -40.0, -46.0, -52.0, -58.0, -64.0, -70.0, -76.0, -81.0, -85.0]
K, R_F = 12, 2.0
RULED = True


X_MIN = 16.0      # the trunk lies at x >= X_MIN (the leg side); sections taken from the whole mesh there


def main_region(W, z):
    from shapely.geometry import box
    s = plates.section(W, 2, z).intersection(box(X_MIN, -60, 120, 60))
    s = s.buffer(0.5).buffer(-1.0).buffer(0.5)
    g = list(s.geoms) if hasattr(s, 'geoms') else [s]
    return orient(max(g, key=lambda p: p.area))


def k_gon(P, k):
    V = np.array(P.simplify(0.3).exterior.coords)[:-1]
    tofit.CUT_WEIGHT = getattr(tofit, 'CUT_WEIGHT', 1.0)
    V = tofit.reduce_sides(list(V), k, P)
    V = np.array(V)
    if Polygon(V).exterior.is_ccw is False:
        V = V[::-1]
    return V


def align(V, prev):
    """The rotation of the ring closest to the previous section's ring (sum of corner distances)."""
    costs = [np.linalg.norm(np.roll(V, -i, axis=0) - prev, axis=1).sum() for i in range(len(V))]
    return np.roll(V, -int(np.argmin(costs)), axis=0)


def radii(V, r):
    """Fillet radius per corner, reduced so each fillet's tangent length is at most 0.45 x its sides."""
    n = len(V); out = []
    for i in range(n):
        a, b, c = V[i - 1], V[i], V[(i + 1) % n]
        u, w = (a - b) / np.linalg.norm(a - b), (c - b) / np.linalg.norm(c - b)
        phi = math.acos(np.clip(u @ w, -1, 1))
        tmax = 0.45 * min(np.linalg.norm(a - b), np.linalg.norm(c - b))
        out.append(min(r, tmax * math.tan(phi / 2)))
    return out


def section_wire(V, z, r=R_F):
    tg = fillet_polygon(V, radii(V, r))
    P3 = lambda p: cq.Vector(float(p[0]), float(p[1]), float(z))
    es = []
    for i in range(len(V)):
        t1, m, t2 = tg[i]
        es.append(cq.Edge.makeThreePointArc(P3(t1), P3(m), P3(t2)))
        n1 = tg[(i + 1) % len(V)][0]
        if np.linalg.norm(n1 - t2) > 1e-6:
            es.append(cq.Edge.makeLine(P3(t2), P3(n1)))
    return cq.Wire.assembleEdges(es)


def sections(W, zs=ZS, k=K):
    out, ref = [], None
    for z in zs:
        V = k_gon(main_region(W, z), k)
        V = V if ref is None else align(V, ref)
        ref = V
        out.append((z, V))
    return out


def loft(secs, r=R_F):
    return cq.Solid.makeLoft([section_wire(V, z, r) for z, V in secs], ruled=RULED)


if __name__ == '__main__':
    W = trimesh.load(STL)
    S = sections(W)
    for z, V in S: print(z, np.round(V, 1).tolist())
    L = loft(S)
    print('loft valid', L.isValid(), round(L.Volume()))
    cq.exporters.export(cq.Workplane().add(L), 'loft.step')
    vs, fs = L.tessellate(0.05, 0.1)
    trimesh.Trimesh(np.array([v.toTuple() for v in vs]), np.array(fs)).export('loft.stl')
