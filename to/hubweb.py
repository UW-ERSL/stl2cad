"""Hubs and web: a TO result made of rings around the original part's holes (hubs) joined by a flat band of
material (web) that does not fill plates of the original (round bars of a truss, an arm or a link).

detect(original, to_mesh) finds:
  hubs: the original's cylindrical holes (all along one axis) that the TO surrounds with a ring; each hub
        is the original within radius ro of the hole's axis, over the hole's extent along the axis, with
        ro the median outer radius of the TO ring (taken at the level along the hub where the TO ring is
        least joined to other material).
  web:  the band [lo, hi] along that axis holding the TO material away from the hubs (median bottom and
        top of the TO columns there), and the TO's shadow along the axis (its outline in the plane).
The web is then fitted like a plate whose TO section is the shadow (voids and notches as cuts).
"""
import math
import numpy as np
import trimesh
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Cylinder

import plates

MIN_RING = 0.6       # a hub needs TO material all round its hole over this fraction of the angles


def holes(original):
    """Cylindrical holes of the original (material outside the cylinder, closed all round), grouped by axis line:
    list of dict(axis, c=(u, v), r, lo, hi)."""
    out = []
    for f in original.Faces():
        ad = BRepAdaptor_Surface(f.wrapped)
        if ad.GetType() != GeomAbs_Cylinder:
            continue
        cyl = ad.Cylinder(); L = cyl.Axis().Location(); D = cyl.Axis().Direction()
        d = np.array([D.X(), D.Y(), D.Z()])
        k = int(np.argmax(np.abs(d)))
        if abs(abs(d[k]) - 1) > 1e-6:
            continue
        P = np.array([L.X(), L.Y(), L.Z()]); r = cyl.Radius()
        # a hole: a point just inside the cylinder (towards its axis) is outside the material
        bb = f.BoundingBox()
        lo_k, hi_k = (bb.xmin, bb.ymin, bb.zmin)[k], (bb.xmax, bb.ymax, bb.zmax)[k]
        u, v = [i for i in range(3) if i != k]
        q = P.copy(); q[k] = 0.5 * (lo_k + hi_k); q[u] += 0.5 * r
        if original.isInside(__import__('cadquery').Vector(*q), 1e-6):
            continue
        span = ad.LastUParameter() - ad.FirstUParameter()
        for h in out:
            if h['axis'] == k and math.hypot(h['c'][0] - P[u], h['c'][1] - P[v]) < 1e-3:
                h['lo'], h['hi'] = min(h['lo'], lo_k), max(h['hi'], hi_k)
                h['r'] = min(h['r'], r)
                h['span'] += span
                break
        else:
            out.append(dict(axis=k, c=(float(P[u]), float(P[v])), r=float(r), lo=float(lo_k), hi=float(hi_k), span=span))
    # a hole closes all round (its faces on one axis line cover 360 deg); a partial cylinder with material outside
    # it is the fillet of a concave corner (the Knuckle's slot), not a hole
    return [{k_: v_ for k_, v_ in h.items() if k_ != 'span'} for h in out if h['span'] >= 2 * math.pi - 1e-6]


def _ring_radius(to, h, nlev=9):
    """Outer radius of the TO ring around hole h: at nlev levels along the hole, the TO section's piece
    around the hole; the level where that piece is smallest (least joined to bars), its exterior's
    median distance from the axis. None if the TO does not surround the hole."""
    k, (cu, cv) = h['axis'], h['c']
    best = None
    for t in np.linspace(0.1, 0.9, nlev):
        lvl = h['lo'] + t * (h['hi'] - h['lo'])
        sec = plates.section(to, k, lvl)
        for g in getattr(sec, 'geoms', [sec]):
            if g.is_empty or not g.contains(Point(cu + h['r'] + 0.5, cv)) and not any(
                    Polygon(i).contains(Point(cu, cv)) for i in g.interiors):
                continue
            ext = np.asarray(g.exterior.coords)
            ang = np.arctan2(ext[:, 1] - cv, ext[:, 0] - cu)
            cover = len(np.unique(np.floor((ang + math.pi) / (2 * math.pi) * 36))) / 36
            if not Polygon(g.exterior).contains(Point(cu, cv)) or cover < MIN_RING:
                continue
            if best is None or g.area < best[0]:
                best = (g.area, float(np.median(np.hypot(ext[:, 0] - cu, ext[:, 1] - cv))), lvl)
    return None if best is None else best[1]


def shadow(to, axis):
    """The TO's outline along axis (union of its projected facets), as a shapely polygon."""
    u, v = [i for i in range(3) if i != axis]
    tri = to.vertices[to.faces][:, :, [u, v]]
    a = 0.5 * np.abs((tri[:, 1, 0] - tri[:, 0, 0]) * (tri[:, 2, 1] - tri[:, 0, 1])
                     - (tri[:, 2, 0] - tri[:, 0, 0]) * (tri[:, 1, 1] - tri[:, 0, 1]))
    polys = [Polygon(t) for t in tri[a > 1e-9]]
    return unary_union(polys).buffer(0)


def detect(original, to, pitch=0.5):
    """Hubs and web of a TO result, or None when the original has no hole the TO rings.
    Returns dict(axis, hubs=[dict(c, ri, ro, lo, hi)], web=(lo, hi), shadow=polygon)."""
    hs = holes(original)
    if not hs:
        return None
    axis = max({h['axis'] for h in hs}, key=lambda k: sum(h['axis'] == k for h in hs))
    hubs = []
    for h in hs:
        if h['axis'] != axis:
            continue
        ro = _ring_radius(to, h)
        if ro is not None and ro > h['r'] + 0.5:
            hubs.append(dict(c=h['c'], ri=h['r'], ro=ro, lo=h['lo'], hi=h['hi']))
    if not hubs:
        return None
    # web band: TO columns along the axis away from the hubs
    V = to.voxelized(pitch).fill()
    idx = np.argwhere(V.matrix)
    pts = trimesh.transformations.transform_points(idx, V.transform)
    u, v = [i for i in range(3) if i != axis]
    far = np.ones(len(pts), bool)
    for hb in hubs:
        far &= np.hypot(pts[:, u] - hb['c'][0], pts[:, v] - hb['c'][1]) > hb['ro'] + 2.0
    q = pts[far]
    if len(q) < 10:
        return None
    key = np.round(q[:, [u, v]] / pitch).astype(int)
    order = np.lexsort((key[:, 1], key[:, 0]))
    key, z = key[order], q[order, axis]
    brk = np.flatnonzero(np.any(np.diff(key, axis=0) != 0, axis=1)) + 1
    zmin = np.minimum.reduceat(z, np.r_[0, brk]); zmax = np.maximum.reduceat(z, np.r_[0, brk])
    lo = float(np.median(zmin) - pitch / 2); hi = float(np.median(zmax) + pitch / 2)
    return dict(axis=axis, hubs=hubs, web=(round(lo * 4) / 4, round(hi * 4) / 4), shadow=shadow(to, axis))
