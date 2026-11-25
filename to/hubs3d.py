"""Hubs on any axis: rings of a TO result around the original part's holes, for 3D frames (frame3d.py).

detect(original, to) returns a list of hubs, dict(axis, c=(u, v), ri, ro, lo, hi), where (u, v) are the two
other coordinates in increasing axis order (axis 0: (y, z); 1: (x, z); 2: (x, y)):
  1. holes: cylindrical faces of the original that close all round (360 deg, possibly in several faces) with
     material outside them; a partial cylinder (a fillet of a concave corner) is not a hole;
  2. coaxial holes are separate hubs when the original has no material between them along the axis (the two
     eyes of a clevis), so each hub covers one stretch of material;
  3. ring: at levels along the hub, the TO section's pieces whose outline encloses the hole; the outermost one
     (a hollow boss gives two nested rings) at the level where its area is smallest (least joined to members);
     ro = the median distance of its outline from the axis. A hole the TO does not enclose is not a hub.
hubweb.detect is unchanged: it handles hubs on one axis with a flat web (swsweep.py, swbars.py).
"""
import math
import numpy as np
from shapely.geometry import Point, Polygon
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Cylinder
import cadquery as cq

import plates

MIN_RING = 0.6       # the ring's outline covers at least this fraction of the directions round the hole
FULL = 2 * math.pi - 1e-6


def other(axis):
    return [i for i in range(3) if i != axis]


def holes(original, gap=0.5):
    """Full cylindrical holes of the original, one entry per stretch of material along the axis:
    list of dict(axis, c, r, lo, hi)."""
    faces = []
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
        u, v = other(k)
        bb = f.BoundingBox()
        lo, hi = (bb.xmin, bb.ymin, bb.zmin)[k], (bb.xmax, bb.ymax, bb.zmax)[k]
        span = ad.LastUParameter() - ad.FirstUParameter()
        faces.append(dict(axis=k, c=(float(P[u]), float(P[v])), r=float(r), lo=float(lo), hi=float(hi), span=span))
    # faces of one cylinder (same axis line and radius, overlapping along the axis) together close all round
    groups = []
    for f in faces:
        for g in groups:
            h = g[0]
            if (h['axis'] == f['axis'] and math.dist(h['c'], f['c']) < 1e-3 and abs(h['r'] - f['r']) < 1e-3
                    and f['lo'] < h['hi'] - 1e-6 and h['lo'] < f['hi'] - 1e-6):
                g.append(f); break
        else:
            groups.append([f])
    out = []
    for g in groups:
        if sum(f['span'] for f in g) < FULL:
            continue                                   # a fillet or a partial cylinder
        h = dict(axis=g[0]['axis'], c=g[0]['c'], r=g[0]['r'], lo=min(f['lo'] for f in g), hi=max(f['hi'] for f in g))
        k = h['axis']; u, v = other(k)
        # material outside the cylinder, half way along it
        q = [0.0, 0.0, 0.0]; q[k] = 0.5 * (h['lo'] + h['hi']); q[u] = h['c'][0] + h['r'] + 0.2; q[v] = h['c'][1]
        if not original.isInside(cq.Vector(*q), 1e-6):
            continue
        out.append(h)
    # coaxial holes touching along the axis (a hole through two faces split by OCC) are one hole
    out.sort(key=lambda h: (h['axis'], round(h['c'][0], 3), round(h['c'][1], 3), h['lo']))
    merged = []
    for h in out:
        p = merged[-1] if merged else None
        if (p is not None and p['axis'] == h['axis'] and math.dist(p['c'], h['c']) < 1e-3
                and h['lo'] <= p['hi'] + gap):
            p['hi'] = max(p['hi'], h['hi']); p['r'] = min(p['r'], h['r'])
        else:
            merged.append(dict(h))
    return merged


def ring(to, h, nlev=9):
    """(ro, level) of the TO ring round hole h: the outermost section piece enclosing the hole, at the level where
    its area is smallest; None if the TO does not enclose the hole."""
    k, (cu, cv) = h['axis'], h['c']
    best = None
    for t in np.linspace(0.1, 0.9, nlev):
        lvl = h['lo'] + t * (h['hi'] - h['lo'])
        sec = plates.section(to, k, lvl)
        enc = []
        for g in getattr(sec, 'geoms', [sec]):
            if g.is_empty:
                continue
            ext = Polygon(g.exterior)
            if not ext.contains(Point(cu, cv)):
                continue
            pts = np.asarray(g.exterior.coords)
            ang = np.arctan2(pts[:, 1] - cv, pts[:, 0] - cu)
            cover = len(np.unique(np.floor((ang + math.pi) / (2 * math.pi) * 36))) / 36
            if cover >= MIN_RING:
                enc.append((ext.area, pts))
        if not enc:
            continue
        area, pts = max(enc, key=lambda e: e[0])           # outermost
        if best is None or area < best[0]:
            best = (area, float(np.median(np.hypot(pts[:, 0] - cu, pts[:, 1] - cv))), float(lvl))
    return None if best is None else (best[1], best[2])


def detect(original, to):
    """Hubs (see the module docstring), each also with 'n' (its axis as a unit vector)."""
    out = []
    for h in holes(original):
        r = ring(to, h)
        if r is None or r[0] <= h['r'] + 0.5:
            continue
        n = [0.0, 0.0, 0.0]; n[h['axis']] = 1.0
        out.append(dict(axis=h['axis'], c=h['c'], ri=h['r'], ro=r[0], lo=h['lo'], hi=h['hi'], n=tuple(n)))
    return out


def point(h, s):
    """Model point on hub h's axis at level s along it."""
    k = h['axis']; u, v = other(k)
    p = [0.0, 0.0, 0.0]; p[k] = s; p[u], p[v] = h['c']
    return np.array(p)


def inside(h, P, margin=0.0):
    """Mask of points P (n x 3) inside hub h's cylinder (radius ro + margin, its stretch along the axis)."""
    k = h['axis']; u, v = other(k)
    return ((np.hypot(P[:, u] - h['c'][0], P[:, v] - h['c'][1]) <= h['ro'] + margin)
            & (P[:, k] >= h['lo'] - margin) & (P[:, k] <= h['hi'] + margin))
