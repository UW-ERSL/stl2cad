"""Rotational symmetry of a plate's cuts: detect it, replace each orbit of cuts by turned copies of one
seed, and record the pattern for the SolidWorks script (seed Cut-Extrudes + a Circular Pattern).

A plate is n-fold symmetric when the original section turned by 360/n deg about its centre is unchanged
(area mismatch < ORIG_TOL) and the TO section nearly so (< TO_TOL). Each cut is matched to the cuts its
turned copies land on (overlap > 50 %); a full orbit of n cuts is replaced by the member whose turned
copies best overlap the others, turned n times. Corner radii are kept (the seed's).
"""
import math
import numpy as np
from shapely import affinity
from shapely.geometry import Polygon

ORIG_TOL = 0.002     # area mismatch of the original section after the turn
TO_TOL = 0.10        # area mismatch of the TO section after the turn
MAX_N = 12


def _turn_xy(q, c, ang):
    ca, sa = math.cos(ang), math.sin(ang)
    x, y = q[0] - c[0], q[1] - c[1]
    return c[0] + ca * x - sa * y, c[1] + sa * x + ca * y


def turn_cut(cut, c, ang):
    """A cut ('poly' or 'arc') turned by ang (rad) about c."""
    if cut[0] == 'poly':
        return ('poly', [(*_turn_xy(v, c, ang), v[2]) for v in cut[1]])
    (cx, cy, r, a0, a1), path = cut[1], cut[2]
    x, y = _turn_xy((cx, cy), c, ang)
    return ('arc', (x, y, r, a0 + ang, a1 + ang), [_turn_xy(q, c, ang) for q in path])


def order(sec):
    """Largest n (2..MAX_N) for which the plate section is n-fold symmetric about the original's centre,
    with the TO mismatch; (1, None, None) if none."""
    O = Polygon(sec.O.exterior, list(sec.O.interiors))
    c = np.array(sec.Oout.centroid.coords[0])
    best = (1, None, None)
    for n in range(2, MAX_N + 1):
        ang = 360.0 / n
        if affinity.rotate(O, ang, origin=tuple(c)).symmetric_difference(O).area > ORIG_TOL * O.area:
            continue
        mis = affinity.rotate(sec.M, ang, origin=tuple(c)).symmetric_difference(sec.M).area / sec.M.area
        if mis < TO_TOL:
            best = (n, c, mis)
    return best


def symmetrize(p):
    """Detect the plate's rotational symmetry; when found, replace each full orbit of cuts by turned copies
    of its seed. Sets p['pattern'] = dict(n, centre, mismatch, orbits=[[cut index, ...], ...]) with the seed
    first in each orbit, and reorders p['cuts'] (and p['refs']) so that all seeds come first, then the
    copies orbit by orbit, then the cuts outside any orbit. Returns n (1: no symmetry)."""
    sec = p['sec']
    n, c, mis = order(sec)
    if n == 1:
        return 1
    cuts = p['cuts']
    polys = [sec.cut_polygon(k) for k in cuts]
    step = 2 * math.pi / n
    used, orbits = set(), []
    for i in range(len(cuts)):
        if i in used:
            continue
        members = [i]
        for k in range(1, n):
            t = affinity.rotate(polys[i], math.degrees(k * step), origin=tuple(c))
            hit = [j for j in range(len(cuts)) if j not in used and j not in members
                   and polys[j].intersection(t).area > 0.5 * min(polys[j].area, t.area)]
            if not hit:
                break
            members.append(hit[0])
        if len(members) == n:
            orbits.append(members)
            used.update(members)

    def fit_of(i):     # how well cut i's turned copies cover the orbit's other members
        return sum(affinity.rotate(polys[i], math.degrees(k * step), origin=tuple(c)).symmetric_difference(
            polys[m]).area for k, m in enumerate(orb) if k)
    new_cuts, new_refs, seeds, copies = [], [], [], []
    refs = p.get('refs') or [None] * len(cuts)
    for orb in orbits:
        # the member whose turned copies best match the others, members re-indexed from it
        ranked = []
        for s0 in range(n):
            rot = orb[s0:] + orb[:s0]
            err = sum(affinity.rotate(polys[rot[0]], math.degrees(k * step), origin=tuple(c)).symmetric_difference(
                polys[rot[k]]).area for k in range(1, n))
            ranked.append((err, rot))
        err, rot = min(ranked)
        seeds.append((cuts[rot[0]], refs[rot[0]]))
        copies.append([turn_cut(cuts[rot[0]], c, k * step) for k in range(1, n)])
    rest = [(cuts[i], refs[i]) for i in range(len(cuts)) if i not in used]
    for s in seeds:
        new_cuts.append(s[0]); new_refs.append(s[1])
    ns = len(seeds)
    for cl in copies:
        new_cuts += cl; new_refs += [None] * len(cl)
    for r in rest:
        new_cuts.append(r[0]); new_refs.append(r[1])
    p['cuts'], p['refs'] = new_cuts, new_refs
    p['pattern'] = dict(n=n, centre=tuple(map(float, c)), mismatch=float(mis), seeds=ns,
                        copies=ns * (n - 1))
    return n
