"""Decompose a topology-optimized solid into flat plates, each a 2.5D feature.

A voxel belongs to a plate normal to axis a when the solid is thin along a there: the run of
filled voxels through it along a is at most TMAX and clearly shorter than along the other two axes.
Voxels thin along the same axis form connected groups; each large group is a plate with a level
range [lo, hi] along its axis (median start and end of the runs through its voxels).
A 2.5D part (a plate with through-cuts) is the special case of one plate spanning the full thickness.
"""
import numpy as np
from scipy import ndimage
import trimesh
from shapely.geometry import Polygon
from shapely.ops import unary_union

TMAX = 8.0      # thickest plate (mm)
RATIO = 0.5     # thin run <= RATIO x the next shortest run
MIN_FRAC = 0.03 # smallest plate, as a fraction of the TO volume


def runs(F, ax):
    """Per voxel: start index, end index (inclusive) and length of its run of filled voxels along ax."""
    G = np.moveaxis(F, ax, -1); n = G.shape[-1]
    start = np.zeros(G.shape, int); end = np.zeros(G.shape, int)
    s = np.zeros(G.shape[:-1], int); prev = np.zeros(G.shape[:-1], bool)
    for i in range(n):
        s = np.where(G[..., i], np.where(prev, s, i), 0); start[..., i] = s; prev = G[..., i]
    e = np.zeros(G.shape[:-1], int); nxt = np.zeros(G.shape[:-1], bool)
    for i in range(n - 1, -1, -1):
        e = np.where(G[..., i], np.where(nxt, e, i), 0); end[..., i] = e; nxt = G[..., i]
    st, en = np.moveaxis(start, -1, ax), np.moveaxis(end, -1, ax)
    return st, en, (en - st + 1) * F


PRISM_TOL = 0.10   # a prism: sections at 25, 50, 75 % of the extent differ by < 10 % of their area


def prism_axis(mesh, bounds):
    """Axis along which the TO is a prism through the whole part (2.5D), or None."""
    best = None
    for a in range(3):
        lv = [bounds[0][a] + f * (bounds[1][a] - bounds[0][a]) for f in (0.25, 0.5, 0.75)]
        S = [section(mesh, a, l) for l in lv]
        if S[1].is_empty or S[1].area == 0: continue
        d = max(S[i].symmetric_difference(S[1]).area for i in (0, 2)) / S[1].area
        if d < PRISM_TOL and (best is None or d < best[1]): best = (a, d)
    return best


def detect(mesh, pitch=0.5, bounds=None, verbose=True):
    if bounds is not None:
        pa = prism_axis(mesh, bounds)
        if pa is not None:
            a = pa[0]
            if verbose:
                print('prism along %s (sections at 25/50/75 %% differ by %.1f %%): one plate through the part'
                      % ('xyz'[a], 100 * pa[1]))
            return [dict(axis=a, lo=float(bounds[0][a]), hi=float(bounds[1][a]), frac=1.0)]
    vg = mesh.voxelized(pitch).fill()
    F = vg.matrix.astype(bool)
    origin = vg.transform[:3, 3]                  # centre of voxel (0, 0, 0)
    R = [runs(F, a) for a in range(3)]
    L = np.stack([r[2] for r in R], -1).astype(float) * pitch
    L[~F] = np.inf
    order = np.sort(L, -1)
    thin = np.argmin(L, -1)
    plate_vox = F & (order[..., 0] <= TMAX) & (order[..., 0] <= RATIO * order[..., 1])
    total = F.sum(); plates = []
    for a in range(3):
        lab, n = ndimage.label(plate_vox & (thin == a))
        for k in range(1, n + 1):
            sel = lab == k
            if sel.sum() < MIN_FRAC * total:
                continue
            st, en = R[a][0][sel], R[a][1][sel]
            lo = origin[a] + (np.median(st) - 0.5) * pitch
            hi = origin[a] + (np.median(en) + 0.5) * pitch
            if bounds is not None:
                lo, hi = max(lo, bounds[0][a]), min(hi, bounds[1][a])
            plates.append(dict(axis=a, lo=float(lo), hi=float(hi), frac=float(sel.sum() / total)))
    # merge plates on the same axis whose ranges overlap
    plates.sort(key=lambda p: (p['axis'], p['lo']))
    merged = []
    for p in plates:
        q = merged[-1] if merged else None
        if q and q['axis'] == p['axis'] and p['lo'] < q['hi'] - 0.5 * pitch:
            q['lo'], q['hi'] = min(q['lo'], p['lo']), max(q['hi'], p['hi']); q['frac'] += p['frac']
        else:
            merged.append(dict(p))
    if verbose:
        print('voxels %d (pitch %.2f mm); in plates %.0f %%' % (total, pitch, 100 * plate_vox.sum() / total))
        for p in merged:
            print('  plate normal to %s: %.2f to %.2f mm (%.1f mm thick), %.0f %% of the TO voxels'
                  % ('xyz'[p['axis']], p['lo'], p['hi'], p['hi'] - p['lo'], 100 * p['frac']))
    return merged


def frame(axis):
    """Section-plane axes (u, v) for a plate normal to axis: the other two axes in increasing order."""
    return [i for i in range(3) if i != axis]


def section(mesh, axis, level):
    """Section of a mesh at a level along axis, as a shapely geometry in (u, v)."""
    u, v = frame(axis)
    o = [0.0, 0.0, 0.0]; o[axis] = level; nrm = [0.0, 0.0, 0.0]; nrm[axis] = 1.0
    s = mesh.section(plane_origin=o, plane_normal=nrm)
    if s is None:
        return Polygon()
    V = s.vertices[:, [u, v]]
    # polygons from the planar section, built directly in (u, v)
    s2 = trimesh.path.Path2D(entities=[e.copy() for e in s.entities], vertices=V)
    return unary_union([Polygon(p.exterior, list(p.interiors)) for p in s2.polygons_full]).buffer(0)
