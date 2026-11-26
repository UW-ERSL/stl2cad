"""Pattern detection and segmentation of a TO result against its design domain (the two STL files only).

A case folder holds the design domain STL and the TO STL (*_finalIsoBinary.stl); nothing else is read.

  1. Frame: principal directions of the domain's area-weighted face normals (sum of A n n^T). A degenerate pair
     (cylinder, disk) or triple (cube) is turned to the nearest global axes. The global axes are kept when the
     frame is within FRAME_SNAP deg of them or does not align ALIGN_GAIN more of the domain's face area with its
     axes (an irregular domain such as a bracket). Both meshes are expressed in the frame.
  2. Voxels of the domain D and the TO V (parity along z-rays), pitch = diagonal / RES (at most MAX_VOX voxels).
  3. Clean-up: TO voxels outside the domain removed (fraction reported); enclosed voids under VOID_F of the TO
     filled; islands under ISLAND_F of the largest removed.
  4. Symmetry about the domain's centroid: mirror planes (frame axes, diagonals of each pair, planes through an
     accepted rotation axis) and rotations C_n (n = 2..MAX_N) about each frame axis. Score = fraction of the
     voxels whose image lies within one voxel of the set, for the domain and for the TO. Accepted: domain >=
     SYM_DOMAIN and TO >= SYM_TO; candidate (reported, not imposed): TO >= SYM_CAND.
  5. Constraint tests along each frame axis, on the columns the TO meets: full (the TO fills the domain column),
     one interval (height field; centred on the domain column or anchored on one of its ends), two intervals.
  6. Retained regions: TO voxels where the TO fills the domain within a window of KEEP_W x diagonal
     (local fill >= KEEP_FILL): hubs, bosses, eye rings, solid blocks. Built later from the domain itself.
  7. Local shape at medial points (ridge of the distance transform), from the principal axes of the TO within
     SHAPE_R x the local radius: rod (s2/s1 < SHAPE_T), plate (s3/s2 < SHAPE_T), else blob; rod axis and plate
     normal kept. Labels spread to the TO voxels from the nearest medial point. Rods and plates are clustered by
     direction (DIR_COS); segments = connected parts per (type, direction); segments under MIN_SEG of the TO are
     merged into the neighbour they touch most; touching segments of one type whose directions differ by less
     than SMOOTH_DEG are merged (a curved shell or a curved rod).
  8. Per segment: volume share, thickness, plane fit (plates: flatness -> plate or shell), line fit (rods:
     straightness -> straight or curved, taper), contact with the domain's faces, neighbours.
  9. Part class: extrude (1) or height field (2a stepped / 2b continuous) from 5; otherwise from the volume
     shares of the segment kinds: rod frame (4), curved rods (5), plates (3), shells (6), solid mass (7), or a
     hybrid of those with share >= HYBRID_SHARE.

Outputs in CASE/cad/: analysis.json, analysis.txt, labels.npz (voxel segment labels, frame, grid).
Usage (from to/):  python analyze.py CASE_DIR [CASE_DIR ...] [--res 300]
                   python analyze.py --all ../casestudies
"""
import argparse, glob, json, math, os, sys, time
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

RES = 300               # voxels along the diagonal of the grid box
MAX_VOX = 1.2e7         # largest grid
FRAME_SNAP = 1.0        # deg
ALIGN_GAIN = 0.05       # a rotated frame is used only if it aligns this much more face area with its axes
DEGEN = 0.05           # eigenvalues of the normal tensor within this fraction of the largest are degenerate
VOID_F = 0.002          # enclosed voids smaller than this x the TO volume are filled
ISLAND_F = 0.01         # islands smaller than this x the largest component are removed
MAX_N = 12
SYM_DOMAIN = 0.97
SYM_TO = 0.93
SYM_CAND = 0.85
SYM_SAMPLES = 200000
FULL = 0.9              # a TO column is full when it fills this fraction of the domain column
EXTRUDE = (0.80, 0.90)  # (full, one interval) fractions of the columns for an extrude
HEIGHT = 0.95           # one-interval fraction for a height field
HEIGHT_POS = 0.9        # ... of which centred or anchored
STEP_LEVELS, STEP_F = 3, 0.9   # stepped (2a): this many depth levels hold this fraction of the columns
KEEP_W, KEEP_FILL, KEEP_MIN = 0.12, 0.95, 0.005
RIDGE_TOL = 0.5         # medial voxels: distance within this (voxels) of the largest in the 3x3x3 neighbourhood
SECTIONS = (0.25, 0.5, 0.75)   # extrude: the TO sections at these fractions of the depth ...
SECTION_IOU = 0.9       # ... overlap each other at least this much
SHAPE_R, SHAPE_T, SHAPE_MINR = 4.0, 0.45, 4.0
MEDIAL_SAMPLES = 6000
SMOOTH_K = 10           # medial samples voting on each sample's kind
DIR_COS = 0.94          # same direction (20 deg)
SMOOTH_DEG = 35.0       # touching segments of one type closer than this in direction are one segment
MIN_SEG = 0.003
FLAT = 0.2              # plate: RMS distance of the mid-surface from its plane < FLAT x thickness
STRAIGHT = 0.04         # rod: RMS distance of the axis points from their line < STRAIGHT x length
HYBRID_SHARE = 0.15
PURE_SHARE = 0.75

KIND_CLASS = {'rod': 4, 'curved rod': 5, 'plate': 3, 'shell': 6, 'mass': 7}


# ---------------------------------------------------------------- input

def load_stl(f):
    d = open(f, 'rb').read()
    n = int(np.frombuffer(d[80:84], np.uint32)[0]) if len(d) >= 84 else -1
    if len(d) == 84 + 50 * n:
        dt = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
        return np.frombuffer(d[84:], dtype=dt)['v'].astype(float)
    v = [l.split()[1:4] for l in d.decode(errors='ignore').splitlines() if l.strip().startswith('vertex')]
    return np.array(v, float).reshape(-1, 3, 3)


def find_pair(case):
    fs = sorted(glob.glob(os.path.join(case, '*.[sS][tT][lL]')))
    to = [f for f in fs if 'finaliso' in os.path.basename(f).lower()]
    dom = [f for f in fs if f not in to]
    if len(to) != 1 or len(dom) != 1:
        raise ValueError('%s: expected one domain STL and one *_finalIsoBinary.stl, found %s' % (case, fs))
    return dom[0], to[0]


# ---------------------------------------------------------------- 1. frame

def frame(tris):
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    a = np.linalg.norm(n, axis=1)
    ok = a > 0
    M = np.einsum('ij,ik->jk', n[ok], n[ok] / a[ok, None]) / 2
    w, E = np.linalg.eigh(M)
    w, E = w[::-1], E[:, ::-1]
    # degenerate groups -> nearest global axes inside the group's span
    groups, i = [], 0
    while i < 3:
        j = i + 1
        while j < 3 and abs(w[j] - w[i]) < DEGEN * w[0]:
            j += 1
        groups.append(list(range(i, j)))
        i = j
    R = np.zeros((3, 3))
    for g in groups:
        S = E[:, g]
        if len(g) == 1:
            R[:, g[0]] = S[:, 0]
            continue
        P = S @ S.T                               # projector on the span
        used = []
        for k in g:
            cand = [(np.linalg.norm(P[:, e]), e) for e in range(3) if e not in used]
            _, e = max(cand)
            v = P[:, e].copy()
            for q in [R[:, kk] for kk in g if kk < k]:
                v -= (v @ q) * q
            R[:, k] = v / np.linalg.norm(v)
            used.append(e)
    # order and sign the columns to match the global axes
    perm = sorted(range(3), key=lambda k: int(np.argmax(np.abs(R[:, k]))))
    if len({int(np.argmax(np.abs(R[:, k]))) for k in range(3)}) == 3:
        R = R[:, perm]
    for k in range(3):
        if R[np.argmax(np.abs(R[:, k])), k] < 0:
            R[:, k] *= -1
    if np.linalg.det(R) < 0:
        R[:, 2] *= -1
    dev = np.degrees(np.arccos(np.clip(np.abs(np.diag(R)), 0, 1)))
    # keep the global axes unless the frame aligns more of the domain's faces with its axes
    un = n[ok] / a[ok, None]

    def aligned(Q):
        return float(a[ok][(np.abs(un @ Q) > math.cos(math.radians(2))).any(1)].sum() / a[ok].sum())
    if dev.max() < FRAME_SNAP or aligned(np.eye(3)) >= aligned(R) - ALIGN_GAIN:
        return np.eye(3), dev * 0
    return R, dev


# ---------------------------------------------------------------- 2. voxels

def voxelize(tris, lo, shape, h):
    """Solid voxels (centres inside) by parity of the crossings of z-rays through the voxel centres."""
    nx, ny, nz = shape
    P = (tris - lo) / h - 0.5
    cols, zs = [], []
    for a, b, c in P:
        x0 = max(int(math.ceil(min(a[0], b[0], c[0]))), 0)
        x1 = min(int(math.floor(max(a[0], b[0], c[0]))), nx - 1)
        y0 = max(int(math.ceil(min(a[1], b[1], c[1]))), 0)
        y1 = min(int(math.floor(max(a[1], b[1], c[1]))), ny - 1)
        if x1 < x0 or y1 < y0:
            continue
        v0, v1 = b[:2] - a[:2], c[:2] - a[:2]
        den = v0[0] * v1[1] - v1[0] * v0[1]
        if abs(den) < 1e-12:
            continue
        gx, gy = np.meshgrid(np.arange(x0, x1 + 1), np.arange(y0, y1 + 1), indexing='ij')
        gx, gy = gx.ravel(), gy.ravel()
        px, py = gx - a[0], gy - a[1]
        u = (px * v1[1] - v1[0] * py) / den
        w = (v0[0] * py - px * v0[1]) / den
        m = (u >= -1e-9) & (w >= -1e-9) & (u + w <= 1 + 1e-9)
        if m.any():
            cols.append(gx[m] * ny + gy[m])
            zs.append(a[2] + u[m] * (b[2] - a[2]) + w[m] * (c[2] - a[2]))
    V = np.zeros(shape, bool)
    if not cols:
        return V, 0
    col = np.concatenate(cols)
    z = np.round(np.concatenate(zs), 6)
    o = np.lexsort((z, col))
    col, z = col[o], z[o]
    keep = np.r_[True, (col[1:] != col[:-1]) | (z[1:] != z[:-1])]
    col, z = col[keep], z[keep]
    start = np.r_[0, np.flatnonzero(col[1:] != col[:-1]) + 1]
    cnt = np.diff(np.r_[start, len(col)])
    rank = np.arange(len(col)) - np.repeat(start, cnt)
    good = np.repeat(cnt % 2 == 0, cnt)
    odd = int((cnt % 2 == 1).sum())
    s = good & (rank % 2 == 0)
    e = np.flatnonzero(s) + 1
    k0 = np.clip(np.ceil(z[s]), 0, nz).astype(int)
    k1 = np.clip(np.ceil(z[e]), 0, nz).astype(int)
    diff = np.zeros((nx * ny, nz + 1), np.int16)
    np.add.at(diff, (col[s], k0), 1)
    np.add.at(diff, (col[s], k1), -1)
    V = (np.cumsum(diff[:, :nz], axis=1) > 0).reshape(shape)
    return V, odd


# ---------------------------------------------------------------- 3. clean-up

def cleanup(V, D):
    out = float((V & ~D).sum() / max(V.sum(), 1))
    V = V & D
    lab, n = ndi.label(V, np.ones((3, 3, 3)))
    sizes = np.bincount(lab.ravel())[1:]
    keep = 1 + np.flatnonzero(sizes >= ISLAND_F * sizes.max())
    islands = int(n - len(keep))
    V = np.isin(lab, keep)
    E, ne = ndi.label(~np.pad(V, 1))
    es = np.bincount(E.ravel())
    outside = E[0, 0, 0]
    filled = 0
    for k in range(1, ne + 1):
        if k != outside and es[k] < VOID_F * V.sum():
            V |= (E[1:-1, 1:-1, 1:-1] == k)
            filled += 1
    voids = int(ne - 1 - filled)
    return V, dict(outside_domain=round(out, 4), islands_removed=islands, voids_filled=filled,
                   voids_kept=voids, components=int(len(keep)))


# ---------------------------------------------------------------- 4. symmetry

def rot(axis, ang):
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K


def mirror(nrm):
    n = nrm / np.linalg.norm(nrm)
    return np.eye(3) - 2 * np.outer(n, n)


class Sym:
    """Scores of a transform g about c: the fraction of a set's voxels whose image lies within one voxel."""

    def __init__(self, A, lo, h, c, rng):
        self.Ad = ndi.binary_dilation(A, np.ones((3, 3, 3)))
        idx = np.argwhere(A)
        if len(idx) > SYM_SAMPLES:
            idx = idx[rng.choice(len(idx), SYM_SAMPLES, replace=False)]
        self.P = lo + (idx + 0.5) * h
        self.lo, self.h, self.c = lo, h, c

    def score(self, g):
        Q = self.c + (self.P - self.c) @ g.T
        I = np.rint((Q - self.lo) / self.h - 0.5).astype(int)
        ok = np.all((I >= 0) & (I < np.array(self.Ad.shape)), axis=1)
        hit = np.zeros(len(I), bool)
        hit[ok] = self.Ad[I[ok, 0], I[ok, 1], I[ok, 2]]
        return float(hit.mean())


def symmetry(D, V, lo, h, rng):
    c = lo + (np.argwhere(D).mean(0) + 0.5) * h
    sd, sv = Sym(D, lo, h, c, rng), Sym(V, lo, h, c, rng)
    res = dict(centre=c.tolist(), mirrors=[], rotations=[], candidates=[])

    def judge(kind, vec, g, extra):
        d, v = sd.score(g), sv.score(g)
        rec = dict(kind=kind, vector=np.round(vec, 4).tolist(), domain=round(d, 3), to=round(v, 3), **extra)
        if d >= SYM_DOMAIN and v >= SYM_TO:
            return rec, 'accepted'
        if d >= SYM_DOMAIN and v >= SYM_CAND:
            return rec, 'candidate'
        return rec, None

    # rotations about the frame axes: the largest n accepted, else the largest candidate
    for ax in range(3):
        a = np.eye(3)[ax]
        best, cand = None, None
        for n in range(2, MAX_N + 1):
            rec, st = judge('rotation', a, rot(a, 2 * math.pi / n), dict(n=n))
            if st == 'accepted':
                ok = all(sv.score(rot(a, 2 * math.pi * k / n)) >= SYM_TO for k in range(2, n // 2 + 1))
                if ok:
                    best = rec
            elif st == 'candidate':
                cand = rec
        if best:
            res['rotations'].append(best)
        elif cand:
            res['candidates'].append(cand)
    # mirror planes: frame axes, diagonals, planes through accepted rotation axes
    normals = [np.eye(3)[k] for k in range(3)]
    for i in range(3):
        for j in range(i + 1, 3):
            for s in (1, -1):
                normals.append((np.eye(3)[i] + s * np.eye(3)[j]) / math.sqrt(2))
    for r in res['rotations']:
        a = np.array(r['vector'])
        b = np.eye(3)[int(np.argmin(np.abs(a)))]
        b = b - (b @ a) * a
        b /= np.linalg.norm(b)
        cc = np.cross(a, b)
        best = max((sv.score(mirror(math.cos(t) * b + math.sin(t) * cc)), t)
                   for t in np.radians(np.arange(0, 180.0 / r['n'], 1.0)))
        for k in range(r['n']):
            t = best[1] + k * math.pi / r['n']
            normals.append(math.cos(t) * b + math.sin(t) * cc)
    seen = []
    for nrm in normals:
        if any(abs(nrm @ q) > 0.999 for q in seen):
            continue
        seen.append(nrm)
        rec, st = judge('mirror', nrm, mirror(nrm), {})
        if st == 'accepted':
            res['mirrors'].append(rec)
        elif st == 'candidate':
            res['candidates'].append(rec)
    return res


# ---------------------------------------------------------------- 5. constraint tests

def runs(W):
    d = np.diff(np.pad(W.astype(np.int8), ((0, 0), (1, 1))), axis=1)
    return (d == 1).sum(1), np.argmax(d == 1, axis=1), np.argmax(d == -1, axis=1)


def constraint_tests(V, D):
    out = []
    for ax in range(3):
        W = np.moveaxis(V, ax, -1).reshape(-1, V.shape[ax])
        Wd = np.moveaxis(D, ax, -1).reshape(-1, V.shape[ax])
        m = W.any(1)
        W, Wd = W[m], Wd[m]
        nr, s0, s1 = runs(W)
        dnr, d0, d1 = runs(Wd)
        dl = Wd.sum(1)
        full = (W & Wd).sum(1) >= FULL * dl
        one = nr == 1
        mid = np.abs((s0 + s1) - (d0 + d1)) <= 3            # in half voxels
        anch = (np.abs(s0 - d0) <= 1) | (np.abs(s1 - d1) <= 1)
        depth = (s1 - s0)[one]
        lv = np.bincount(depth).astype(float) if len(depth) else np.zeros(1)
        held, rest = 0.0, lv.copy()
        for _ in range(STEP_LEVELS):                        # depth levels +-1 voxel, greedily
            cv = np.convolve(rest, np.ones(3), 'same')
            i = int(np.argmax(cv))
            held += cv[i]
            rest[max(i - 1, 0):i + 2] = 0
        stepped = float(held / max(lv.sum(), 1))
        # sections across the axis at fractions of the domain's extent along it
        ext = np.flatnonzero(D.any(axis=tuple(k for k in range(3) if k != ax)))
        secs = [np.take(V, int(ext[0] + f * (ext[-1] - ext[0])), ax) for f in SECTIONS]
        sec_iou = min(float((a & b).sum() / max((a | b).sum(), 1)) for i, a in enumerate(secs) for b in secs[i + 1:])
        # two layers: two TO intervals in a column where the domain has one, each on one end of it
        two_l = (nr == 2) & (dnr == 1)
        if two_l.any():
            last = W.shape[1] - 1 - np.argmax(W[two_l][:, ::-1], axis=1)
            two_l_f = float(((np.abs(s0[two_l] - d0[two_l]) <= 1) & (np.abs(last - (d1[two_l] - 1)) <= 1)).sum()
                            / len(W))
        else:
            two_l_f = 0.0
        out.append(dict(axis='xyz'[ax], full=round(float(full.mean()), 3), one=round(float(one.mean()), 3),
                        two=round(float((nr == 2).mean()), 3), two_layers=round(two_l_f, 3),
                        section_iou=round(sec_iou, 3),
                        centred=round(float(mid[one].mean()) if one.any() else 0.0, 3),
                        anchored=round(float(anch[one].mean()) if one.any() else 0.0, 3),
                        step_levels=round(stepped, 3)))
    return out


def route(tests):
    ex = [t for t in tests if t['full'] >= EXTRUDE[0] and t['one'] >= EXTRUDE[1] and t['section_iou'] >= SECTION_IOU]
    if ex:
        t = max(ex, key=lambda t: t['full'])
        return dict(route='extrude', axis=t['axis'], cls='1')
    hf = [t for t in tests if t['one'] >= HEIGHT and max(t['centred'], t['anchored']) >= HEIGHT_POS]
    if hf:
        t = max(hf, key=lambda t: t['one'])
        kind = 'centred' if t['centred'] >= t['anchored'] else 'anchored'
        return dict(route='height field', axis=t['axis'], kind=kind,
                    cls='2a' if t['step_levels'] >= STEP_F else '2b')
    two = [t for t in tests if t['two_layers'] >= 0.6]
    if two:
        return dict(route='two layers', axis=max(two, key=lambda t: t['two_layers'])['axis'], cls=None)
    return dict(route='3D', axis=None, cls=None)


# ---------------------------------------------------------------- 6. retained regions

def retained(V, D, h, diag):
    w = max(5, int(round(KEEP_W * diag / h)) | 1)
    f = ndi.uniform_filter(V.astype(np.float32), w) / np.maximum(ndi.uniform_filter(D.astype(np.float32), w), 1e-6)
    K = V & (f >= KEEP_FILL)
    lab, n = ndi.label(K)
    sizes = np.bincount(lab.ravel())[1:]
    keep = 1 + np.flatnonzero(sizes >= KEEP_MIN * V.sum())
    return np.isin(lab, keep) if len(keep) else np.zeros_like(V), w


# ---------------------------------------------------------------- 7. local shape

def cluster_dirs(U):
    """Greedy clustering of unit directions (sign ignored): cluster id per row, mean direction per cluster."""
    lab = -np.ones(len(U), int)
    means = []
    C = np.abs(U @ U.T) if len(U) <= 6000 else None
    while (lab < 0).any():
        free = np.flatnonzero(lab < 0)
        if C is not None:
            cnt = (C[np.ix_(free, free)] >= DIR_COS).sum(1)
        else:
            cnt = (np.abs(U[free] @ U[free].T) >= DIR_COS).sum(1)
        seed = free[int(np.argmax(cnt))]
        mem = free[np.abs(U[free] @ U[seed]) >= DIR_COS]
        S = np.einsum('ij,ik->jk', U[mem], U[mem])
        m = np.linalg.eigh(S)[1][:, -1]
        mem = free[np.abs(U[free] @ m) >= DIR_COS]
        if len(mem) == 0:
            mem = np.array([seed])
        lab[mem] = len(means)
        means.append(m)
    return lab, np.array(means)


def local_shape(V, K, rng):
    DT = ndi.distance_transform_edt(np.pad(V, 1))[1:-1, 1:-1, 1:-1]
    ridge = V & (DT >= ndi.maximum_filter(DT, 3) - RIDGE_TOL) & (DT >= 1.0) & ~K
    R = np.argwhere(ridge)
    if len(R) == 0:
        return None
    if len(R) > MEDIAL_SAMPLES:
        R = R[rng.choice(len(R), MEDIAL_SAMPLES, replace=False)]
    rad = np.maximum(SHAPE_R * DT[tuple(R.T)], SHAPE_MINR)
    kind, vec, ratio = [], [], []
    for p, r in zip(R, rad):
        # the TO inside the ball that is connected to the medial point (not a parallel wall nearby)
        ri = int(math.ceil(r))
        a = np.maximum(p - ri, 0)
        b = np.minimum(p + ri + 1, V.shape)
        sub = V[a[0]:b[0], a[1]:b[1], a[2]:b[2]]
        g = np.indices(sub.shape).reshape(3, -1).T + a - p
        ball = (np.sum(g * g, 1) <= r * r).reshape(sub.shape) & sub
        lab, _ = ndi.label(ball, np.ones((3, 3, 3)))
        c = tuple(p - a)
        q = np.argwhere(lab == lab[c]) + a - p if lab[c] else np.zeros((1, 3))
        e, Ev = np.linalg.eigh(np.cov(q.T) if len(q) > 3 else np.eye(3))
        s = np.sqrt(np.maximum(e[::-1], 1e-12))
        Ev = Ev[:, ::-1]
        if s[1] / s[0] < SHAPE_T:
            kind.append('rod')
        elif s[2] / s[1] < SHAPE_T:
            kind.append('plate')
        else:
            kind.append('blob')
        vec.append((Ev[:, 0], Ev[:, 2]))           # rod axis, plate normal
        ratio.append((s[1] / s[0], s[2] / s[1]))
    kind = np.array(kind)
    ax_, nm_ = np.array([v[0] for v in vec]), np.array([v[1] for v in vec])
    # majority of the SMOOTH_K nearest medial samples (isolated votes at edges and junctions)
    nb = cKDTree(R).query(R, min(SMOOTH_K, len(R)))[1]
    codes = np.searchsorted(['blob', 'plate', 'rod'], kind)
    votes = np.stack([(codes[nb] == c).sum(1) for c in range(3)], 1)
    kind = np.array(['blob', 'plate', 'rod'])[np.argmax(votes, 1)]
    vec = np.where((kind == 'rod')[:, None], ax_, np.where((kind == 'plate')[:, None], nm_, 0.0))
    return dict(DT=DT, R=R, kind=kind, vec=vec, ratio=np.array(ratio))


# ---------------------------------------------------------------- 7. segments

def touching(L):
    """Shared-face counts between labels (> 0) of a label volume: dict {(a, b): count}, a < b."""
    pairs = {}
    for ax in range(3):
        a = np.moveaxis(L, ax, 0)
        x, y = a[:-1].ravel(), a[1:].ravel()
        m = (x != y) & (x > 0) & (y > 0)
        if not m.any():
            continue
        lo_, hi_ = np.minimum(x[m], y[m]), np.maximum(x[m], y[m])
        key, cnt = np.unique(lo_.astype(np.int64) * 100000 + hi_, return_counts=True)
        for k, c in zip(key, cnt):
            pairs[(int(k // 100000), int(k % 100000))] = pairs.get((int(k // 100000), int(k % 100000)), 0) + int(c)
    return pairs


def segments(V, K, sh):
    """Label volume: 1 = retained, 2.. = segments; per segment (type, direction)."""
    R, kind, vec = sh['R'], sh['kind'], sh['vec']
    code = np.zeros(len(R), int)                   # type and direction cluster -> one code
    dirs = {}
    nxt = 1
    for t in ('rod', 'plate'):
        m = np.flatnonzero(kind == t)
        if len(m):
            lab, means = cluster_dirs(vec[m])
            for k in range(len(means)):
                code[m[lab == k]] = nxt
                dirs[nxt] = (t, means[k])
                nxt += 1
    mb = np.flatnonzero(kind == 'blob')
    if len(mb):
        code[mb] = nxt
        dirs[nxt] = ('blob', np.zeros(3))
    W = V & ~K
    idx = np.argwhere(W)
    near = cKDTree(R).query(idx)[1]
    C = np.zeros(V.shape, np.int32)
    C[tuple(idx.T)] = code[near]
    # connected parts per code
    L = np.zeros(V.shape, np.int32)
    L[K] = 1
    info = {1: dict(type='retained', dir=np.zeros(3))}
    nl = 2
    for cde, (t, d) in dirs.items():
        lab, n = ndi.label(C == cde, np.ones((3, 3, 3)))
        if n:
            L[lab > 0] = lab[lab > 0] + nl - 1
            for k in range(n):
                info[nl + k] = dict(type=t, dir=d)
            nl += n
    tot = V.sum()
    # merge small segments into the neighbour they touch most
    for _ in range(20):
        vol = np.bincount(L.ravel(), minlength=nl)
        small = [k for k in info if k != 1 and 0 < vol[k] < MIN_SEG * tot]
        if not small:
            break
        pairs = touching(L)
        remap = np.arange(nl)
        for k in small:
            nb = [(c, b if a == k else a) for (a, b), c in pairs.items() if k in (a, b)]
            nb = [x for x in nb if remap[x[1]] == x[1] and x[1] not in small]
            if nb:
                remap[k] = max(nb)[1]
        if (remap == np.arange(nl)).all():
            break
        L = remap[L]
        for k in list(info):
            if remap[k] != k:
                info.pop(k)
    # merge touching segments of one type with close directions (curved shell, curved rod)
    pairs = touching(L)
    parent = {k: k for k in info}

    def find(k):
        while parent[k] != k:
            k = parent[k]
        return k
    cos = math.cos(math.radians(SMOOTH_DEG))
    for (a, b), c in sorted(pairs.items(), key=lambda x: -x[1]):
        if a in info and b in info and info[a]['type'] == info[b]['type'] in ('rod', 'plate'):
            if abs(info[a]['dir'] @ info[b]['dir']) >= cos:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[rb] = ra
    remap = np.arange(L.max() + 1)
    for k in info:
        remap[k] = find(k)
    L = remap[L]
    info = {k: v for k, v in info.items() if find(k) == k}
    # renumber
    ks = sorted(info)
    rn = np.zeros(L.max() + 1, np.int32)
    for i, k in enumerate(ks):
        rn[k] = i + 1
    L = rn[L]
    return L, [info[k] for k in ks]


# ---------------------------------------------------------------- 8. segment properties

def describe(L, info, V, D, sh, lo, h):
    DT, R = sh['DT'], sh['R']
    dd = ndi.distance_transform_edt(np.pad(D, 1))[1:-1, 1:-1, 1:-1]
    surf = V & ~ndi.binary_erosion(np.pad(V, 1))[1:-1, 1:-1, 1:-1]
    on_face = surf & (dd <= 1.5)
    tot = V.sum()
    vol = np.bincount(L.ravel(), minlength=len(info) + 1)
    sv = np.bincount(L[surf], minlength=len(info) + 1)
    fv = np.bincount(L[on_face], minlength=len(info) + 1)
    ridge = V & (DT >= ndi.maximum_filter(DT, 3) - RIDGE_TOL) & (DT >= 1.0)
    rl = L[ridge]
    rp = np.argwhere(ridge)
    rr = DT[ridge]
    pairs = touching(L)
    out = []
    for k, inf in enumerate(info, start=1):
        if vol[k] == 0:
            continue
        m = rl == k
        P, r = rp[m] * h, rr[m] * h
        seg = dict(id=k, type=inf['type'], share=round(float(vol[k] / tot), 3),
                   face_contact=round(float(fv[k] / max(sv[k], 1)), 2),
                   neighbours=sorted([b if a == k else a for (a, b) in pairs if k in (a, b)]))
        if len(P) >= 3:
            seg['thickness'] = round(float(2 * np.median(r)), 4)
            c = P.mean(0)
            e, Ev = np.linalg.eigh(np.cov((P - c).T))
            e, Ev = e[::-1], Ev[:, ::-1]
            if inf['type'] == 'plate':
                n = Ev[:, 2]
                rms = float(np.sqrt(np.mean(((P - c) @ n) ** 2)))
                seg['normal'] = np.round(n * np.sign(n[np.argmax(np.abs(n))]), 3).tolist()
                seg['flatness'] = round(rms / max(seg['thickness'], h), 3)
                seg['kind'] = 'plate' if seg['flatness'] < FLAT else 'shell'
            elif inf['type'] == 'rod':
                a = Ev[:, 0]
                t = (P - c) @ a
                length = float(t.max() - t.min())
                rms = float(np.sqrt(np.mean(np.sum(((P - c) - np.outer(t, a)) ** 2, 1))))
                seg['axis'] = np.round(a * np.sign(a[np.argmax(np.abs(a))]), 3).tolist()
                seg['length'] = round(length, 4)
                seg['straightness'] = round(rms / max(length, h), 3)
                seg['kind'] = 'rod' if seg['straightness'] < STRAIGHT else 'curved rod'
                if len(t) > 5 and np.ptp(t) > 0:
                    k1, k0 = np.polyfit(t, r, 1)
                    ends = sorted([abs(k0 + k1 * t.min()), abs(k0 + k1 * t.max())])
                    seg['taper'] = round(ends[1] / max(ends[0], h / 2), 2)
            elif inf['type'] == 'blob':
                seg['kind'] = 'mass'               # set to 'node' / 'junction' below
            else:
                seg['kind'] = 'retained'
        else:
            seg['kind'] = {'retained': 'retained', 'blob': 'mass'}.get(inf['type'], inf['type'])
        out.append(seg)
    # a blob between members is a joint: a node of rods or a junction of plates and shells
    by = {s['id']: s for s in out}
    for s in out:
        if s['kind'] != 'mass':
            continue
        nb = [by[k] for k in s['neighbours'] if k in by]
        rods = [n for n in nb if n['kind'] in ('rod', 'curved rod')]
        walls = [n for n in nb if n['kind'] in ('plate', 'shell')]
        for grp, name in ((rods, 'node'), (walls, 'junction')):
            if len(grp) >= 2 and s['share'] <= sum(n['share'] for n in grp):
                s['kind'] = name
                s['joins'] = max(grp, key=lambda n: n['share'])['kind']
                break
    return out


def part_class(rt, segs):
    if rt['cls']:
        return rt['cls'], 'route %s along %s' % (rt['route'], rt['axis'])
    share = {}
    for s in segs:
        k = s['kind']
        if k in ('retained',):
            continue
        if k in ('node', 'junction'):
            k = s['joins']              # joints count with the members they join
        share[k] = share.get(k, 0) + s['share']
    tot = sum(share.values())
    if tot == 0:
        return '7', 'retained regions only'
    share = {k: v / tot for k, v in share.items()}
    top = max(share, key=share.get)
    if share[top] >= PURE_SHARE:
        return str(KIND_CLASS[top]), '%s %.0f %%' % (top, 100 * share[top])
    parts = sorted([k for k, v in share.items() if v >= HYBRID_SHARE], key=lambda k: -share[k])
    return '8/9 (%s)' % ' + '.join(str(KIND_CLASS[k]) for k in parts), \
        ', '.join('%s %.0f %%' % (k, 100 * share[k]) for k in parts)


# ---------------------------------------------------------------- driver

def analyze(case, res=RES, log=print):
    t0 = time.time()
    rng = np.random.default_rng(0)
    fd, ft = find_pair(case)
    dom, to = load_stl(fd), load_stl(ft)
    Rf, dev = frame(dom)
    dom, to = dom @ Rf, to @ Rf
    pd, pt = dom.reshape(-1, 3), to.reshape(-1, 3)
    lo = np.minimum(pd.min(0), pt.min(0))
    hi = np.maximum(pd.max(0), pt.max(0))
    diag = float(np.linalg.norm(pd.max(0) - pd.min(0)))
    h = max(np.linalg.norm(hi - lo) / res, (np.prod(hi - lo) / MAX_VOX) ** (1 / 3))
    lo = lo - h
    shape = tuple(int(x) for x in np.ceil((hi - lo) / h) + 1)
    D, od = voxelize(dom, lo, shape, h)
    V, ot = voxelize(to, lo, shape, h)
    V, clean = cleanup(V, D)
    clean.update(odd_columns_domain=od, odd_columns_to=ot)
    sym = symmetry(D, V, lo, h, rng)
    tests = constraint_tests(V, D)
    rt = route(tests)
    K, kw = retained(V, D, h, diag)
    sh = local_shape(V, K, rng)
    if sh is not None:
        L, info = segments(V, K, sh)
        segs = describe(L, info, V, D, sh, lo, h)
    else:
        L, segs = K.astype(np.int32), [dict(id=1, type='retained', kind='retained', share=1.0)]
    cls, why = part_class(rt, segs)
    rep = dict(case=os.path.abspath(case), domain_stl=os.path.basename(fd), to_stl=os.path.basename(ft),
               frame=dict(R=np.round(Rf, 6).tolist(), deviation_deg=np.round(dev, 2).tolist()),
               grid=dict(lo=lo.tolist(), pitch=h, shape=list(shape)),
               domain=dict(size=np.round(pd.max(0) - pd.min(0), 6).tolist(),
                           fill_of_box=round(float(D.sum() / np.prod(np.ceil((pd.max(0) - pd.min(0)) / h))), 3)),
               volume_fraction=round(float(V.sum() / D.sum()), 3), cleanup=clean, symmetry=sym,
               tests=tests, route=rt, retained=dict(share=round(float(K.sum() / V.sum()), 3), window=kw),
               segments=segs, part_class=cls, reason=why, seconds=round(time.time() - t0, 1))
    od_ = os.path.join(case, 'cad')
    os.makedirs(od_, exist_ok=True)
    with open(os.path.join(od_, 'analysis.json'), 'w') as f:
        json.dump(rep, f, indent=1)
    np.savez_compressed(os.path.join(od_, 'labels.npz'), labels=L.astype(np.int16), lo=lo, pitch=h, R=Rf)
    txt = report(rep)
    with open(os.path.join(od_, 'analysis.txt'), 'w') as f:
        f.write(txt)
    log(txt)
    return rep


def report(r):
    s = r['symmetry']
    L = ['%s  (domain %s, TO %s)' % (r['case'], r['domain_stl'], r['to_stl']),
         'class %s: %s' % (r['part_class'], r['reason']),
         'frame deviation from global axes (deg): %s' % r['frame']['deviation_deg'],
         'domain %s, fill of its box %.3f; TO volume fraction %.3f; pitch %.4g; %.1f s' % (
             r['domain']['size'], r['domain']['fill_of_box'], r['volume_fraction'], r['grid']['pitch'],
             r['seconds']),
         'clean-up: %s' % r['cleanup'],
         'symmetry: rotations %s; mirrors %s' % (
             ['C%d about %s (%.2f/%.2f)' % (x['n'], x['vector'], x['domain'], x['to']) for x in s['rotations']],
             ['%s (%.2f/%.2f)' % (x['vector'], x['domain'], x['to']) for x in s['mirrors']]),
         '  candidates (not imposed): %s' % [
             ('C%d ' % x['n'] if x['kind'] == 'rotation' else 'mirror ') + '%s (%.2f/%.2f)' % (
                 x['vector'], x['domain'], x['to']) for x in s['candidates']],
         'tests (full / one / two / two layers / section overlap / centred / anchored / step levels):']
    for t in r['tests']:
        L.append('  %s: %.2f / %.2f / %.2f / %.2f / %.2f / %.2f / %.2f / %.2f' % (
            t['axis'], t['full'], t['one'], t['two'], t['two_layers'], t['section_iou'], t['centred'],
            t['anchored'], t['step_levels']))
    L.append('route: %s' % r['route'])
    L.append('retained regions: %.1f %% of the TO' % (100 * r['retained']['share']))
    L.append('segments (share, kind, thickness, direction, flatness/straightness, face contact, neighbours):')
    for g in sorted(r['segments'], key=lambda g: -g['share']):
        d = g.get('normal') or g.get('axis') or ''
        f = g.get('flatness', g.get('straightness', ''))
        L.append('  %2d %5.1f %%  %-10s t=%-8s %-22s %-6s face %.2f  nb %s' % (
            g['id'], 100 * g['share'], g['kind'], g.get('thickness', ''), d, f, g['face_contact'],
            g['neighbours']))
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cases', nargs='*')
    ap.add_argument('--all', help='analyse every folder below this one that holds a *_finalIsoBinary.stl')
    ap.add_argument('--res', type=int, default=RES)
    a = ap.parse_args()
    cases = list(a.cases)
    if a.all:
        for f in sorted(glob.glob(os.path.join(a.all, '**', '*'), recursive=True)):
            if 'finaliso' in os.path.basename(f).lower() and f.lower().endswith('.stl'):
                cases.append(os.path.dirname(f))
    for c in cases:
        try:
            analyze(c, a.res)
        except Exception as e:                    # keep going through a batch
            print('%s: FAILED %s: %s' % (c, type(e).__name__, e), file=sys.stderr)


if __name__ == '__main__':
    main()
