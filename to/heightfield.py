"""Height-field TO results (class 2b): nested symmetric plates.

A TO result made with a draw (casting) constraint and mid-plane symmetry is a height field: every column along
the draw axis holds one solid interval centred on the mid-plane, so the part is its thickness map t(u, v). The
section at a distance h from the mid-plane is {t >= 2h}: the sections shrink as h grows, and the area A(h)
falls from the footprint (h = 0) to the full-depth region (h = H, half the domain's depth).

The part is modelled as K nested plates, all centred on the mid-plane. Plate j spans |z - mid| <= b_j
(0 < b_1 < ... < b_K = H) and its section is the TO section at a height r_j in (b_{j-1}, b_j]. Between
b_{j-1} and b_j the model's section is then that one section (the thinner plates are wider), so the model is
a piecewise-constant approximation of the sections in h. The breaks b_j are chosen by dynamic programming to
minimise the squared error of the area, sum over h of (A(h) - mean A over its interval)^2, and r_j is the
height in the interval whose area is closest to the interval's mean area (the volume of each layer is kept).

levels(to, bounds) -> (axis, mid, H, [(b_j, r_j, area)]); plates(...) -> plate dicts for run_plates.py
(sec_at: the level of the TO section to fit; order: the thickest plate is built first, so that the cuts of a
thick plate never remove a thinner plate's material: each thinner plate is added after it). A section in several
separate regions (the full-depth chords of a cantilever) gives one plate per region (island: the region fitted);
such parts are built additively in SolidWorks (one Extruded Boss of each plate's fitted region, per_cut=False),
since the cuts of one region's plate would remove the other regions of the same layer.
"""
import numpy as np
import plates as _plates

N_LEVELS = {1: 2, 2: 3, 3: 3, 4: 4, 5: 5}   # plates per detail level
SAMPLES = 40                               # heights sampled between the mid-plane and the face
MIN_AREA = 0.01                            # a layer whose section is below this x the footprint is dropped


def draw_axis(to, bounds, pitch=1.0):
    """The axis along which most TO columns hold one solid interval, and that fraction."""
    F = to.voxelized(pitch).fill().matrix.astype(bool)
    best = None
    for a in range(3):
        W = np.moveaxis(F, a, -1).reshape(-1, F.shape[a])
        W = W[W.any(1)]
        d = np.diff(np.pad(W.astype(np.int8), ((0, 0), (1, 1))), axis=1)
        one = float(((d == 1).sum(1) == 1).mean())
        if best is None or one > best[1]:
            best = (a, one)
    return best


def _partition(A, K):
    """Breaks of the best K-piece constant fit of A (least squares): end indices of the pieces."""
    n = len(A)
    c1, c2 = np.r_[0, np.cumsum(A)], np.r_[0, np.cumsum(A * A)]

    def cost(i, j):                       # samples i..j-1
        s, s2, m = c1[j] - c1[i], c2[j] - c2[i], j - i
        return s2 - s * s / m
    E = np.full((K + 1, n + 1), np.inf)
    B = np.zeros((K + 1, n + 1), int)
    E[0, 0] = 0
    for k in range(1, K + 1):
        for j in range(k, n + 1):
            for i in range(k - 1, j):
                e = E[k - 1, i] + cost(i, j)
                if e < E[k, j]:
                    E[k, j], B[k, j] = e, i
    ends, j = [], n
    for k in range(K, 0, -1):
        ends.append(j)
        j = B[k, j]
    return ends[::-1]


def levels(to, bounds, K, axis=None):
    a = draw_axis(to, bounds)[0] if axis is None else axis
    lo, hi = float(bounds[0][a]), float(bounds[1][a])
    mid, H = 0.5 * (lo + hi), 0.5 * (hi - lo)
    h = (np.arange(SAMPLES) + 0.5) * H / SAMPLES
    A = np.array([_plates.section(to, a, mid + x).area for x in h])
    foot = A.max()
    ends = _partition(A, min(K, SAMPLES))
    out, i0 = [], 0
    for j in ends:
        seg = slice(i0, j)
        b = H if j == SAMPLES else 0.5 * (h[j - 1] + h[j])
        m = A[seg].mean()
        r = h[seg][int(np.argmin(np.abs(A[seg] - m)))]
        if m >= MIN_AREA * foot:
            out.append((float(b), float(r), float(m)))
        i0 = j
    return a, mid, H, out


def plates(to, bounds, detail, axis=None, log=print):
    a, mid, H, L = levels(to, bounds, N_LEVELS[detail], axis)
    log('height field along %s: mid-plane %.2f, half-depth %.2f mm; %d nested plate(s)' % ('xyz'[a], mid, H, len(L)))
    foot = _plates.section(to, a, mid).area
    P = []
    for n, (b, r, m) in enumerate(L[::-1]):            # thickest first
        S = _plates.section(to, a, mid + r)
        parts = [g for g in getattr(S, 'geoms', [S]) if g.area >= MIN_AREA * foot]
        log('  plate %.2f mm thick, section at %.2f mm from the mid-plane (area %.0f mm^2)%s'
            % (2 * b, r, m, ', %d separate regions: one plate each' % len(parts) if len(parts) > 1 else ''))
        for g in parts:
            P.append(dict(axis=a, lo=mid - b, hi=mid + b, kind='Plate', frac=0.0, sec_at=mid + r, order=n,
                          island=g if len(parts) > 1 else None))
    return P
