"""Corner radii from a few global values.

Every filleted corner of every cut of a part takes one of at most three radius values (levels).
A corner can take a level only if its fillet fits: the tangent points stay on its two sides
(shared with the neighbouring fillets) and an acute tip is shortened by at most TIP_LOSS. Each
corner takes the feasible level whose local fit to the TO outline is best. The fewest levels are
used for which every profile stays within tolerance: max error <= max(ACCEPT, its error with
fitted radii + SLACK).
"""
import itertools
import numpy as np
import lf2
from lf2 import rounded, report, dist_to, dense
from scipy.spatial import cKDTree

GRID = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
R_MIN = 1.0      # smallest radius, except at corners too tight for it (acute tips, short sides)
R_ARC = 10.0     # fitted radii above this are curved sides of the shape, not fillets: kept as fitted
TIP_LOSS = 0.5
ACCEPT = 1.0
SLACK = 0.25


def _geom(V, k):
    n = len(V); a, b, c = V[k - 1], V[k], V[(k + 1) % n]
    u = (b - a) / np.linalg.norm(b - a); w = (c - b) / np.linalg.norm(c - b)
    return a, b, c, u, w, np.arccos(np.clip(u @ w, -1, 1))


def _cap(verts, k):
    """Largest radius corner k can take: tangent points on its sides (half of a side when the
    neighbour is filleted too),."""
    V = np.array([v[:2] for v in verts], float); R = [v[2] for v in verts]; n = len(V)
    a, b, c, u, w, th = _geom(V, k)
    fa = 0.5 if R[k - 1] > 0 else 0.95
    fb = 0.5 if R[(k + 1) % n] > 0 else 0.95
    tmax = min(fa * np.linalg.norm(b - a), fb * np.linalg.norm(c - b))
    sa = np.sin((np.pi - th) / 2)
    return tmax / max(np.tan(th / 2), 1e-9), TIP_LOSS * sa / max(1e-9, 1 - sa)


def _corner_error(verts, k, r, ref):
    """Max deviation (two-way) of the profile with corner k at radius r, other corners as fitted."""
    kind, pts, nb, lf2.OUTLINE = ref
    V = [tuple(v) for v in verts]; V[k] = (V[k][0], V[k][1], r)
    if kind == 'closed':
        m, x, ok = report(pts, rounded(V), closed=True)
    else:
        m, x, ok = report(pts, rounded(V[:nb], closed=False), closed=False)
    return x if ok else np.inf


def quantize(plates, max_levels=3, verbose=True, accept=None):
    """Replace the fitted corner radii of all cuts of all plates by at most max_levels values.
    plates: list of dicts with 'cuts' and 'refs' (TO points per cut, None for arc notches)."""
    corners = []                                   # (plate, cut, vertex, cap, fitted r)
    for ip, p in enumerate(plates):
        for ic, (cut, ref) in enumerate(zip(p['cuts'], p['refs'])):
            if cut[0] != 'poly' or ref is None:
                continue
            nb = ref[2] or len(cut[1])
            for k, v in enumerate(cut[1][:nb]):
                if 0 < v[2] <= R_ARC and (ref[0] == 'closed' or 0 < k < nb - 1):
                    corners.append((ip, ic, k, *_cap(cut[1], k), v[2]))
    if not corners:
        return dict(levels=[], counts={})
    E = np.full((len(corners), len(GRID)), np.inf)
    for i, (ip, ic, k, cap, rtip, r0) in enumerate(corners):
        verts, ref = plates[ip]['cuts'][ic][1], plates[ip]['refs'][ic]
        for j, R in enumerate(GRID):
            if R <= cap + 1e-9 and (R >= R_MIN or min(cap, rtip) < R_MIN):
                E[i, j] = _corner_error(verts, k, R, ref)
        if not np.isfinite(E[i]).any():            # not even the smallest level fits: smallest anyway
            E[i, 0] = _corner_error(verts, k, GRID[0], ref)

    def profile_errors(assign):
        out = {}
        for (ip, ic) in {(c[0], c[1]) for c in corners}:
            verts = [list(v) for v in plates[ip]['cuts'][ic][1]]
            for i, c in enumerate(corners):
                if c[0] == ip and c[1] == ic:
                    verts[c[2]][2] = assign[i]
            kind, pts, nb, lf2.OUTLINE = plates[ip]['refs'][ic]
            if kind == 'closed':
                m, x, ok = report(pts, rounded([tuple(v) for v in verts]), closed=True)
            else:
                m, x, ok = report(pts, rounded([tuple(v) for v in verts[:nb]], closed=False), closed=False)
            out[(ip, ic)] = x if ok else np.inf
        return out

    base = profile_errors([c[5] for c in corners])
    lim = {}
    for i, c in enumerate(corners):
        key = (c[0], c[1]); lim.setdefault(key, max(ACCEPT if accept is None else accept, base[key] + SLACK))
    L = np.array([lim[(c[0], c[1])] for c in corners])
    best = None
    for nl in range(1, max_levels + 1):
        cands = []
        for combo in itertools.combinations(range(len(GRID)), nl):
            sub = E[:, combo]
            over = np.maximum(sub.min(1) - L, 0)          # corners that miss the tolerance on their own
            cands.append(((over > 0).sum(), over.sum(), np.minimum(sub, 50).min(1).sum(), combo))
        cands.sort(key=lambda c: c[:3])
        for *_, combo in cands[:5]:                      # the best few, checked with all corners changed together
            pick = [GRID[combo[int(np.argmin(E[i, list(combo)]))]] for i in range(len(corners))]
            errs = profile_errors(pick)
            worst = max(errs[key] - lim[key] for key in errs)
            if best is None or worst < best[0] - 1e-9 and (best[0] > 0):
                best = (worst, pick, nl, errs)
            if worst <= 0:
                break
        if verbose > 1:
            bad = sorted(((best[3][k] - lim[k], k) for k in best[3]), reverse=True)[:3]
            print('  %d level(s) %s: worst excess %.2f mm; %s' % (nl, sorted(set(best[1])), best[0], ', '.join('plate %d cut %d +%.2f' % (k[0], k[1], e) for e, k in bad if e > 0)))
        if best[0] <= 0:
            break
    worst, pick, nl, errs = best
    for i, c in enumerate(corners):
        cut = plates[c[0]]['cuts'][c[1]]
        v = cut[1][c[2]]; cut[1][c[2]] = (v[0], v[1], float(pick[i]))
    levels = sorted(set(pick))
    counts = {R: pick.count(R) for R in levels}
    if verbose:
        print('corner radii: %d corners from %d global value(s): %s; profile max error %.2f mm (fitted radii %.2f mm)%s'
              % (len(corners), len(levels), ', '.join('%g mm x %d' % (R, counts[R]) for R in levels),
                 max(errs.values()), max(base.values()), '' if worst <= 0 else '; tolerance not met at %d values' % nl))
    return dict(levels=levels, counts=counts, max_error=max(errs.values()), base_error=max(base.values()))
