"""Fit the material removed from an original section O by a TO section M as cuts made of
straight sides with tangent fillets (closed voids) and outline notches (open, closed outside O).

Section(O, M).fit() returns the cuts in their final form:
    ('poly', [(x, y, r), ...])                 closed polygon with a fillet radius per vertex
    ('arc', (cx, cy, r, a0, a1), [(x, y), ...]) arc notch, closed by a polyline outside O
Coordinates are the 2D coordinates of the section plane (mm).
"""
import numpy as np
from shapely.geometry import Polygon, Point, LinearRing
from shapely.ops import unary_union
import lf2
from lf2 import (denoise, fit_closed, fit_open, fit_open_arc, rounded, report, dist_to,
                 inside_depth)

COARSE = dict(tol=0.5, Lmin=4.0, accept=1.0)
FINE = dict(tol=0.15, Lmin=2.0, accept=0.35)
MAX_SIDES = 0           # most sides of a void (0: no limit)
CUT_MIN = 0.0           # % of the section area: smaller voids and notches of the TO are not cut

# Detail level 1-5: 1 = fewest dimensions (larger model error), 5 = smallest model error.
# tol: line fit tolerance (mm), Lmin: shortest side (mm), accept: largest profile error (mm),
# cut_min: smallest void or notch cut, % of the section area (voids under 1 mm^2 are never cut), min_flat: shortest straight between two fillets (mm),
# max_sides: most sides of a void (a notch: max_sides - 1 corners inside the part; 0 = no limit),
# radii: number of global corner radii (0 = radius per corner), digits: significant digits of dimensions
DETAIL = {
    1: dict(tol=2.0, Lmin=6.0, accept=3.0, cut_min=3.0, max_sides=3, min_flat=3.0, radii=1, digits=2),
    2: dict(tol=1.2, Lmin=5.0, accept=2.0, cut_min=1.5, max_sides=4, min_flat=2.0, radii=2, digits=2),
    3: dict(tol=0.8, Lmin=4.5, accept=1.5, cut_min=0.5, max_sides=5, min_flat=1.5, radii=2, digits=2),
    4: dict(tol=0.5, Lmin=4.0, accept=1.0, cut_min=0.0, max_sides=0, min_flat=1.0, radii=3, digits=3),
    5: dict(tol=0.25, Lmin=2.5, accept=0.5, cut_min=0.0, max_sides=0, min_flat=0.5, radii=0, digits=3),
}


def set_detail(level):
    """Fit settings of a detail level (1-5); returns the level's settings."""
    global CUT_MIN, MAX_SIDES
    d = DETAIL[level]
    COARSE.update(tol=d['tol'], Lmin=d['Lmin'], accept=d['accept'])
    FINE.update(tol=0.3 * d['tol'], Lmin=0.5 * d['Lmin'], accept=0.35 * d['accept'])
    CUT_MIN = d['cut_min']
    MAX_SIDES = d['max_sides']
    lf2.MIN_FLAT = d['min_flat']
    return d


def nseg(verts):
    return len(verts) + sum(v[2] > 1e-6 for v in verts)


def _line_x(p, q, r, t):
    """Intersection of line p-q with line r-t (None when nearly parallel)."""
    d1, d2 = q - p, t - r
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-9 * np.linalg.norm(d1) * np.linalg.norm(d2):
        return None
    s_ = ((r - p)[0] * d2[1] - (r - p)[1] * d2[0]) / den
    return p + s_ * d1


CUT_WEIGHT = 4.0
HOLE_WALL = 2.0         # mm: TO material this close to an original hole is kept by simplified cuts


def reduce_sides(V, n, ref, closed=True, cap=np.inf, keep=None):
    """Corner points V (k x 2) reduced to n corners (closed) or n inner corners (open: the two ends stay).
    Each step removes one corner, either by joining its neighbours with a straight side or by extending
    the two neighbouring sides to meet (one side removed), whichever keeps the polygon closest to ref
    (area cut from TO material x CUT_WEIGHT + area of TO void left as material; open profiles are closed by
    their end chord). Candidates that cut into keep are not taken. Stops early when the next step would put the outline more than cap from ref's
    (Hausdorff distance between the two boundaries)."""
    V = [np.asarray(v, float) for v in V]

    def cost(W):                                 # cutting TO material counts CUT_WEIGHT times leaving material
        P = Polygon(W).buffer(0)
        if not P.is_valid or P.area <= 0:
            return np.inf
        if keep is not None and P.intersection(keep).area > 1e-3:
            return np.inf                        # cuts into material that must stay (around the original's holes)
        return CUT_WEIGHT * P.difference(ref).area + ref.difference(P).area
    while (len(V) if closed else len(V) - 2) > n:
        m = len(V); best = None
        for k in (range(m) if closed else range(1, m - 1)):
            cands = [V[:k] + V[k + 1:]]                              # drop corner k
            j = (k + 1) % m                                          # remove side k-j: extend its neighbours
            if closed or j < m - 1:
                X = _line_x(V[k - 1], V[k], V[(j + 1) % m], V[j])
                if X is not None:
                    W = [X if i == k else v for i, v in enumerate(V) if i != j]
                    cands.append(W)
            for W in cands:
                c = cost(W)
                if best is None or c < best[0]:
                    best = (c, W)
        if best is None or not np.isfinite(best[0]):
            break
        if Polygon(best[1]).exterior.hausdorff_distance(ref.boundary) > cap:
            break
        V = best[1]
    return V


def _largest(M):
    if M.geom_type == 'Polygon':
        return M, []
    parts = sorted(M.geoms, key=lambda g: -g.area)
    return parts[0], parts[1:]


class Section:
    """Original section O (outline + its holes) and TO section M in one plane."""

    def __init__(self, O, M):
        self.O = O if O.geom_type == 'Polygon' else _largest(O)[0]
        self.Oout = Polygon(self.O.exterior)
        # the TO section is fitted unclipped (clipping puts corners on the outline that smoothing turns into notches)
        self.M, self.islands = _largest(M.buffer(0))
        orig = [Polygon(h) for h in self.O.interiors]
        self.voids = [h for h in self.M.interiors if Polygon(h).area > max(1.0, 0.01 * CUT_MIN * self.Oout.area) and not any(
            Polygon(h).symmetric_difference(o).area < 0.5 * o.area for o in orig)]
        self.R = self.Oout.buffer(4, join_style=2).exterior
        # TO material within HOLE_WALL of the original's holes: simplified cuts may not remove it
        holes = unary_union([Polygon(h) for h in self.O.interiors])
        self.keep = holes.buffer(HOLE_WALL).difference(holes).intersection(self.M) if not holes.is_empty else None

    # ---- closing a notch outside the original outline --------------------
    def closing(self, p_end, p_start, body):
        R = self.R; L = R.length
        se, ss = R.project(Point(p_end)), R.project(Point(p_start))
        rc = np.array(R.coords)[:-1]; pos = np.array([R.project(Point(q)) for q in rc])
        best = None
        for direc in (1, -1):
            span = (ss - se) % L if direc == 1 else (se - ss) % L
            rel = ((pos - se) % L) if direc == 1 else ((se - pos) % L)
            mid = [rc[i] for i in np.argsort(rel) if 1e-6 < rel[i] < span - 1e-6]
            path = [np.array(R.interpolate(se).coords[0])] + mid + [np.array(R.interpolate(ss).coords[0])]
            poly = Polygon(list(map(tuple, body)) + [tuple(q) for q in path]).buffer(0)
            a = poly.intersection(self.M).area     # the side that cuts the least TO material
            if best is None or a < best[0]:
                best = (a, path)
        return best[1]

    def finalize(self, raw):
        """Raw fit results -> cuts in final form."""
        cuts = []
        for v in raw['voids']:
            cuts.append(('poly', [tuple(map(float, x)) for x in v]))
        for n in raw['notches']:
            if n[0] == 'arc':
                _, c, r, a0, a1 = n
                body = c + r * np.c_[np.cos(np.linspace(a0, a1, 200)), np.sin(np.linspace(a0, a1, 200))]
                ring = self.Oout.buffer(1.5, join_style=2).exterior
                x1 = np.array(ring.interpolate(ring.project(Point(body[-1]))).coords[0])
                x0 = np.array(ring.interpolate(ring.project(Point(body[0]))).coords[0])
                path = [x1] + self.closing(x1, x0, np.vstack([body, x1])) + [x0]
                cuts.append(('arc', tuple(map(float, (c[0], c[1], r, a0, a1))), [tuple(map(float, q)) for q in path]))
            else:
                verts = [tuple(map(float, x)) for x in n]
                body = rounded(verts, closed=False)
                path = self.closing(np.array(verts[-1][:2]), np.array(verts[0][:2]), body)
                cuts.append(('poly', verts + [(float(q[0]), float(q[1]), 0.0) for q in path]))
        return cuts

    @staticmethod
    def cut_polygon(cut):
        if cut[0] == 'poly':
            return Polygon(rounded(cut[1])).buffer(0)
        (cx, cy, r, a0, a1), path = cut[1], cut[2]
        t = np.linspace(a0, a1, 200)
        return Polygon(list(zip(cx + r * np.cos(t), cy + r * np.sin(t))) + list(path)).buffer(0)

    def material(self, cuts):
        return Polygon(self.O.exterior, list(self.O.interiors)).difference(
            unary_union([self.cut_polygon(c) for c in cuts]))

    # ---- fitting -----------------------------------------------------------
    def _fit(self, fine=(), verbose=False):
        lf2.OUTLINE = self.Oout
        res = {'voids': [], 'notches': [], 'refs': []}; segs = 0; worst = 0; means = []
        log = []
        for iv, h in enumerate(self.voids):
            st = FINE if ('v', iv) in fine else COARSE
            lf2.AMIN = max(8.0, 30 * st['tol']); lf2.TOLM = 2 * st['tol']
            best = None
            for lm in (st['Lmin'], 0.75 * st['Lmin'], 0.5 * st['Lmin'], None):
                if lm is None:      # round or oval void: corners of a simplified outline, then fillets
                    raw, Qd = denoise(h)
                    P = Polygon(Qd).simplify(st['tol'] * 2)
                    pts = np.asarray(P.exterior.coords)[:-1]
                    if len(pts) < 3: continue
                    verts = lf2.fit_radii(Qd, [(x0, y0, 0.0) for x0, y0 in pts], True)
                else:
                    try: raw, verts = fit_closed(h, st['tol'], lm)
                    except (IndexError, ValueError): continue
                    if len(verts) < 3: continue
                C = rounded(verts)
                if len(C) < 4: continue
                m, x, ok = report(raw, C)
                bad = not ok or x > st['accept']
                cand = (not ok, bad, x if bad else nseg(verts), nseg(verts), m, verts, x, raw)
                if best is None or cand[:4] < best[:4]: best = cand
                if not bad and ok: break
            if MAX_SIDES and len(best[5]) > MAX_SIDES:     # fewer sides: corners reduced, then fillets
                raw = best[7]; ref = Polygon(h).buffer(0); _, Qd = denoise(h)
                V = reduce_sides([v[:2] for v in best[5]], MAX_SIDES, ref, cap=2 * st['accept'], keep=self.keep)
                verts = lf2.fit_radii(Qd, [(float(x), float(y), 0.0) for x, y in V], True)
                m, x, ok = report(raw, rounded(verts))
                best = (not ok, x > st['accept'], 0, nseg(verts), m, verts, x, raw)
            res['voids'].append(best[5]); res['refs'].append(('closed', best[7], None, self.Oout)); segs += best[3]; worst = max(worst, best[6]); means.append(best[4])
            log.append('void: %d segments, mean %.3f max %.3f%s' % (best[3], best[4], best[6], '' if not best[0] else ' (INVALID)'))
        raw, Q = denoise(self.M.exterior)
        ins = inside_depth(self.Oout, Q) > 0.05
        if ins.all():
            raise NotImplementedError('TO section does not reach the original outline: needs an outer cut')
        s0 = int(np.argmin(ins)); Q = np.roll(Q, -s0, 0); ins = np.roll(ins, -s0)
        k = 0; n = len(Q); jn = -1
        while k < n:
            if not ins[k]: k += 1; continue
            a = k
            while k < n and ins[k]: k += 1
            if inside_depth(self.Oout, Q[a:k]).max() < max(0.3, 2 * COARSE['tol']): continue   # shallower dips are surface noise
            if CUT_MIN:                          # small notch, below the detail level: not cut
                gone = self.Oout.difference(Polygon(Q).buffer(0))
                mid = Point(Q[(a + k) // 2]).buffer(0.5)
                if sum(g.area for g in getattr(gone, 'geoms', [gone]) if g.intersects(mid)) < 0.01 * CUT_MIN * self.Oout.area:
                    continue
            seg = Q[max(0, a - 8):min(n, k + 8)]; jn += 1
            st = FINE if ('n', jn) in fine else COARSE
            lf2.AMIN = max(8.0, 30 * st['tol']); lf2.TOLM = 2 * st['tol']
            rs = raw[dist_to(seg, raw, closed=False) < 0.5]; inner = rs[inside_depth(self.Oout, rs) > 0.3]
            if len(inner) < 3:
                inner = seg[inside_depth(self.Oout, seg) > 0.3]
            arc, adev = fit_open_arc(seg, self.Oout)
            best = None
            for lm in (st['Lmin'], 0.75 * st['Lmin'], 0.5 * st['Lmin']):
                try: verts = fit_open(seg, self.Oout, st['tol'], lm)
                except IndexError: continue
                C = rounded(verts, closed=False); m, x, ok = report(inner, C, closed=False)
                bad = not ok or x > st['accept']
                cand = (not ok, bad, x if bad else nseg(verts) - 2, nseg(verts) - 2, m, verts, x)
                if best is None or cand[:4] < best[:4]: best = cand
                if not bad and ok: break
            if arc is not None:
                ca, ra, b0, b1 = arc; tt = np.linspace(b0, b1, 400)
                adev = report(inner, ca + ra * np.c_[np.cos(tt), np.sin(tt)], closed=False)[1]
            if arc is not None and (adev <= st['accept'] or best is None or ((best[0] or best[1]) and adev < best[6])):
                c, r, a0, a1 = arc; t = np.linspace(a0, a1, 400)
                m, x, _ = report(inner, c + r * np.c_[np.cos(t), np.sin(t)], closed=False)
                res['notches'].append(('arc', c, r, a0, a1)); res['refs'].append(None); segs += 1; worst = max(worst, x); means.append(m)
                log.append('notch: single arc r %.1f, mean %.3f max %.3f' % (r, m, x))
            else:
                if False and MAX_SIDES and len(best[5]) - 2 > MAX_SIDES - 1:   # notches: no side limit (one notch can run
                                                                          # around several features, e.g. both arm tops)
                    ends = np.array([best[5][0][:2], best[5][-1][:2]])
                    ref = Polygon(np.vstack([seg, ends[1], ends[0]])).buffer(0)
                    V = reduce_sides([v[:2] for v in best[5]], MAX_SIDES - 1, ref, closed=False, cap=2 * st['accept'],
                                     keep=self.keep)
                    verts = lf2.fit_radii(seg, [(float(x), float(y), 0.0) for x, y in V], False)
                    verts = [(*verts[0][:2], 0.0)] + list(verts[1:-1]) + [(*verts[-1][:2], 0.0)]
                    m, x, ok = report(inner, rounded(verts, closed=False), closed=False)
                    best = (not ok, x > st['accept'], 0, nseg(verts) - 2, m, verts, x)
                res['notches'].append(best[5]); res['refs'].append(('open', inner, len(best[5]), self.Oout)); segs += best[3]; worst = max(worst, best[6]); means.append(best[4])
                log.append('notch: %d segments, mean %.3f max %.3f%s' % (best[3], best[4], best[6], '' if not best[0] else ' (INVALID)'))
        if verbose:
            print('\n'.join('    ' + s for s in log))
        return res, segs, (np.mean(means) if means else 0.0), worst

    def thin_strips(self, P, w=1.0, amin=0.05):
        import shapely
        P = shapely.set_precision(P.buffer(0), 1e-6)        # snapped: avoids GEOS non-noded intersections
        t = P.difference(shapely.set_precision(P.buffer(-w / 2, join_style=2).buffer(w / 2, join_style=2), 1e-6))
        return [g for g in (t.geoms if hasattr(t, 'geoms') else [t]) if g.area > amin]

    def fit(self, verbose=False):
        """Coarse fit; profiles that leave material strips < 1 mm (not already thin in the TO) are refit finely."""
        inTO = self.thin_strips(self.M)
        fine = set()
        for it in range(3):
            raw, segs, mean, worst = self._fit(fine)
            cuts = self.finalize(raw); self.refs = raw['refs']
            keys = [('v', i) for i in range(len(raw['voids']))] + [('n', j) for j in range(len(raw['notches']))]
            polys = [self.cut_polygon(c) for c in cuts]
            strips = [g for g in self.thin_strips(self.material(cuts)) if not any(g.distance(h) < 0.5 for h in inTO)]
            new = {keys[int(np.argmin([p.exterior.distance(g) for p in polys]))] for g in strips} if polys else set()
            if not (new - fine): break
            fine |= new
        if verbose:
            self._fit(fine, verbose=True)
        self.stats = dict(segments=segs, mean=mean, max=worst, fine=sorted(fine), thin_strips=len(strips),
                          islands=len(self.islands))
        return cuts
