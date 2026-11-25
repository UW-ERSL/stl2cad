"""Members of a truss-like TO web: curved centre lines between joints, with a section that varies along them.

members(shadow, hubs, to, axis, detail) works in the plane normal to axis:
  1. medial axis of the web region (skeleton.py's raster, pruning and graph);
  2. joints = branch points (pixel clusters within MERGE merged) and the points where chains meet a hub
     ring; a member = one pixel chain between two joints, not split at bends;
  3. path = polyline (Douglas-Peucker, TOL) with a tangent arc at each bend (radius fitted to the chain),
     so the path is lines and arcs, tangent-continuous, ready to sweep;
  4. section = rectangle w x d normal to the path, centred on the mid-plane; w and d vary linearly from the
     member's start to its end (least squares over stations along the chain): w from the distance
     transform, scaled so that w x length matches the member's share of the shadow (area-equivalent);
     d from rays through the TO.
Returns dict(joints=[(u, v)], hub_of={joint: hub}, members=[dict(i, j, path, length, w0, w1, d0, d1)], mid).
"""
import math
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
from shapely.geometry import Point, LineString
from skimage.morphology import skeletonize

import skeleton as S

TOL = {1: 2.0, 2: 1.5, 3: 1.0, 4: 0.7, 5: 0.5}      # path tolerance (mm) by detail level
MERGE = 3.0           # branch pixels within this form one joint (mm)
HUB_REACH = 3.0       # a chain end within this of a hub's ring + margin ends on the ring (mm)
MIN_MEMBER = 2.0      # chains shorter than this between two joints are dropped (mm)
SHORT_MEMBER = 5.0    # a member shorter than this merges its two joints (mm)
HUB_BAND = 4.0        # a hub member whose chain lies mostly within this of the hub's face goes straight to it (mm)


def _medial(shadow, hubs, pix):
    region = shadow
    for h in hubs:
        region = region.difference(Point(*h['c']).buffer(h['ro'] + S.HUB_MARGIN))
    region = region.buffer(-S.OPEN).buffer(S.OPEN).buffer(0)
    mask, (x0, y0) = S._raster(region, pix)
    edt = ndimage.distance_transform_edt(mask) * pix
    sk = skeletonize(mask)
    xy = lambda p: np.array([x0 + (p[0] + 0.5) * pix, y0 + (p[1] + 0.5) * pix])
    disks = [Point(*h['c']).buffer(h['ro']) for h in hubs]
    for _ in range(5):                                   # prune short free spurs (not at a hub)
        nodes, chains, deg = S._graph(sk)
        cut = False
        for ch in chains:
            free = [deg[ch[0]] == 1, deg[ch[-1]] == 1]
            near_hub = any(min(Point(*xy(ch[k])).distance(d) for d in disks) < 2.0 for k in (0, -1)) if disks else False
            if any(free) and not all(free) and len(ch) * pix < S.SPUR and not near_hub:
                for p in (ch[:-1] if free[0] else ch[1:]):
                    sk[p] = False
                cut = True
        if not cut:
            break
    return sk, xy, edt, (x0, y0), region


def _aim(V, aim):
    """Ends on a hub (aim: {0 or -1: (centre, ro)}) on its ring, on the line from the neighbouring vertex to the
    hub's centre."""
    V = [np.asarray(v, float) for v in V]
    for end, (c, ro) in (aim or {}).items():
        prev = V[1] if end == 0 else V[-2]
        V[end] = c + (prev - c) / np.linalg.norm(prev - c) * ro
    return V


def _fillet_path(V, chain, tol, rmin=0.0, aim=None):
    """Polyline V (ends fixed) with a tangent arc at each interior vertex: the radius (largest that fits
    between the neighbouring bends) chosen on a grid to bring the path closest to the chain points, and
    not below rmin (a swept section of half-width h folds over itself on an arc of radius < h). A bend
    with no room for rmin is dropped (its neighbours joined) and the path made again.
    Returns edges [('line', p, q) | ('arc', s, m, e)] head to tail."""
    V = _aim(V, aim)
    n = len(V)
    for k in range(1, n - 1):
        a, b = V[k] - V[k - 1], V[k + 1] - V[k]
        La, Lb = np.linalg.norm(a), np.linalg.norm(b)
        turn = math.acos(float(np.clip(a @ b / (La * Lb), -1, 1)))
        tmax = min(0.5 * La if k > 1 else 0.9 * La, 0.5 * Lb if k < n - 2 else 0.9 * Lb)
        if turn > 1e-3 and tmax < rmin * math.tan(turn / 2):
            return _fillet_path(V[:k] + V[k + 1:], chain, tol, rmin, aim)   # ends aimed again
    if n == 2:
        return [('line', V[0], V[1])]
    L = [np.linalg.norm(V[k + 1] - V[k]) for k in range(n - 1)]
    tree = chain
    out_t = []                                           # per interior vertex: (t, r, turn) tangent length
    for k in range(1, n - 1):
        a, b = V[k] - V[k - 1], V[k + 1] - V[k]
        a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
        turn = math.acos(float(np.clip(a @ b, -1, 1)))
        if turn < 1e-3:
            out_t.append((0.0, 0.0, 0.0)); continue
        share = lambda m: L[m] * (1.0 if m in (0, n - 2) and False else 0.5)
        tmax = min(0.5 * L[k - 1] if k > 1 else 0.9 * L[0], 0.5 * L[k] if k < n - 2 else 0.9 * L[-1])
        # local chain points around the vertex
        near = chain[np.linalg.norm(chain - V[k], axis=1) < max(tmax, 1e-6) * 1.2]
        best = None
        tmin = rmin * math.tan(turn / 2)
        for f in np.linspace(0.1, 1.0, 10):
            t = tmin + f * (tmax - tmin); r = t / math.tan(turn / 2)
            p, q = V[k] - a * t, V[k] + b * t
            bis = (b - a) / np.linalg.norm(b - a)
            c = V[k] + bis * r / math.cos(turn / 2)
            dev = float(np.max(np.abs(np.linalg.norm(near - c, axis=1) - r))) if len(near) else 0.0
            # points beyond the tangent points are on the lines, not the arc: judge only the arc's sector
            if best is None or dev < best[0]:
                best = (dev, t, r)
        out_t.append((best[1], best[2], turn))
    edges, start = [], V[0]
    for k in range(1, n - 1):
        t, r, turn = out_t[k - 1]
        a, b = V[k] - V[k - 1], V[k + 1] - V[k]
        a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
        if t <= 0:
            edges.append(('line', start, V[k])); start = V[k]; continue
        p, q = V[k] - a * t, V[k] + b * t
        bis = (b - a) / np.linalg.norm(b - a)
        c = V[k] + bis * r / math.cos(turn / 2)
        m = c + (V[k] - c) / np.linalg.norm(V[k] - c) * r
        if np.linalg.norm(p - start) > 1e-6:
            edges.append(('line', start, p))
        edges.append(('arc', p, m, q)); start = q
    edges.append(('line', start, V[-1]))
    return [tuple([e[0]] + [tuple(map(float, x)) for x in e[1:]]) for e in edges]


def path_points(path, h=0.25):
    """Points along a path (mm), spacing about h, and their arc length."""
    P = []
    for e in path:
        if e[0] == 'line':
            p, q = np.array(e[1]), np.array(e[2])
            k = max(2, int(np.linalg.norm(q - p) / h) + 1)
            P += list(p + np.linspace(0, 1, k)[:, None] * (q - p))
        else:
            s, m, t = (np.array(x) for x in e[1:])
            c = np.array(S_circle(s, m, t)); r = np.linalg.norm(s - c)
            a0, am, a1 = (math.atan2(*(v - c)[::-1]) for v in (s, m, t))
            sw = (a1 - a0) % (2 * math.pi)
            if not ((am - a0) % (2 * math.pi) < sw):
                sw -= 2 * math.pi
            k = max(3, int(abs(sw) * r / h) + 1)
            ang = a0 + np.linspace(0, sw, k)
            P += list(c + r * np.column_stack([np.cos(ang), np.sin(ang)]))
    P = np.array(P)
    keep = np.r_[True, np.linalg.norm(np.diff(P, axis=0), axis=1) > 1e-9]
    P = P[keep]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    return P, s


def S_circle(p, q, r):
    ax, ay = p; bx, by = q; cx, cy = r
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    return ux, uy


def members(shadow, hubs, to=None, axis=2, detail=3, pix=S.PIX):
    tol = TOL[detail]
    sk, xy, edt, (x0, y0), region = _medial(shadow, hubs, pix)
    nodes, chains, deg = S._graph(sk)
    # joints: branch-pixel clusters (degree >= 3) merged within MERGE; free ends are joints too
    branch = [xy(p) for p in nodes if deg[p] >= 3]
    J, owner = [], {}
    for p in nodes:
        q = xy(p)
        for k, c in enumerate(J):
            if np.linalg.norm(q - c) < MERGE and deg[p] >= 3:
                owner[p] = k; break
        else:
            owner[p] = len(J); J.append(q)
    # joint position: mean of its pixels
    acc = {}
    for p, k in owner.items():
        acc.setdefault(k, []).append(xy(p))
    J = [np.mean(acc[k], axis=0) for k in range(len(J))]
    # ends near a hub go onto its ring
    hub_of = {}
    for k, c in enumerate(J):
        for hi, h in enumerate(hubs):
            d = np.linalg.norm(c - h['c'])
            if d < h['ro'] + S.HUB_MARGIN + HUB_REACH:
                J[k] = np.asarray(h['c']) + (c - h['c']) / d * h['ro']; hub_of[k] = hi
    raw = []
    for ch in chains:
        a, b = owner.get(ch[0]), owner.get(ch[-1])
        if a is None or b is None or a == b:
            continue
        P = np.array([xy(p) for p in ch])
        raw.append([a, b, P, ch])
    # members shorter than SHORT_MEMBER join their two joints (a hub joint stays where it is)
    rep = list(range(len(J)))
    find = lambda k: k if rep[k] == k else find(rep[k])
    for a, b, P, ch in sorted(raw, key=lambda r: len(r[2])):
        ra, rb = find(a), find(b)
        if ra == rb or (ra in hub_of and rb in hub_of):
            continue
        if float(np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1))) < SHORT_MEMBER:
            keep, gone = (ra, rb) if ra in hub_of or (rb not in hub_of and ra < rb) else (rb, ra)
            if keep not in hub_of:
                J[keep] = 0.5 * (J[keep] + J[gone])
            rep[gone] = keep
    M = []
    for a, b, P, ch in raw:
        a, b = find(a), find(b)
        if a == b:
            continue
        P = np.vstack([J[a], P, J[b]])                     # from joint to joint
        Lc = float(np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1)))
        if Lc < MIN_MEMBER:
            continue
        V = np.asarray(LineString(P).simplify(tol).coords)
        # a chain that mostly hugs a hub (the TO's ring is thicker than the hub there) becomes a straight member
        # from its other joint to the hub (it wrapped round the hub as a flat flap in SolidWorks); dropped below
        # if that joint has another member to the hub
        wrap = None
        for end, j, other in ((0, a, b), (-1, b, a)):
            if j in hub_of and other not in hub_of:
                h = hubs[hub_of[j]]
                near = np.linalg.norm(P - np.asarray(h['c']), axis=1) < h['ro'] + HUB_BAND
                if near.mean() > 0.5:
                    V = np.array([J[other], J[j]]) if end == -1 else np.array([J[j], J[other]])
                    wrap = (other, hub_of[j])
        # an end on a hub: the last segment aims at the hub's centre and ends on its ring (a chain meeting the
        # ring at a glancing angle gave end sections standing outside the hub in SolidWorks)
        aim = {end: (np.asarray(hubs[hub_of[j]]['c'], float), hubs[hub_of[j]]['ro'])
               for end, j in ((0, a), (-1, b)) if j in hub_of}
        path = _fillet_path(V, P, tol, aim=aim)
        M.append(dict(i=a, j=b, chain=P, V=V, aim=aim, path=path, ij=[tuple(x) for x in ch], wrap=wrap))
    # a wrapped member whose joint already reaches the same hub by another member is dropped
    reach = {(m['i'] if m['j'] in hub_of else m['j'], hub_of[m['j'] if m['j'] in hub_of else m['i']])
             for m in M if m['wrap'] is None and (m['i'] in hub_of) != (m['j'] in hub_of)}
    M = [m for m in M if m['wrap'] is None or m['wrap'] not in reach]
    # free ends that are not on a hub: leave them (a free member end is a cantilever stub)
    # ---- sections -----------------------------------------------------------------------
    used = sorted({m['i'] for m in M} | {m['j'] for m in M})
    for m in M:
        Pp, s = path_points(m['path'])
        m['length'] = float(s[-1]); m['pts'], m['s'] = Pp, s
    # area-equivalent: shadow outside the hub disks shared out to the nearest path point
    reg = shadow
    for h in hubs:
        reg = reg.difference(Point(*h['c']).buffer(h['ro']))
    mask, (X0, Y0) = S._raster(reg, pix)
    ij = np.argwhere(mask)
    Q = np.column_stack([X0 + (ij[:, 0] + 0.5) * pix, Y0 + (ij[:, 1] + 0.5) * pix])
    allp = np.vstack([m['pts'] for m in M]); lab = np.concatenate([[k] * len(m['pts']) for k, m in enumerate(M)])
    _, nn = cKDTree(allp).query(Q)
    area = np.bincount(lab[nn], minlength=len(M)) * pix * pix
    lo, hi = (to.bounds[0][axis] - 1, to.bounds[1][axis] + 1) if to is not None else (0, 0)
    u_, v_ = [i for i in range(3) if i != axis]
    for k, m in enumerate(M):
        # width: chord of the shadow (outside the hub disks) across the path at stations; the ends are
        # left out (chords there cross the other members of the joint), outliers trimmed, linear fit
        P2, s2 = m['pts'], m['s']
        ts, ws = [], []
        hub_end = (m['i'] in hub_of, m['j'] in hub_of)
        for tt in np.linspace(0.02, 0.98, 25):
            if (tt < 0.15 and not hub_end[0]) or (tt > 0.85 and not hub_end[1]):
                continue
            n_ = np.searchsorted(s2, tt * m['length']).clip(1, len(P2) - 1)
            c, d = P2[n_], P2[n_] - P2[n_ - 1]
            d = d / np.linalg.norm(d); nrm = np.array([-d[1], d[0]])
            seg = LineString([c - 15 * nrm, c + 15 * nrm]).intersection(reg)
            parts = [g for g in getattr(seg, 'geoms', [seg]) if g.length > 0]
            hit = [g for g in parts if g.distance(Point(*c)) < 1e-6]
            if hit:
                ts.append(tt); ws.append(hit[0].length)
        ts, ws = np.array(ts), np.array(ws)
        if len(ws) >= 3 and np.ptp(ts) >= 0.4:
            A = np.column_stack([1 - ts, ts])
            keep = np.ones(len(ws), bool)
            for _ in range(2):
                cf = np.linalg.lstsq(A[keep], ws[keep], rcond=None)[0]
                res = ws - A @ cf; mad = np.median(np.abs(res[keep])) + 1e-6
                keep = np.abs(res) < 3 * 1.4826 * mad
            med = float(np.median(ws))
            w0, w1 = np.clip(cf, 0.6 * med, 1.6 * med)       # a fit over few stations can run away at the ends
        else:
            w0 = w1 = float(np.median(ws)) if len(ws) else S.MIN_WIDTH
        mean = 0.5 * (w0 + w1)
        f = float(np.clip((area[k] / m['length']) / mean, 0.8, 1.5)) if mean > 0 else 1.0
        m['w0'], m['w1'] = max(S.MIN_WIDTH, f * w0), max(S.MIN_WIDTH, f * w1)
        m['area_factor'] = f
        if m['wrap'] is not None:                       # straight to the hub: no taper from the ring's fan
            m['w0'] = m['w1'] = min(m['w0'], m['w1'])
        # depth and centre from rays at stations
        dz, cz, ts = [], [], []
        if to is not None:
            for tt in np.linspace(0.1, 0.9, 9):
                p2 = m['pts'][np.searchsorted(m['s'], tt * m['length']).clip(0, len(m['pts']) - 1)]
                o = np.zeros(3); o[u_], o[v_], o[axis] = p2[0], p2[1], lo
                d3 = np.zeros(3); d3[axis] = 1
                locs = to.ray.intersects_location([o], [d3])[0]
                if len(locs) >= 2:
                    z = np.sort(locs[:, axis]); dz.append(z[-1] - z[0]); cz.append(0.5 * (z[0] + z[-1])); ts.append(tt)
        if m['wrap'] is not None and dz:
            m['d0'] = m['d1'] = float(np.median(dz)); m['centre'] = float(np.mean(cz))
        elif len(dz) >= 2:
            A = np.column_stack([1 - np.array(ts), ts])
            med = float(np.median(dz))
            m['d0'], m['d1'] = (float(x) for x in np.clip(np.linalg.lstsq(A, np.array(dz), rcond=None)[0],
                                                           max(S.MIN_WIDTH, 0.6 * med), 1.6 * med))
            m['centre'] = float(np.mean(cz))
        else:
            m['d0'] = m['d1'] = float(np.mean(dz)) if dz else 1.0
            m['centre'] = float(np.mean(cz)) if cz else 0.0
    # paths again, with arcs no tighter than the section's half-width (and 0.2 mm)
    for m in M:
        m['path'] = _fillet_path(m['V'], m['chain'], tol, 0.5 * max(m['w0'], m['w1']) + 0.2, m['aim'])
        m['length'] = float(path_points(m['path'])[1][-1])
    L = np.array([m['length'] for m in M])
    mid = float(np.average([m['centre'] for m in M], weights=L))
    for m in M:
        for key in ('chain', 'ij', 'pts', 's', 'V', 'wrap', 'aim'):
            m.pop(key)
    return dict(joints=[tuple(map(float, c)) for c in J], hub_of=hub_of, members=M, mid=mid,
                used=used, region=region)
