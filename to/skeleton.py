"""Skeleton of a truss-like TO web: joints, straight bars and bar sizes.

fit(shadow, hubs, to, axis) works in the plane normal to axis:
  1. the web region = the TO's shadow minus the hub disks, rasterized at PIX mm;
  2. its medial axis (skimage skeletonize) as a pixel graph; free branches shorter than SPUR mm are pruned;
  3. each chain between branch points is simplified to straight pieces (Douglas-Peucker, SIMPLIFY mm);
     joints closer than MERGE mm merge;
  4. a chain end within a hub's reach ends on that hub (on the line to its centre, at the ring radius);
  5. per bar: width = 2 x the mean distance to the region's edge along the bar (distance transform),
     depth = the TO's mean extent along the axis over the bar (voxel columns).
Returns dict(joints=[(u, v)], bars=[dict(i, j, width, depth, centre)], hub_of={joint: hub index}).
"""
import math
import numpy as np
import shapely
from scipy import ndimage
from shapely.geometry import Point, LineString
from skimage.morphology import skeletonize

PIX = 0.25       # mm
SPUR = 4.0       # free branches shorter than this are pruned (mm)
SIMPLIFY = 1.0   # Douglas-Peucker tolerance for the bar polylines (mm)
MERGE = 3.0      # joints closer than this merge (mm)
HUB_MARGIN = 1.5 # the hub disk taken out of the web region: ring radius + this (mm)
OPEN = 0.75      # opening radius that removes slivers (mm)
MIN_WIDTH = 1.0  # smallest bar width (mm)
STRAIGHT = 1.0   # a two-bar joint within this of the line through its neighbours is removed (mm)
SHORT = 4.0      # a bar shorter than this contracts into its better-connected end (mm)
REACH = 8.0      # a free bar end connects to a hub or joint within this distance (mm)


def _raster(region, pix):
    x0, y0, x1, y1 = region.bounds
    x0 -= 2 * pix; y0 -= 2 * pix
    nx, ny = int((x1 - x0) / pix) + 4, int((y1 - y0) / pix) + 4
    X, Y = np.meshgrid(x0 + (np.arange(nx) + 0.5) * pix, y0 + (np.arange(ny) + 0.5) * pix, indexing='ij')
    return shapely.contains_xy(region, X, Y), (x0, y0)


def _graph(sk):
    """Pixel skeleton -> nodes (pixels with != 2 neighbours) and chains of pixels between them."""
    pts = set(map(tuple, np.argwhere(sk)))
    nb = lambda p: [(p[0] + dx, p[1] + dy) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
                    if (dx or dy) and (p[0] + dx, p[1] + dy) in pts]
    deg = {p: len(nb(p)) for p in pts}
    nodes = {p for p in pts if deg[p] != 2}
    chains, seen = [], set()
    for n in nodes:
        for q in nb(n):
            if (n, q) in seen:
                continue
            ch = [n, q]; prev = n
            while ch[-1] not in nodes:
                nxt = [r for r in nb(ch[-1]) if r != prev and r not in ch[-2:]]
                if not nxt:
                    break
                # prefer a node if one is adjacent (stops at junction clusters)
                nn = [r for r in nxt if r in nodes]
                prev = ch[-1]; ch.append(nn[0] if nn else nxt[0])
            seen.add((n, ch[1])); seen.add((ch[-1], ch[-2]))
            chains.append(ch)
    # a closed loop without nodes
    return nodes, chains, deg


def _dev(P, a, b):
    """Largest distance of the points P from the line through a and b."""
    w = b - a
    return max(abs(w[0] * (p - a)[1] - w[1] * (p - a)[0]) for p in P) / max(1e-9, np.linalg.norm(w))


def _clean(joints, bars, hub_of, hubs, shadow):
    """Graph clean-up: free ends reach a hub or joint, straight two-bar joints and short bars go."""
    J = [np.array(p, float) for p in joints]
    E = {frozenset(b) for b in bars if b[0] != b[1]}
    nbrs = lambda k: [next(iter(e - {k})) for e in E if k in e]
    inside = lambda a, b: LineString([a, b]).intersection(shadow).length / max(1e-9, np.linalg.norm(a - b))
    # 1. free ends (one bar, not on a hub) connect to the nearest hub ring or joint in the shadow
    for k in range(len(J)):
        if k in hub_of or len(nbrs(k)) != 1:
            continue
        cand = []
        for h_i, h in enumerate(hubs):
            c = np.array(h['c'], float); u = (J[k] - c) / np.linalg.norm(J[k] - c)
            q = c + u * h['ro']
            cand.append((np.linalg.norm(q - J[k]), ('hub', h_i, q)))
        for m in range(len(J)):
            if m != k and m not in nbrs(k) and any(m in e for e in E):
                cand.append((np.linalg.norm(J[m] - J[k]), ('joint', m, J[m])))
        for d, (kind, m, q) in sorted(cand, key=lambda t: t[0]):
            if d > REACH:
                break
            if inside(J[k], q) > 0.9:
                if kind == 'hub':
                    J.append(q); hub_of[len(J) - 1] = m; E.add(frozenset((k, len(J) - 1)))
                else:
                    E.add(frozenset((k, m)))
                break
    # 2./3. repeat: remove straight two-bar joints, contract short bars
    under_ = {}      # bar -> the removed joints it stands for
    changed = True
    while changed:
        changed = False
        for k in range(len(J)):
            nb = nbrs(k)
            if k in hub_of or len(nb) != 2:
                continue
            a, b = J[nb[0]], J[nb[1]]
            e1, e2 = frozenset((k, nb[0])), frozenset((k, nb[1]))
            under = under_.get(e1, []) + under_.get(e2, []) + [J[k]]   # every point the new bar replaces
            d = _dev(under, a, b)
            if d < STRAIGHT and frozenset(nb) not in E and inside(a, b) > 0.9:
                E -= {e1, e2}; E.add(frozenset(nb)); under_[frozenset(nb)] = under; changed = True
        for e in sorted(E, key=lambda e: np.linalg.norm(J[min(e)] - J[max(e)])):
            a, b = sorted(e)
            if np.linalg.norm(J[a] - J[b]) >= SHORT or (a in hub_of and b in hub_of):
                continue
            # keep the hub joint, else the end with more bars; the other one moves into it
            keep, gone = (a, b) if (a in hub_of or (b not in hub_of and len(nbrs(a)) >= len(nbrs(b)))) else (b, a)
            E2, U2, ok = set(), {}, True
            for f in E:
                if f == e:
                    continue
                g = frozenset(keep if x == gone else x for x in f)
                if len(g) == 2:
                    E2.add(g); U2[g] = under_.get(f, []) + ([(J[min(f)] + J[max(f)]) / 2] if g != f else [])
                    if g != f:            # a moved bar must still pass its old mid-point and removed joints
                        p, q = (J[x] for x in g)
                        ok &= _dev(U2[g], p, q) < STRAIGHT
            if not ok:
                continue
            E, under_ = E2, U2; changed = True
            break
    # renumber the joints still in use
    used = sorted({x for e in E for x in e})
    new = {k: i for i, k in enumerate(used)}
    return ([tuple(J[k]) for k in used], [sorted(new[x] for x in e) for e in E],
            {new[k]: h for k, h in hub_of.items() if k in new})


def fit(shadow, hubs, to=None, axis=2, pix=PIX):
    disks = [Point(*h['c']).buffer(h['ro']) for h in hubs]
    region = shadow
    for h in hubs:               # the hub ring and a little more: the shadow's ring is slightly larger
        region = region.difference(Point(*h['c']).buffer(h['ro'] + HUB_MARGIN))
    # an opening removes slivers left along the hub disks
    region = region.buffer(-OPEN).buffer(OPEN).buffer(0)
    mask, (x0, y0) = _raster(region, pix)
    edt = ndimage.distance_transform_edt(mask) * pix
    sk = skeletonize(mask)
    xy = lambda p: (x0 + (p[0] + 0.5) * pix, y0 + (p[1] + 0.5) * pix)
    # prune short free spurs repeatedly
    for _ in range(5):
        nodes, chains, deg = _graph(sk)
        cut = False
        for ch in chains:
            ends_free = [deg[ch[0]] == 1, deg[ch[-1]] == 1]
            L = len(ch) * pix
            near_hub = any(min(Point(*xy(ch[k])).distance(d) for d in disks) < 2.0 for k in (0, -1)) if disks else False
            if any(ends_free) and not all(ends_free) and L < SPUR and not near_hub:
                for p in (ch[:-1] if ends_free[0] else ch[1:]):
                    sk[p] = False
                cut = True
        if not cut:
            break
    nodes, chains, deg = _graph(sk)
    # chains -> simplified polylines in mm
    polys = []
    for ch in chains:
        P = np.array([xy(p) for p in ch])
        if len(P) < 2:
            continue
        ls = LineString(P).simplify(SIMPLIFY)
        polys.append(np.asarray(ls.coords))
    # joints: all polyline vertices, merged within MERGE
    joints, bars = [], []

    def joint(p):
        for k, q in enumerate(joints):
            if math.dist(p, q) < MERGE:
                return k
        joints.append(tuple(p)); return len(joints) - 1
    for P in polys:
        ids = [joint(p) for p in P]
        for a, b in zip(ids[:-1], ids[1:]):
            if a != b and not any({a, b} == {e[0], e[1]} for e in bars):
                bars.append([a, b])
    # chain ends near a hub end on the hub's ring
    hub_of = {}
    for k, p in enumerate(joints):
        for h_i, h in enumerate(hubs):
            d = math.dist(p, h['c'])
            if d < h['ro'] + HUB_MARGIN + MERGE:
                u = (np.array(p) - h['c']) / max(d, 1e-9)
                joints[k] = tuple(np.array(h['c']) + u * h['ro'])
                hub_of[k] = h_i
    joints, bars, hub_of = _clean(joints, bars, hub_of, hubs, shadow)
    # bar sizes
    out = []
    for a, b in bars:
        p, q = np.array(joints[a]), np.array(joints[b])
        L = np.linalg.norm(q - p)
        if L < 1e-6:
            continue
        ts = np.linspace(0.15, 0.85, 15)
        S = p + ts[:, None] * (q - p)
        ij = np.clip(((S - (x0, y0)) / pix).astype(int), 0, np.array(edt.shape) - 1)
        # width from the samples clear of the hub disks (the web region has a hole there)
        clear = [all(math.dist(s, h['c']) > h['ro'] + HUB_MARGIN + OPEN for h in hubs) for s in S]
        e = edt[ij[:, 0], ij[:, 1]][np.array(clear)] if any(clear) else edt[ij[:, 0], ij[:, 1]]
        w = max(MIN_WIDTH, 2 * float(np.median(e)))
        dep, cen = None, None
        if to is not None:
            u_, v_ = [i for i in range(3) if i != axis]
            lo, hi = to.bounds[0][axis] - 1, to.bounds[1][axis] + 1
            ex = []
            for s in S[::3]:
                o = np.zeros(3); o[u_], o[v_], o[axis] = s[0], s[1], lo
                d3 = np.zeros(3); d3[axis] = 1
                locs = to.ray.intersects_location([o], [d3])[0]
                if len(locs) >= 2:
                    z = np.sort(locs[:, axis]); ex.append((z[0], z[-1]))
            if ex:
                ex = np.array(ex); dep = float(np.mean(ex[:, 1] - ex[:, 0])); cen = float(np.mean(ex.mean(1)))
        out.append(dict(i=a, j=b, width=w, depth=dep, centre=cen, length=float(L)))
    _area_widths(out, joints, shadow, hubs, pix)
    return dict(joints=joints, bars=out, hub_of=hub_of, region=region)


def _area_widths(bars, joints, shadow, hubs, pix):
    """Area-equivalent width of each bar: the shadow outside the hub disks is shared out to the nearest
    bar, and the bar's share divided by its length. It counts what the medial width misses: widened
    joints and fans where bars meet a hub."""
    reg = shadow
    for h in hubs:
        reg = reg.difference(Point(*h['c']).buffer(h['ro']))
    mask, (x0, y0) = _raster(reg, pix)
    ij = np.argwhere(mask)
    P = np.column_stack([x0 + (ij[:, 0] + 0.5) * pix, y0 + (ij[:, 1] + 0.5) * pix])
    D = np.empty((len(P), len(bars)))
    for k, b in enumerate(bars):
        a, q = np.array(joints[b['i']]), np.array(joints[b['j']])
        d = q - a; t = np.clip(((P - a) @ d) / (d @ d), 0, 1)
        D[:, k] = np.linalg.norm(P - (a + t[:, None] * d), axis=1)
    owner = np.argmin(D, axis=1)
    area = np.bincount(owner, minlength=len(bars)) * pix * pix
    for k, b in enumerate(bars):
        b['width_area'] = float(area[k] / b['length'])
