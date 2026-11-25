"""3D frame model of a TO result: hubs on any axis (hubs3d.py), straight tapered round members between joints,
a sphere at every joint that is not on a hub, the original's holes cut last.

fit(original, to, detail) works on voxels of the TO (pitch PITCH x the part size, at least MIN_PITCH):
  1. hubs (hubs3d.detect); the web = the TO voxels outside the hub cylinders;
  2. mirror symmetry: each axis-aligned plane through the original's box centre about which the TO voxels are
     symmetric (overlap >= SYM); the voxels are made symmetric (majority over the mirror images);
  3. the web's medial curves (skimage skeletonize, 3D) as a graph: junction voxels in clusters are one joint,
     chains of voxels between joints are edges; free branches shorter than SPUR are pruned, joints of degree 2
     dissolved, joints closer than the sum of their radii (x MERGE) merged;
  4. a chain end whose inscribed tube comes within REACH of a hub's cylinder ends on that hub; a chain passing
     that close elsewhere gets a joint there that ends on the hub; such a hub joint (anchor) is moved to the
     nearest point inside the hub at least ANCHOR_M x its inscribed radius from every face; anchors closer than
     ANCHOR_MERGE are one;
  5. each chain is simplified to straight pieces (Douglas-Peucker, tolerance by detail level); a bend is a joint;
  6. a symmetric part is made exactly symmetric: joints near a mirror plane put on it, the members in the
     fundamental domain kept and copied by the mirrors (orbits of joints share their dimensions in SolidWorks);
  7. member radii: the TO's distance transform along the piece, a straight-line fit giving the radius at each
     end, scaled so that the member's volume matches the TO voxels nearest to it (area-equivalent), capped at
     the depth of the member's axis inside the original; radii in classes (detail level), mirror images equal;
  8. sphere at each free joint: 1.1 x the widest member end there, at least the TO's own inscribed radius there
     (at most NODE_CAP x the widest end). Members stop TRIM_F x the node radius short of a joint's centre; at a
     hub they end at their anchor, a point inside the hub at least the member's end radius from every face of
     the hub (members sharing an anchor stop HUB_SHARE_F x its depth short of it), so no two member ends
     coincide and every end section lies inside a node or a hub.
Returns dict(hubs, joints=[xyz], kind=['free'|'hub'], hub_of={joint: hub}, members=[dict(i, j, r0, r1, s0, s1)],
nodes=[(joint, R)], sym=[axes], centre, orbit=[orbit id per joint], detail).
"""
import math
import numpy as np
import networkx as nx
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage.morphology import skeletonize
import cadquery as cq
from cadquery import Vector

import hubs3d

PITCH = 1 / 200          # voxel pitch as a fraction of the box diagonal
MIN_PITCH = 0.25         # mm
SYM = 0.95               # mirror symmetry: overlap of the TO voxels with their mirror image
REACH = 1.5              # mm (or 3 voxels), a chain end this close to a hub's cylinder ends on the hub
MERGE = 0.8              # joints closer than MERGE x the sum of their inscribed radii are one joint
ANCHOR_MERGE = 1.0       # mm, anchors on one hub closer than this are one
ANCHOR_M = 1.5           # an anchor is this x its inscribed radius (+ 0.3 mm) inside every face of its hub
NODE_F = 1.1             # node radius over the widest member end radius there
NODE_CAP = 1.3           # the node radius at most this x the widest member end (the TO's own size there)
TRIM_F = 0.3             # members stop this x the node radius short of a joint's centre
HUB_SHARE_F = 0.25       # members sharing an anchor stop this x the anchor's depth in the hub short of it
AREA_CAP = 2.0           # area-equivalent scaling of a member's radii at most this factor
MIN_R = 0.6              # mm, thinnest member radius
STUBBY, STUBBY_SHARE = 1.5, 0.5   # warn when over half the members are shorter than 1.5 x the sum of their end radii
# detail level: Douglas-Peucker tolerance (x the mean member radius), spur length (x the junction radius),
# member radius classes, node radius classes (0: every value its own), significant digits
DETAIL = {1: dict(tol=2.0, spur=3.0, rclass=2, nclass=1, digits=2),
          2: dict(tol=1.5, spur=2.5, rclass=3, nclass=2, digits=2),
          3: dict(tol=1.0, spur=2.0, rclass=4, nclass=2, digits=3),
          4: dict(tol=0.7, spur=1.5, rclass=5, nclass=3, digits=3),
          5: dict(tol=0.4, spur=1.2, rclass=0, nclass=0, digits=3)}

OFFS = np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)])


def sig(x, n):
    if x == 0:
        return 0.0
    return round(x, max(0, n - 1 - int(math.floor(math.log10(abs(x))))))


# ---- voxels ----------------------------------------------------------------------------------------------
def voxelize(to, pitch):
    v = to.voxelized(pitch).fill()
    return v.matrix.copy(), np.asarray(v.transform[:3, 3], float)


def mirror_index(n, o, pitch, c):
    """Index map i -> i' of the mirror about the plane at c along one axis (nearest voxel), -1 outside."""
    i = np.arange(n)
    j = np.rint((2 * c - (o + i * pitch) - o) / pitch).astype(int)
    j[(j < 0) | (j >= n)] = -1
    return j


def mirrored(M, ax, j):
    out = np.zeros_like(M)
    ok = j >= 0
    src = np.moveaxis(M, ax, 0)
    dst = np.moveaxis(out, ax, 0)
    dst[ok] = src[j[ok]]
    return out


def symmetry(M, o, pitch, centre):
    """Axes about which the voxels are mirror symmetric (plane through centre), and the symmetrized voxels."""
    axes, maps = [], {}
    for ax in range(3):
        j = mirror_index(M.shape[ax], o[ax], pitch, centre[ax])
        Mm = mirrored(M, ax, j)
        ov = (M & Mm).sum() / max(M.sum(), 1)
        if ov >= SYM:
            axes.append(ax); maps[ax] = j
    if not axes:
        return axes, M
    imgs = [M.astype(np.int16)]
    for k in range(1, 2 ** len(axes)):              # every element of the group the mirrors generate
        im = M
        for b, ax in enumerate(axes):
            if k >> b & 1:
                im = mirrored(im, ax, maps[ax])
        imgs.append(im.astype(np.int16))
    S = np.sum(imgs, axis=0)
    return axes, S * 2 >= len(imgs)


# ---- skeleton graph --------------------------------------------------------------------------------------
def skeleton_graph(S):
    """Voxel graph of a skeleton (26-connected): networkx graph on voxel index tuples."""
    idx = [tuple(p) for p in np.argwhere(S)]
    pset = set(idx)
    G = nx.Graph()
    G.add_nodes_from(idx)
    for p in idx:
        for d in OFFS:
            q = (p[0] + d[0], p[1] + d[1], p[2] + d[2])
            if q in pset:
                G.add_edge(p, q)
    return G


def chains(G):
    """Contracted graph: joints = clusters of voxels of degree != 2 (ends and junction clusters); edges = voxel
    chains between them. Returns (joint voxel lists, edges [(a, b, [voxels from a to b])])."""
    special = {p for p in G if G.degree(p) != 2}
    H = G.subgraph(special)
    clusters = [sorted(c) for c in nx.connected_components(H)]
    cid = {p: k for k, c in enumerate(clusters) for p in c}
    edges, seen = [], set()
    for k, c in enumerate(clusters):
        for p in c:
            for q in G.neighbors(p):
                if q in special or (p, q) in seen:
                    continue
                path = [p, q]; prev, cur = p, q
                seen.add((p, q))
                while cur not in special:
                    nb = [w for w in G.neighbors(cur) if w != prev]
                    if not nb:
                        break
                    prev, cur = cur, nb[0]
                    path.append(cur)
                if cur in special:
                    seen.add((cur, prev))
                    edges.append((k, cid[cur], path))
    # pure loops (no special voxel): one voxel made a joint
    done = {p for _, _, pa in edges for p in pa} | special
    for comp in nx.connected_components(G):
        if comp & done:
            continue
        cyc = list(comp)
        p0 = cyc[0]
        clusters.append([p0]); k = len(clusters) - 1
        order = list(nx.dfs_preorder_nodes(G.subgraph(comp), p0)) + [p0]
        edges.append((k, k, order))
    return clusters, edges


class Frame:
    """Contracted skeleton graph in model mm: joints (position, inscribed radius, hub index or None) and chains
    (joint a, joint b, points from a to b, inscribed radius at each point)."""

    def __init__(self):
        self.P, self.R, self.hub = {}, {}, {}
        self.E = {}                                   # edge id -> [a, b, points (n x 3), radii (n)]
        self._n = 0
        self._e = 0

    def add_joint(self, p, r, hub=None):
        k = self._n; self._n += 1
        self.P[k], self.R[k], self.hub[k] = np.asarray(p, float), float(r), hub
        return k

    def add_edge(self, a, b, pts, rr):
        e = self._e; self._e += 1
        self.E[e] = [a, b, np.asarray(pts, float), np.asarray(rr, float)]
        return e

    def deg(self, k):
        return sum((a == k) + (b == k) for a, b, _, _ in self.E.values())

    def edges_of(self, k):
        return [e for e, (a, b, _, _) in self.E.items() if a == k or b == k]

    @staticmethod
    def length(pts):
        return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum()) if len(pts) > 1 else 0.0

    def oriented(self, e, k):
        """(points, radii) of edge e starting at joint k."""
        a, b, pts, rr = self.E[e]
        return (pts, rr) if a == k else (pts[::-1], rr[::-1])

    def drop_joint(self, k):
        for d in (self.P, self.R, self.hub):
            d.pop(k, None)

    def set_end(self, k):
        """Chain ends at joint k moved to its position."""
        for e in self.edges_of(k):
            a, b, pts, rr = self.E[e]
            pts = pts.copy()
            if a == k:
                pts[0] = self.P[k]
            if b == k:
                pts[-1] = self.P[k]
            self.E[e][2] = pts

    def dissolve(self):
        """Joints of degree 2 (not on a hub) removed: their two chains joined."""
        changed = True
        while changed:
            changed = False
            for k in list(self.P):
                es = self.edges_of(k)
                if self.hub[k] is not None or len(es) != 2 or es[0] == es[1]:
                    continue
                e1, e2 = es
                (p1, r1), (p2, r2) = self.oriented(e1, k), self.oriented(e2, k)
                a = self.E[e1][0] if self.E[e1][1] == k else self.E[e1][1]
                b = self.E[e2][1] if self.E[e2][0] == k else self.E[e2][0]
                del self.E[e1], self.E[e2]
                self.add_edge(a, b, np.vstack([p1[::-1], p2[1:]]), np.concatenate([r1[::-1], r2[1:]]))
                self.drop_joint(k)
                changed = True

    def merge(self, a, b):
        """Joint b merged into a: edges between them dropped, b's edges moved to a."""
        for e in list(self.E):
            x, y, pts, rr = self.E[e]
            if {x, y} == {a, b}:
                del self.E[e]; continue
            if x == b:
                self.E[e][0] = a
            if y == b:
                self.E[e][1] = a
        self.drop_joint(b)
        self.set_end(a)

    def split(self, e, i, hub=None):
        """Edge e split at its point i: a new joint there."""
        a, b, pts, rr = self.E[e]
        k = self.add_joint(pts[i], rr[i], hub)
        del self.E[e]
        self.add_edge(a, k, pts[:i + 1], rr[:i + 1])
        self.add_edge(k, b, pts[i:], rr[i:])
        return k


def build_frame(S, o, pitch, D):
    """Frame of a skeleton mask S (voxel origin o, pitch), inscribed radii from the distance map D (mm)."""
    G = skeleton_graph(S)
    clusters, edges = chains(G)
    F = Frame()
    ids = []
    for c in clusters:
        c = np.array(c)
        ids.append(F.add_joint(o + c.mean(axis=0) * pitch, float(np.max(D[tuple(c.T)]))))
    for a, b, path in edges:
        path = np.array(path)
        pts = o + path.astype(float) * pitch
        rr = D[tuple(path.T)]
        pts[0], pts[-1] = F.P[ids[a]], F.P[ids[b]]
        if ids[a] == ids[b] and len(pts) < 4:
            continue
        F.add_edge(ids[a], ids[b], pts, rr)
    return F


def dist_cyl(h, P):
    """Distance of points P (n x 3) from hub h's cylinder (negative inside)."""
    k = h['axis']; u, v = hubs3d.other(k)
    rad = np.hypot(P[:, u] - h['c'][0], P[:, v] - h['c'][1]) - h['ro']
    ax = np.maximum(h['lo'] - P[:, k], P[:, k] - h['hi'])
    out = np.hypot(np.maximum(rad, 0), np.maximum(ax, 0))
    return np.where((rad <= 0) & (ax <= 0), np.maximum(rad, ax), out)


def attach_hubs(F, hubs, reach):
    """Contacts with the hubs, by the tube of inscribed radius round the skeleton (g = distance from the hub's
    cylinder minus the inscribed radius): a chain end with g <= reach ends on the hub; a chain passing within
    reach elsewhere (the arms meeting under a boss) gets a joint at its closest point and a spoke to the hub."""
    for k in list(F.P):
        if F.deg(k) != 1 or F.hub[k] is not None:
            continue
        g = [float(dist_cyl(h, F.P[k][None])[0]) - F.R[k] for h in hubs]
        n = int(np.argmin(g))
        if g[n] <= reach:
            F.hub[k] = n
    for n, h in enumerate(hubs):
        for e in list(F.E):
            if e not in F.E:
                continue
            a, b, pts, rr = F.E[e]
            if F.hub[a] == n or F.hub[b] == n:
                continue
            g = dist_cyl(h, pts) - rr
            near = g <= reach
            if not near.any():
                continue
            i = int(np.argmin(g))
            L = len(pts)
            if i <= 1 or i >= L - 2:                  # at an end joint: the spoke starts there
                k = a if i <= 1 else b
            else:
                k = F.split(e, i)
            q = F.P[k]
            hk = F.add_joint(q, F.R[k], n)            # moved onto the axis by anchor()
            F.add_edge(k, hk, np.vstack([q, q]), np.array([F.R[k], F.R[k]]))


def prune(F, spur):
    """Free branches (an end joint not on a hub) shorter than spur x the radius at their junction removed, small
    loops and pieces with no hub removed; repeated until nothing changes."""
    changed = True
    while changed:
        changed = False
        for e in list(F.E):
            if e not in F.E:
                continue
            a, b, pts, rr = F.E[e]
            L = F.length(pts)
            if a == b:
                if L < 4 * math.pi * max(F.R[a], 1.0):
                    del F.E[e]; changed = True
                continue
            for end, oth in ((a, b), (b, a)):
                if F.deg(end) == 1 and F.hub[end] is None and (F.deg(oth) > 1 or F.hub[oth] is not None) \
                        and L < spur * max(F.R[oth], 1.0):
                    del F.E[e]; F.drop_joint(end); changed = True
                    break
        for k in list(F.P):
            if F.deg(k) == 0:
                F.drop_joint(k); changed = True
        F.dissolve()
    G = nx.Graph()
    G.add_nodes_from(F.P)
    G.add_edges_from((a, b) for a, b, _, _ in F.E.values())
    for comp in list(nx.connected_components(G)):
        if not any(F.hub[k] is not None for k in comp):
            for e in [e for e, v in F.E.items() if v[0] in comp]:
                del F.E[e]
            for k in comp:
                F.drop_joint(k)


def merge_close(F, f=MERGE):
    """Joints joined by a chain shorter than f x the sum of their radii are one joint (a hub joint wins, then the
    joint of higher degree, which keeps its position)."""
    changed = True
    while changed:
        changed = False
        for e in sorted(F.E, key=lambda e: F.length(F.E[e][2])):
            a, b, pts, rr = F.E[e]
            if a == b or F.length(pts) >= f * (F.R[a] + F.R[b]):
                continue
            if F.hub[a] is not None and F.hub[b] is not None:
                continue
            keep, gone = (a, b) if (F.hub[a] is not None, F.deg(a)) >= (F.hub[b] is not None, F.deg(b)) else (b, a)
            F.R[keep] = max(F.R[keep], F.R[gone])
            F.merge(keep, gone)
            changed = True
            break
        F.dissolve()


def inner(h, q, m):
    """Point of hub h's cylinder shrunk by m on every face (radius ro - m, ends lo + m and hi - m) nearest to q."""
    k = h['axis']; u, v = hubs3d.other(k)
    m = min(m, 0.9 * h['ro'], 0.45 * (h['hi'] - h['lo']))
    p = np.array(q, float)
    p[k] = min(max(p[k], h['lo'] + m), h['hi'] - m)
    d = np.array([p[u] - h['c'][0], p[v] - h['c'][1]])
    rho = float(np.linalg.norm(d))
    if rho > h['ro'] - m:
        d = d / rho * (h['ro'] - m)
        p[u], p[v] = h['c'][0] + d[0], h['c'][1] + d[1]
    return p


def anchor(F, hubs):
    """Hub joints moved inside their hub: the nearest point at least ANCHOR_M x their inscribed radius (plus
    0.3 mm) from every face of the hub, so a member's end section there is inside the hub whatever its
    direction; hub joints on one hub closer than ANCHOR_MERGE merged."""
    for k in list(F.P):
        n = F.hub.get(k)
        if n is None:
            continue
        F.P[k] = inner(hubs[n], F.P[k], ANCHOR_M * F.R[k] + 0.3)
        F.set_end(k)
    changed = True
    while changed:
        changed = False
        ks = [k for k in F.P if F.hub[k] is not None]
        for i, a in enumerate(ks):
            for b in ks[i + 1:]:
                if F.hub[a] == F.hub[b] and np.linalg.norm(F.P[a] - F.P[b]) < ANCHOR_MERGE:
                    F.P[a] = 0.5 * (F.P[a] + F.P[b]); F.R[a] = max(F.R[a], F.R[b])
                    F.merge(a, b)
                    changed = True
                    break
            if changed:
                break
    for e in list(F.E):                               # a chain from a hub to itself
        a, b, _, _ = F.E[e]
        if a == b and F.hub[a] is not None:
            del F.E[e]


def dp(pts, tol):
    """Douglas-Peucker: indices of the kept points (first and last always)."""
    keep = [0, len(pts) - 1]
    stack = [(0, len(pts) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, b = pts[i], pts[j]
        d = b - a
        L = np.linalg.norm(d)
        seg = pts[i + 1:j]
        if L < 1e-9:
            dist = np.linalg.norm(seg - a, axis=1)
        else:
            dist = np.linalg.norm(np.cross(seg - a, d / L), axis=1)
        m = int(np.argmax(dist))
        if dist[m] > tol:
            k = i + 1 + m
            keep.append(k)
            stack += [(i, k), (k, j)]
    return sorted(set(keep))


def straighten(F, tolf, pitch):
    """Members: each chain cut into straight pieces (Douglas-Peucker, tolerance tolf x the chain's mean radius,
    at least 1.5 voxels); bends become joints. Returns members [dict(i, j, pts, rr)]."""
    members = []
    for e, (a, b, pts, rr) in list(F.E.items()):
        tol = max(tolf * float(np.mean(rr)), 1.5 * pitch)
        idx = dp(pts, tol)
        if a == b and len(idx) < 4:                   # a loop needs at least a triangle
            n = len(pts)
            idx = sorted(set([0, n // 3, 2 * n // 3, n - 1]))
        js = [a]
        for i in idx[1:-1]:
            js.append(F.add_joint(pts[i], rr[i]))
        js.append(b)
        for s, (i0, i1) in enumerate(zip(idx[:-1], idx[1:])):
            members.append(dict(i=js[s], j=js[s + 1], pts=pts[i0:i1 + 1], rr=rr[i0:i1 + 1]))
    return members


def contract_short(F, members, f=MERGE):
    """Members shorter than f x the sum of their joints' inscribed radii contracted (bends close together, a
    bend next to a junction); then members with the same two joints joined (area kept)."""
    changed = True
    while changed:
        changed = False
        for m in sorted(members, key=lambda m: np.linalg.norm(F.P[m['i']] - F.P[m['j']])):
            a, b = m['i'], m['j']
            if a == b or np.linalg.norm(F.P[a] - F.P[b]) >= f * (F.R[a] + F.R[b]):
                continue
            if F.hub[a] is not None and F.hub[b] is not None:
                continue
            da = sum((x['i'] == a) + (x['j'] == a) for x in members)
            db = sum((x['i'] == b) + (x['j'] == b) for x in members)
            keep, gone = (a, b) if (F.hub[a] is not None, da) >= (F.hub[b] is not None, db) else (b, a)
            members.remove(m)
            for x in members:
                if x['i'] == gone:
                    x['i'] = keep
                if x['j'] == gone:
                    x['j'] = keep
            F.R[keep] = max(F.R[keep], F.R[gone])
            F.drop_joint(gone)
            changed = True
            break
    members[:] = [m for m in members if m['i'] != m['j']]
    out = {}
    for m in members:
        key = tuple(sorted((m['i'], m['j'])))
        if key in out:
            o = out[key]
            pts = m['pts'] if m['i'] == o['i'] else m['pts'][::-1]
            rr = m['rr'] if m['i'] == o['i'] else m['rr'][::-1]
            o['extra'] = o.get('extra', []) + [(pts, rr)]
        else:
            out[key] = m
    members[:] = list(out.values())


def ops(axes):
    """Elements of the group the mirrors generate: tuples of the flipped axes (the identity first)."""
    return [tuple(ax for b, ax in enumerate(axes) if k >> b & 1) for k in range(2 ** len(axes))]


def apply(g, p, centre):
    q = np.array(p, float)
    for ax in g:
        q[ax] = 2 * centre[ax] - q[ax]
    return q


def mirror_regenerate(F, members, axes, centre, tol, hubs):
    """Exact mirror symmetry: joints within tol of a mirror plane put on it; the members in the fundamental domain
    (midpoint on the low side of every mirror plane) kept, a member crossing a plane made symmetric (its far end
    the mirror image of its near end); every member then copied by every mirror of the group, joints at one
    position made one. Returns {joint: (orbit, op)}: op takes the orbit's representative to the joint."""
    used = {m['i'] for m in members} | {m['j'] for m in members}
    if not axes:
        return {k: (k, ()) for k in used}
    G = ops(axes)
    c = np.asarray(centre, float)
    for k in used:
        for ax in axes:
            if abs(F.P[k][ax] - c[ax]) < tol:
                F.P[k][ax] = c[ax]
    low = lambda p: all(p[ax] <= c[ax] + 1e-6 for ax in axes)
    keep = [m for m in members if low(0.5 * (F.P[m['i']] + F.P[m['j']]))]
    for m in keep:                                    # a member crossing a plane: ends mirror images
        for ax in axes:
            a, b = F.P[m['i']], F.P[m['j']]
            if (a[ax] - c[ax]) * (b[ax] - c[ax]) < 0:
                lo_, hi_ = (m['i'], m['j']) if a[ax] < b[ax] else (m['j'], m['i'])
                im = apply((ax,), F.P[lo_], c)
                if np.linalg.norm(im - F.P[hi_]) < 2 * tol:
                    F.P[hi_] = im
    P0, R0, H0 = dict(F.P), dict(F.R), dict(F.hub)
    F.P, F.R, F.hub = {}, {}, {}
    orb, where = {}, []

    def joint(k, g):
        p = apply(g, P0[k], c)
        for q, kk in where:
            if np.linalg.norm(q - p) < 1e-6:
                return kk
        h = None
        if H0[k] is not None:
            h = min(range(len(hubs)), key=lambda n: float(dist_cyl(hubs[n], p[None])[0]))
        kk = F.add_joint(p, R0[k], h)
        where.append((p, kk))
        orb[kk] = (k, g)
        return kk
    out, seen = [], set()
    for g in G:
        for m in keep:
            a, b = joint(m['i'], g), joint(m['j'], g)
            key = tuple(sorted((a, b)))
            if a == b or key in seen:
                continue
            seen.add(key)
            mm = dict(m, i=a, j=b, pts=np.array([apply(g, p, c) for p in m['pts']]))
            if 'extra' in m:
                mm['extra'] = [(np.array([apply(g, p, c) for p in pts]), rr) for pts, rr in m['extra']]
            out.append(mm)
    members[:] = out
    # orbits numbered by their representative
    reps = sorted({r for r, _ in orb.values()})
    num = {r: n for n, r in enumerate(reps)}
    return {k: (num[r], g) for k, (r, g) in orb.items()}


def round_joints(F, orb, hubs, axes, centre, grid):
    """Joint coordinates on a grid (distances from the mirror planes along symmetric axes, so mirror images stay
    images); a hub joint keeps its axis line and is rounded along the axis only."""
    reps = {}
    for k, (n, g) in orb.items():
        p = apply(g, F.P[k], centre)                  # back to the orbit's representative
        reps.setdefault(n, p)
    for n, p in reps.items():
        q = p.copy()
        for ax in range(3):
            q[ax] = centre[ax] + np.round((p[ax] - centre[ax]) / grid) * grid if ax in axes else np.round(p[ax] / grid) * grid
        reps[n] = q
    for k, (n, g) in orb.items():
        F.P[k] = apply(g, reps[n], centre)


def coalesce(F, members, orb, tol):
    """Joints that symmetrizing and rounding put at one point (within tol) are one joint (a hub joint wins);
    members between them dropped, members with the same two joints joined (contract_short)."""
    ks = sorted(orb, key=lambda k: (F.hub[k] is None, k))
    gone = {}
    for a_, a in enumerate(ks):
        if a in gone:
            continue
        for b in ks[a_ + 1:]:
            if b not in gone and np.linalg.norm(F.P[a] - F.P[b]) < tol:
                if F.hub[a] is not None and F.hub[b] is not None and F.hub[a] != F.hub[b]:
                    continue
                gone[b] = a
                F.R[a] = max(F.R[a], F.R[b])
    for m in members:
        m['i'] = gone.get(m['i'], m['i'])
        m['j'] = gone.get(m['j'], m['j'])
    for b in gone:
        F.drop_joint(b); orb.pop(b, None)
    contract_short(F, members, f=0.0)


def _fit_line(s, r):
    if len(s) >= 3 and np.ptp(s) > 1e-6:
        b, a = np.polyfit(s, r, 1)
        return a, b
    return float(np.mean(r)), 0.0


def size_members(F, members, hubs, web_pts, vox_vol, depth, pitch):
    """Radii r0 (at joint i) and r1 (at joint j): inscribed radii along the member, straight-line fit; scaled to
    the volume of the TO voxels nearest the member (area-equivalent, factor in [0.8, AREA_CAP]); taper at most 2;
    capped at the depth of the member's axis inside the original."""
    samp, lab = [], []
    for n, m in enumerate(members):
        A, B = F.P[m['i']], F.P[m['j']]
        L = np.linalg.norm(B - A)
        t = np.linspace(0, 1, max(2, int(L / (0.5 * pitch)) + 1))
        samp.append(A + t[:, None] * (B - A)); lab.append(np.full(len(t), n))
    samp, lab = np.vstack(samp), np.concatenate(lab)
    _, nn = cKDTree(samp).query(web_pts)
    vol = np.bincount(lab[nn], minlength=len(members)) * vox_vol
    for n, m in enumerate(members):
        A, B = F.P[m['i']], F.P[m['j']]
        L = float(np.linalg.norm(B - A)); t = (B - A) / L
        ss, rr = [], []
        for pts, r_ in [(m['pts'], m['rr'])] + m.get('extra', []):
            s = (pts - A) @ t
            ok = (s > 0.5 * F.R[m['i']]) & (s < L - 0.5 * F.R[m['j']])
            if ok.sum() < 2:
                ok = np.ones(len(s), bool)
            ss.append(s[ok]); rr.append(r_[ok])
        a, b = _fit_line(np.concatenate(ss), np.concatenate(rr))
        r0, r1 = max(a, MIN_R), max(a + b * L, MIN_R)
        if len(m.get('extra', [])):                   # parallel chains joined: their areas add
            k = math.sqrt(1 + len(m['extra']))
            r0, r1 = k * r0, k * r1
        lo_, hi_ = sorted((r0, r1))
        if hi_ > 2 * lo_:
            mean = 0.5 * (lo_ + hi_)
            small, big = 2 * mean / 3, 4 * mean / 3
            r0, r1 = (small, big) if r0 < r1 else (big, small)
        # length outside the hubs, and the depth of the axis inside the original there
        P = A + np.linspace(0, 1, 41)[:, None] * (B - A)
        outside = np.all([dist_cyl(h, P) > 0 for h in hubs], axis=0) if hubs else np.ones(len(P), bool)
        Leff = L * outside.mean()
        vcur = math.pi / 3 * Leff * (r0 * r0 + r0 * r1 + r1 * r1)
        f = math.sqrt(vol[n] / vcur) if vcur > 0 and vol[n] > 0 else 1.0
        f = min(max(f, 0.8), AREA_CAP)
        r0, r1 = f * r0, f * r1
        m['r0'], m['r1'], m['L'] = r0, r1, L
        if depth is None:
            continue
        dep = depth(P)
        half = np.arange(len(P)) < len(P) // 2
        for end, sel in ((0, half & outside), (1, ~half & outside)):
            if sel.any():
                cap = max(float(dep[sel].min()), MIN_R)
                if end == 0:
                    r0 = min(r0, cap)
                else:
                    r1 = min(r1, cap)
        m['r0'], m['r1'], m['L'] = r0, r1, L


def push_inside(F, members, orb, axes, centre, om, rounds=3):
    """Free joints moved into the original until a node of NODE_F x the widest member end there fits (TO members
    lie flat against the original's faces; a round member on their medial curve stuck out or was thinned).
    One representative per orbit is moved, along the inward normal of the nearest face of the original
    (within a mirror plane for a joint on it); the others are its mirror images."""
    import trimesh
    reps = {n: k for k, (n, g) in orb.items() if g == () and F.hub[k] is None}
    for _ in range(rounds):
        moved = 0
        for n, k in reps.items():
            ends = [m['r0'] for m in members if m['i'] == k] + [m['r1'] for m in members if m['j'] == k]
            if not ends:
                continue
            need = NODE_F * max(ends)
            p = F.P[k]
            C, dist, tri = trimesh.proximity.closest_point(om, p[None])
            inside = om.contains(p[None])[0]
            dep = dist[0] if inside else -dist[0]
            if dep >= need - 1e-3:
                continue
            nin = -om.face_normals[tri[0]]
            if inside and dist[0] > 1e-3:
                nin = (p - C[0]) / dist[0]
            for ax in axes:                           # a joint on a mirror plane stays on it
                if abs(p[ax] - centre[ax]) < 1e-9:
                    nin[ax] = 0.0
            if np.linalg.norm(nin) < 1e-6:
                continue
            nin = nin / np.linalg.norm(nin)
            F.P[k] = p + (need - dep) * nin
            moved += 1
        for k, (n, g) in orb.items():                 # mirror images follow their representative
            r = reps.get(n)
            if r is not None and k != r and F.hub[k] is None:
                F.P[k] = apply(g, F.P[r], centre)
        if not moved:
            break


def symmetric_radii(F, members, orb, axes, centre):
    """Member end radii equal across mirror images (mean over each class of (member, end))."""
    if not axes:
        return
    G = ops(axes)
    pos = {k: F.P[k] for k in orb}
    keys = {}
    for n, m in enumerate(members):
        keys[tuple(sorted((m['i'], m['j'])))] = n
    X = np.array([pos[k] for k in orb]); kk = list(orb)
    tree = cKDTree(X)
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for n, m in enumerate(members):
        for g in G[1:]:
            di, ii = tree.query(apply(g, pos[m['i']], centre))
            dj, jj = tree.query(apply(g, pos[m['j']], centre))
            if di > 1e-6 or dj > 1e-6:
                continue
            gi, gj = kk[ii], kk[jj]
            n2 = keys.get(tuple(sorted((gi, gj))))
            if n2 is None:
                continue
            parent[find((n, m['i']))] = find((n2, gi))
            parent[find((n, m['j']))] = find((n2, gj))
    cls = {}
    for n, m in enumerate(members):
        for k, key in ((m['i'], 'r0'), (m['j'], 'r1')):
            cls.setdefault(find((n, k)), []).append((n, key))
    for c in cls.values():
        v = np.mean([members[n][key] for n, key in c])
        for n, key in c:
            members[n][key] = v


def nodes(F, members, depth, D):
    """Node radius at each free joint: NODE_F x the widest member end there, at least the TO's inscribed radius
    there (at most NODE_CAP x the widest end), at most the depth inside the original (member ends then thinned)."""
    out = {}
    for k in F.P:
        if F.hub[k] is not None:
            continue
        ends = [m['r0'] for m in members if m['i'] == k] + [m['r1'] for m in members if m['j'] == k]
        if not ends:
            continue
        rmax = max(ends)
        R = max(NODE_F * rmax, min(F.R[k], NODE_CAP * rmax))
        dep = max(float(depth(F.P[k][None])[0]), MIN_R * NODE_F)
        if R > dep:
            R = dep
            for m in members:
                if m['i'] == k:
                    m['r0'] = min(m['r0'], R / NODE_F)
                if m['j'] == k:
                    m['r1'] = min(m['r1'], R / NODE_F)
        out[k] = R
    return out


def hub_depth(h, p):
    """Distance of point p inside hub h's cylinder from its nearest face."""
    k = h['axis']; u, v = hubs3d.other(k)
    return min(h['ro'] - math.hypot(p[u] - h['c'][0], p[v] - h['c'][1]), p[k] - h['lo'], h['hi'] - p[k])


def classes(vals, w, k, digits):
    import swbars
    x = swbars.classes(np.asarray(vals, float), np.asarray(w, float), k)
    return [sig(v, digits) for v in x]


def fit(original, to, detail=3, log=print):
    """The frame model of a TO result (see the module docstring)."""
    DT = DETAIL[detail]
    hubs = hubs3d.detect(original, to)
    if not hubs:
        raise ValueError('no hubs: the original has no hole that the TO surrounds')
    log('hubs: ' + '; '.join('axis %s at (%.2f, %.2f), r %.2f-%.2f, %.2f to %.2f mm' % (
        'xyz'[h['axis']], *h['c'], h['ri'], h['ro'], h['lo'], h['hi']) for h in hubs))
    diag = float(np.linalg.norm(np.ptp(to.bounds, axis=0)))
    pitch = max(MIN_PITCH, round(diag * PITCH * 4) / 4)
    M, o = voxelize(to, pitch)
    bb = original.BoundingBox()
    centre = np.array([(bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2, (bb.zmin + bb.zmax) / 2])
    axes, Ms = symmetry(M, o, pitch, centre)
    log('voxels: pitch %.2f mm, %d filled; mirror symmetry about %s' % (
        pitch, M.sum(), ', '.join('%s = %.4g' % ('xyz'[a], centre[a]) for a in axes) or 'none'))
    I = np.argwhere(Ms); P = o + I * pitch
    inhub = np.zeros(len(P), bool)
    for h in hubs:
        inhub |= hubs3d.inside(h, P, margin=1.5 * pitch)
    W = Ms.copy(); W[tuple(I[inhub].T)] = False
    D = ndimage.distance_transform_edt(np.pad(Ms, 1))[1:-1, 1:-1, 1:-1] * pitch
    S = skeletonize(W) > 0
    F = build_frame(S, o, pitch, D)
    reach = max(REACH, 3 * pitch)
    attach_hubs(F, hubs, reach)
    for _ in range(2):
        prune(F, DT['spur'])
        merge_close(F)
    anchor(F, hubs)
    members = straighten(F, DT['tol'], pitch)
    contract_short(F, members)
    tol = max(2.0, 3 * pitch)
    orb = mirror_regenerate(F, members, axes, centre, tol, hubs)
    grid = 0.5 if DT['digits'] <= 2 else 0.1
    round_joints(F, orb, hubs, axes, centre, grid)
    coalesce(F, members, orb, 0.5 * grid + 1e-9)
    # depth inside the original, from a fine tessellation of it
    import trimesh
    vs, fs = original.tessellate(0.02, 0.1)
    om = trimesh.Trimesh(np.array([v.toTuple() for v in vs]), np.array(fs))
    depth = lambda Q: trimesh.proximity.signed_distance(om, np.asarray(Q, float))
    web_pts = o + np.argwhere(W) * pitch
    corr = to.volume / (Ms.sum() * pitch ** 3)       # filled voxels overstate the volume (boundary voxels)
    size_members(F, members, hubs, web_pts, pitch ** 3 * corr, None, pitch)
    symmetric_radii(F, members, orb, axes, centre)
    push_inside(F, members, orb, axes, centre, om)
    round_joints(F, orb, hubs, axes, centre, grid)
    size_members(F, members, hubs, web_pts, pitch ** 3 * corr, depth, pitch)
    symmetric_radii(F, members, orb, axes, centre)
    R0 = nodes(F, members, depth, D)
    symmetric_radii(F, members, orb, axes, centre)
    # radius classes (length-weighted), then node radii from the class radii
    vals = [m[k] for m in members for k in ('r0', 'r1')]
    w = np.repeat([m['L'] for m in members], 2)
    cv = classes(vals, w, DT['rclass'], DT['digits'])
    for n, m in enumerate(members):
        m['r0'], m['r1'] = cv[2 * n], cv[2 * n + 1]
    req = {}
    for k in R0:
        ends = [m['r0'] for m in members if m['i'] == k] + [m['r1'] for m in members if m['j'] == k]
        req[k] = NODE_F * max(ends)
    ks = sorted(R0)
    raw = [max(R0[k], req[k]) for k in ks]
    cn = classes(raw, np.ones(len(raw)), DT['nclass'], DT['digits'])
    levels = sorted(set(cn))
    node_r = {}
    for k, v in zip(ks, cn):
        if v < req[k] - 1e-9:
            v = next((c for c in levels if c >= req[k] - 1e-9), sig(req[k] * 1.0001, DT['digits']))
            if v < req[k]:
                v = sig(req[k], DT['digits']) + 10 ** (math.floor(math.log10(req[k])) - DT['digits'] + 1)
        node_r[k] = v
    # renumber the joints: free joints and anchors used by members
    used = sorted({m['i'] for m in members} | {m['j'] for m in members})
    new = {k: n for n, k in enumerate(used)}
    joints = [tuple(float(x) for x in F.P[k]) for k in used]
    kind = ['hub' if F.hub[k] is not None else 'free' for k in used]
    hub_of = {new[k]: F.hub[k] for k in used if F.hub[k] is not None}
    orbit = [orb[k][0] for k in used]
    out_m = []
    share = {}
    for m in members:
        for k in (m['i'], m['j']):
            share[k] = share.get(k, 0) + 1
    for m in members:
        i, j = new[m['i']], new[m['j']]
        s = []
        for k, rk in ((m['i'], 'r0'), (m['j'], 'r1')):
            if F.hub[k] is not None:
                dep = hub_depth(hubs[F.hub[k]], F.P[k])
                sk = 0.0 if share[k] == 1 else sig(HUB_SHARE_F * dep, DT['digits'])
                if m[rk] + sk > dep:                  # the end section inside the hub
                    m[rk] = max(MIN_R, sig(dep - sk, DT['digits']) - 10 ** (-DT['digits'] + 1))
                s.append(sk)
            else:
                s.append(sig(TRIM_F * node_r[k], DT['digits']))
        out_m.append(dict(i=i, j=j, r0=float(m['r0']), r1=float(m['r1']), s0=float(s[0]), s1=float(s[1])))
    node_list = [(new[k], float(node_r[k])) for k in ks if k in new]
    R = dict(hubs=hubs, joints=joints, kind=kind, hub_of=hub_of, members=out_m, nodes=node_list, sym=axes,
             centre=tuple(float(c) for c in centre), orbit=orbit, ops=[orb[k][1] for k in used], detail=detail,
             pitch=pitch)
    log('frame: %d joints (%d on hubs), %d members, %d nodes, %d joint orbits' % (
        len(joints), len(hub_of), len(out_m), len(node_list), len(set(orbit))))
    # a truss has slender members; a thick wall (Knuckle_TO_2) gives members about as long as they are wide, and the
    # frame becomes a chain of nodes
    J = np.array(joints)
    slender = np.array([np.linalg.norm(J[m['j']] - J[m['i']]) / (m['r0'] + m['r1']) for m in out_m])
    R['stubby'] = float(np.mean(slender < STUBBY)) if len(slender) else 0.0
    if R['stubby'] > STUBBY_SHARE:
        log('WARNING: the TO is not truss-like: %.0f %% of the members are shorter than %.1f x the sum of their end '
            'radii; the frame model is a chain of nodes (a thick wall needs another method)' % (100 * R['stubby'], STUBBY))
    return R


# ---- CadQuery model --------------------------------------------------------------------------------------
EPS = 1e-3


def mass_props(shape):
    c = shape.Center()
    return float(shape.Volume()), (c.x, c.y, c.z)


def hub_solid(h, r, ext=0.0):
    k = h['axis']
    p = hubs3d.point(h, h['lo'] - ext)
    d = [0.0, 0.0, 0.0]; d[k] = 1.0
    return cq.Solid.makeCylinder(r, h['hi'] - h['lo'] + 2 * ext, Vector(*p), Vector(*d))


def member_geom(R, m):
    """(start centre, end centre, direction, r0, r1) of a member's frustum."""
    A, B = np.array(R['joints'][m['i']]), np.array(R['joints'][m['j']])
    t = (B - A) / np.linalg.norm(B - A)
    return A + m['s0'] * t, B - m['s1'] * t, t, m['r0'], m['r1']


def member_solid(R, m):
    p0, p1, t, r0, r1 = member_geom(R, m)
    h = float(np.linalg.norm(p1 - p0))
    if abs(r0 - r1) < 1e-9:
        return cq.Solid.makeCylinder(r0, h, Vector(*p0), Vector(*t))
    return cq.Solid.makeCone(r0, r1, h, Vector(*p0), Vector(*t))


NODE_AXIS = (0.2672612, 0.5345225, 0.8017837)     # sphere axis in the CadQuery model: a pole or the seam on an
                                                  # intersection curve made STEP read back the wrong face


def node_solid(R, j, r):
    return cq.Solid.makeSphere(r, Vector(*R['joints'][j]), Vector(*NODE_AXIS), angleDegrees1=-90, angleDegrees2=90)


def steps(R):
    """Features in build order (hubs, nodes, members, holes) with the volume and centre of mass after each."""
    out, union = [], None

    def combine(solid, cut=False):
        if union is None:
            return solid
        n0, v0, vs = len(union.Solids()), union.Volume(), solid.Volume()
        for tol in (None, 1e-4, 1e-3):
            try:
                u = (union.cut(solid, tol=tol) if cut else union.fuse(solid, tol=tol)).clean()
            except Exception:
                continue
            v = u.Volume()
            e = 1e-4 * (v0 + vs) + 1e-2
            ok_v = (v0 - vs - e <= v <= v0 + e) if cut else (max(v0, vs) - e <= v <= v0 + vs + e)
            if u.isValid() and len(u.Solids()) <= n0 + (0 if cut else 1) and ok_v:
                return u
        return None

    def add(kind, name, solid, cut=False, **kw):
        nonlocal union
        u = combine(solid, cut)
        if u is None:
            raise RuntimeError('%s: no valid Boolean' % name)
        union = u
        v, c = mass_props(union)
        out.append(dict(kind=kind, name=name, volume=v, com=np.array(c), own=float(solid.Volume()), **kw))

    for n, h in enumerate(R['hubs']):
        add('hub', 'Hub %d' % (n + 1), hub_solid(h, h['ro']), hub=n)
    for n, (j, r) in enumerate(R['nodes']):
        add('node', 'Node %d' % (n + 1), node_solid(R, j, r), joint=j, r=r)
    for n, m in enumerate(R['members']):
        add('member', 'Member %d' % (n + 1), member_solid(R, m), m=n)
    for n, h in enumerate(R['hubs']):
        add('hole', 'Hole %d' % (n + 1), hub_solid(h, h['ri'], EPS), cut=True, hub=n)
    return out, union
