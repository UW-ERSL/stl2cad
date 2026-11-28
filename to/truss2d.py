"""Planar trusses on the faces of a mirror-symmetric TO (cantilever Case 7: class 8/9, rods + plate).

The TO is a planar truss on each mirror face of the domain (mirror images about the mid-plane), planar structures
in face layers on other faces (Case 7: plate, bar and fork on y = 0), and webs between the mirror faces. Each is
fitted in its own plane by _planar; fit() finds the layers (see its comments). Method of the face trusses:
  1. Mirror axis k (heightfield.mirror_axis) and the break b between web and truss: the two-piece fit of the
     section area A(h) (heightfield.levels, K = 2). The truss is the TO with |h| >= b, the web the TO with |h| < b.
  2. Truss graph: projection of the truss half (h >= 0) onto the face, skeleton (2D), graph of pixel paths,
     spurs pruned (not those ending on the outline: supports), each path cut into straight members (Douglas-
     Peucker, tolerance by detail), close joints merged in clusters; end joints within reach of the outline put
     on it. Member lines fitted to their skeleton pixels; joints at the least-squares intersection of their
     members' lines, in the directions the lines determine (supports where their member's line meets the outline).
  3. Member widths: twice the median distance-transform value on the member's skeleton pixels away from its
     joints (the in-plane width of the TO member), in width classes (by detail). Member depth: one global value t,
     the truss volume of one face over the section's area (flat bars of equal volume; swtruss.py sets it from
     the section it builds).
  4. Chords: members along which the TO reaches the outline (CHORD_SHARE of the length); a chord is a strip
     flush with the outline line nearest its middle, its width the TO's depth from the outline; a chord end with
     no other member runs on along the outline as far as the TO does (to a support face, or the TO's own end);
     a chord end within one width of the outline runs on to it.
  5. Webs: the regions of the projection of the web's core (|h| < WEB_CORE x b, clear of the truss members) along
     the outline, each a strip flush with the outline (depth: mean section area over length), between the two
     truss faces.

fit(to, original, detail) -> the truss (joints, members with widths and chord edges, depth, webs); the models are
built from it by swtruss.py (member-driven CadQuery and SolidWorks scripts).
"""
import numpy as np
from scipy import ndimage
from shapely.geometry import Polygon, LineString, Point, MultiPoint
from shapely.ops import unary_union
from skimage.morphology import skeletonize, remove_small_holes
import heightfield as hf
import plates as _plates

PITCH = 0.25                                    # mm, projection pixels
DP_TOL = {1: 3.0, 2: 2.0, 3: 1.5, 4: 1.0, 5: 0.6}  # mm, straight members
WIDTH_CLASSES = {1: 1, 2: 2, 3: 2, 4: 3, 5: 0}   # 0: one width per member
SPUR = 3.0                                      # spurs shorter than SPUR x the local half-width are pruned
MERGE = 2.5                                     # joints closer than MERGE x the widest member half-width are merged
CHORD_GAP = 1.0                                 # mm: the TO reaches the outline where outline distance <= half-width + this
CHORD_SHARE = 0.6                               # a member is a chord when the TO reaches the outline along this share
LAYER = 3.0                                     # a face layer's section is this x larger than beyond it
LAYER_CORE = 0.1                                # ... and holds this share of its volume near the mid-plane
HOLE_MIN = 4.0                                  # mm^2: smaller holes in a projection are filled
WSTEP = 1.8                                     # a member is split where its width steps by this factor
LAYER_SPAN = 0.25                               # face layers are looked for within this share of the extent
WEB_CORE = 0.5                                  # webs are found in |h| < WEB_CORE x b
WEB_MIN = 0.002                                 # web regions below this share of the outline area are dropped
FILLET_SHARE = 0.4                              # corner radius: this x the narrowest member
COND = 0.1                                      # joint fit: a direction is fixed by the lines when its eigenvalue > this
EDGE_REACH = 1.0                                # a joint within EDGE_REACH x its half-width of the outline goes on it


# ---------------------------------------------------------------- skeleton graph
def _graph(S):
    """Pixel skeleton S (bool 2D) -> nodes {id: (i, j)}, edges [(n0, n1, [pixels])]."""
    nb = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    pts = set(map(tuple, np.argwhere(S)))
    deg = {p: sum((p[0] + a, p[1] + b) in pts for a, b in nb) for p in pts}
    node_px = {p for p in pts if deg[p] != 2}
    # clusters of adjacent node pixels are one node
    lab, nid = {}, 0
    for p in node_px:
        if p in lab:
            continue
        stack = [p]; lab[p] = nid
        while stack:
            q = stack.pop()
            for a, b in nb:
                r = (q[0] + a, q[1] + b)
                if r in node_px and r not in lab:
                    lab[r] = nid; stack.append(r)
        nid += 1
    nodes = {k: np.mean([p for p in node_px if lab[p] == k], axis=0) for k in range(nid)}
    edges, seen = [], set()
    for p in node_px:
        for a, b in nb:
            r = (p[0] + a, p[1] + b)
            if r not in pts or r in node_px or (p, r) in seen:
                continue
            path, prev, cur = [p, r], p, r
            while cur not in node_px:
                nxt = [(cur[0] + a2, cur[1] + b2) for a2, b2 in nb if (cur[0] + a2, cur[1] + b2) in pts
                       and (cur[0] + a2, cur[1] + b2) != prev and (cur[0] + a2, cur[1] + b2) not in path[-3:]]
                if not nxt:
                    break
                prev, cur = cur, nxt[0]; path.append(cur)
            if cur in node_px:
                seen.add((cur, path[-2]))
                edges.append((lab[p], lab[cur], path))
    return nodes, edges


def _dp(P, tol):
    """Douglas-Peucker indices of a polyline P (n x 2)."""
    keep = [0, len(P) - 1]

    def rec(i, j):
        if j <= i + 1:
            return
        a, b = P[i], P[j]; ab = b - a; L = np.linalg.norm(ab)
        Q = P[i + 1:j] - a
        d = np.abs(ab[0] * Q[:, 1] - ab[1] * Q[:, 0]) / L if L > 0 else np.linalg.norm(P[i + 1:j] - a, axis=1)
        k = int(np.argmax(d))
        if d[k] > tol:
            keep.append(i + 1 + k); rec(i, i + 1 + k); rec(i + 1 + k, j)
    rec(0, len(P) - 1)
    return sorted(keep)


def _width_steps(path, idx, DT):
    """Members split where the width steps (a plate running into a bar): at the pixel that best separates two
    pieces whose mean half-widths differ by more than WSTEP x, both pieces at least as long as the wider half-width."""
    out = [idx[0]]
    for s in range(len(idx) - 1):
        stack, cuts = [(idx[s], idx[s + 1])], []
        while stack:
            i0, i1 = stack.pop()
            d = np.array([DT[tuple(path[i])] for i in range(i0, i1 + 1)])
            best = None
            for j in range(2, len(d) - 2):
                m1, m2 = d[:j].mean(), d[j:].mean()
                hi = max(m1, m2)
                if hi > WSTEP * min(m1, m2) and min(j, len(d) - j) * PITCH >= hi:
                    sc = abs(m1 - m2) * min(j, len(d) - j)
                    if best is None or sc > best[0]:
                        best = (sc, i0 + j)
            if best is not None:
                cuts.append(best[1]); stack += [(i0, best[1]), (best[1], i1)]
        out += sorted(cuts) + [idx[s + 1]]
    return out


def _classes(w, n):
    """Width classes: 1D k-means of the widths into n values (each member gets its class value)."""
    w = np.asarray(w, float)
    if n == 0 or n >= len(set(np.round(w, 2))):
        return w.copy()
    c = np.quantile(w, np.linspace(0, 1, n + 2)[1:-1])
    for _ in range(50):
        a = np.argmin(np.abs(w[:, None] - c[None]), 1)
        c = np.array([w[a == k].mean() if (a == k).any() else c[k] for k in range(n)])
    return c[np.argmin(np.abs(w[:, None] - c[None]), 1)]


# ---------------------------------------------------------------- truss fit
def _outline(original, axis, level, cut=()):
    """Section of the original normal to axis at level, rounded to 1e-3 mm (the STL's float noise off the exact
    outline), minus the bands cut = [(axis index in the section frame, lo, hi)] (face layers fitted apart)."""
    O = _plates.section(original, axis, level)
    O = O if O.geom_type == 'Polygon' else max(O.geoms, key=lambda q: q.area)
    O = Polygon(np.round(np.array(O.exterior.coords), 3), [np.round(np.array(i.coords), 3) for i in O.interiors])
    for i, lo, hi in cut:
        x0, y0, x1, y1 = O.bounds
        band = Polygon([(lo, y0 - 1), (hi, y0 - 1), (hi, y1 + 1), (lo, y1 + 1)]) if i == 0 else \
            Polygon([(x0 - 1, lo), (x1 + 1, lo), (x1 + 1, hi), (x0 - 1, hi)])
        O = O.difference(band)
        O = O if O.geom_type == 'Polygon' else max(O.geoms, key=lambda q: q.area)
        O = Polygon(np.round(np.array(O.exterior.coords), 3))
    return O


def _break(A):
    """Two-piece constant fit of A (least squares): end index of the first piece and the two means."""
    e = hf._partition(np.asarray(A, float), 2)[0]
    return e, float(np.mean(A[:e])), float(np.mean(A[e:]))


def fit(to, original, detail=4, log=print):
    """Structures of a truss-like TO with a mirror plane: the planar trusses on the two mirror faces (pair), the
    face layers on other faces (each a planar truss in its own plane), and the webs between the mirror faces.
    Returns dict(axis, mid, H, structures=[...], webs, ...)."""
    ob = original.bounds
    k = hf.mirror_axis(to, ob)
    mid = 0.5 * (ob[0][k] + ob[1][k]); H = 0.5 * (ob[1][k] - ob[0][k])
    g = to.voxelized(PITCH).fill(); F = np.asarray(g.matrix); o = np.asarray(g.transform)[:3, 3]
    vs = to.volume / F.sum()                            # voxel volume scaled to the mesh's (the fill counts
    ax = [o[a] + PITCH * np.arange(F.shape[a]) for a in range(3)]   # boundary voxels whole: +12 % on 4 mm bars)
    # face layers on the faces other than the mirror faces: the section area parallel to the face, A(s) at a depth
    # s from it, is LAYER x larger in a first piece than beyond it
    layers, L = [], np.zeros(F.shape, bool)
    for a in range(3):
        if a == k:
            continue
        ext = ob[1][a] - ob[0][a]
        A = F.sum(axis=tuple(i for i in range(3) if i != a)).astype(float)
        for side in (0, 1):
            face = ob[side][a]
            dep = (ax[a] - face) * (1 if side == 0 else -1)
            sel = np.where((dep >= -PITCH) & (dep <= LAYER_SPAN * ext))[0]
            sel = sel[np.argsort(dep[sel])]
            if len(sel) < 8:
                continue
            e, m1, m2 = _break(A[sel])
            if m1 > LAYER * m2 and e < len(sel) - 2:
                d = float(dep[sel][e - 1] + PITCH / 2)
                band = [slice(None)] * 3
                band[a] = (dep >= -PITCH) & (dep < d)
                hk = (ax[k] - mid).reshape([-1 if i == k else 1 for i in range(3)])
                core = (F[tuple(band)] & (np.abs(hk) < H / 2)).sum() / max(F[tuple(band)].sum(), 1)
                if core < LAYER_CORE:       # material only near the mirror faces: chords of the face trusses
                    log('no face layer on %s = %.2f: %.0f %% of the band near the mid-plane' % ('xyz'[a], face, 100 * core))
                    continue
                layers.append(dict(axis=a, side=side, face=float(face), depth=d))
                sl = [slice(None)] * 3
                sl[a] = (dep >= -PITCH) & (dep < d)
                L[tuple(sl)] = True
                log('face layer on %s = %.2f: %.2f mm deep (section area %.0f, beyond it %.0f mm^2)'
                    % ('xyz'[a], face, d, m1 * PITCH ** 2, m2 * PITCH ** 2))
    # mirror faces: web / truss break b from the two-piece fit of A(h), material outside the face layers
    u, v = _plates.frame(k)
    G = np.moveaxis(F & ~L, (u, v, k), (0, 1, 2)); o2 = o[[u, v]]
    h = ax[k] - mid
    Ah = np.array([G[:, :, np.abs(np.abs(h) - x) < PITCH / 2].sum() for x in np.arange(PITCH / 2, H, PITCH)], float)
    e, _, _ = _break(Ah)
    b = (e + 0.5) * PITCH
    lstructs = []
    for ly in layers:                                    # each face layer: a planar truss in its own plane
        a = ly['axis']; uu, vv = _plates.frame(a)
        dep = (ax[a] - ly['face']) * (1 if ly['side'] == 0 else -1)
        Gl = np.moveaxis(F & L & (dep < ly['depth']).reshape([-1 if i == a else 1 for i in range(3)]),
                         (uu, vv, a), (0, 1, 2))
        Ol = _outline(original, a, ly['face'] + (0.5 if ly['side'] == 0 else -0.5) * ly['depth'])
        Sl = _planar(Gl.any(2), o[[uu, vv]], Ol, detail, log, 'layer %s = %.2f' % ('xyz'[a], ly['face']))
        vl = Gl.sum() * vs
        Sl.update(axis=a, frame=(uu, vv), name='layer', v=vl, t=vl / Sl['poly'].area,
                  faces=[(ly['face'], 1 if ly['side'] == 0 else -1)], O=Ol, layer=ly)
        ly['t'] = Sl['t']
        lstructs.append(Sl)
    # the face trusses stand on the layers: their outline starts at each layer's thickness t
    cut = [(([i for i in (u, v)].index(ly['axis'])), *sorted((ly['face'], ly['face'] + (ly['t'] if ly['side'] == 0
            else -ly['t'])))) for ly in layers]
    O = _outline(original, k, mid, cut)
    proj = G[:, :, h >= b].any(2)                       # truss of the + face, projected
    webm = G[:, :, np.abs(h) < WEB_CORE * b]            # core of the web slab (clear of the truss members)
    v_truss = 0.5 * G[:, :, np.abs(h) >= b].sum() * vs
    v_web = G[:, :, np.abs(h) < b].sum() * vs
    log('truss: mirror axis %s, mid-plane %.2f, web |h| < %.2f mm, truss volume %.0f mm^3 per face, web %.0f mm^3'
        % ('xyz'[k], mid, b, v_truss, v_web))
    S = _planar(proj, o2, O, detail, log, 'faces')
    S.update(axis=k, frame=(u, v), name='faces', v=v_truss, t=v_truss / S['poly'].area,
             faces=[(mid - H, 1), (mid + H, -1)], O=O)
    structures = [S] + lstructs
    webs = web_strips(webm, o2, O, v_web, vs, log)
    return dict(axis=k, mid=mid, H=H, b=b, structures=structures, webs=webs, v_web=v_web, layers=layers,
                J=S['J'], members=S['members'], t=S['t'], O=O, poly=S['poly'], v_truss=v_truss, on_edge=S['on_edge'])


def _planar(proj, o2, O, detail, log, label):
    """Planar truss from the projection proj (pixels, origin o2) inside the outline O (a polygon in the same
    plane): joints, members (with chords and run-on), on_edge, the idealized section poly."""
    o = np.r_[o2, 0.0]
    ring = O.exterior
    # skeleton graph in mm
    proj = ndimage.binary_opening(ndimage.binary_closing(proj, iterations=2), iterations=1)
    proj = remove_small_holes(proj, max_size=int(HOLE_MIN / PITCH ** 2))   # pinholes of the TO, not openings
    DT = ndimage.distance_transform_edt(proj) * PITCH
    S = skeletonize(proj)
    nodes, edges = _graph(S)
    X = lambda ij: o[:2] + PITCH * np.asarray(ij, float)

    def near_edge(n):                   # a skeleton end on the outline is a support (anchor), not a spur
        p = nodes[n]
        return ring.distance(Point(X(p))) < DT[tuple(np.round(p).astype(int))] + 2 * PITCH
    # prune spurs (an edge to a degree-1 node shorter than SPUR x its half-width), repeatedly
    while True:
        deg = {n: 0 for n in nodes}
        for a, c, _ in edges:
            deg[a] += 1; deg[c] += 1
        drop = [e for e in edges if e[0] != e[1] and any(deg[n] == 1 and not near_edge(n) for n in e[:2])
                and len(e[2]) * PITCH < SPUR * max(DT[tuple(p)] for p in e[2])]
        if not drop:
            break
        edges = [e for e in edges if e not in drop[:1]]
        edges = _join_deg2(nodes, edges)
    # straight members
    J = {n: X(nodes[n]) for n in nodes}
    members = []
    for a, c, path in edges:
        P = X(np.array(path)); idx = _dp(P, DP_TOL[detail])
        idx = _width_steps(path, idx, DT)
        ids = [a] + ['%d_%d_%d' % (a, c, i) for i in idx[1:-1]] + [c]
        for q, i in zip(ids[1:-1], idx[1:-1]):
            J[q] = P[i]
        for s in range(len(idx) - 1):
            seg = path[idx[s]:idx[s + 1] + 1]
            n0 = max(1, len(seg) // 5)
            core = seg[n0:len(seg) - n0] or seg
            w = 2 * float(np.median([DT[tuple(p)] for p in core]))
            members.append([ids[s], ids[s + 1], w, X(np.array(core)), np.array([DT[tuple(p)] for p in core])])
    # merge close joints: clusters of joints closer than MERGE x the narrower half-width, each at its joint of the
    # highest degree (a skeleton junction rather than a bend), in a fixed order
    used = sorted({m[0] for m in members} | {m[1] for m in members}, key=str)
    J = {q: J[q] for q in used}
    wmax = {q: max(m[2] for m in members if q in m[:2]) for q in J}
    deg = {q: sum(q in m[:2] for m in members) for q in J}
    root = {q: q for q in used}
    find = lambda q: q if root[q] == q else find(root[q])
    for i, p in enumerate(used):
        for q in used[i + 1:]:
            if np.linalg.norm(J[p] - J[q]) < MERGE * 0.5 * min(wmax[p], wmax[q]):
                root[find(q)] = find(p)
    for r in sorted({find(q) for q in used}, key=str):
        cl = [q for q in used if find(q) == r]
        top = max(deg[q] for q in cl)
        x = np.mean([J[q] for q in cl if deg[q] == top], axis=0)
        keep = cl[0]
        for q in cl:
            if q != keep:
                del J[q]
        J[keep] = x; wmax[keep] = max(wmax[q] for q in cl)
        for m in members:
            m[0] = keep if m[0] in cl else m[0]; m[1] = keep if m[1] in cl else m[1]
    members = [m for m in members if m[0] != m[1]]
    # end joints (one member) near the outline are supports: on it
    on_edge = set()
    deg = {q: sum(q in m[:2] for m in members) for q in J}
    for q in J:
        if deg[q] == 1 and (ring.distance(Point(J[q])) < EDGE_REACH * 0.5 * wmax[q] + 0.5 or not O.contains(Point(J[q]))):
            J[q] = np.array(ring.interpolate(ring.project(Point(J[q]))).coords[0]); on_edge.add(q)
    # member lines: least-squares line through each member's skeleton pixels
    lines = {}
    for m in members:
        c0 = m[3].mean(0); _, _, Vt = np.linalg.svd(m[3] - c0)
        lines[id(m)] = (c0, Vt[0])
    # chords: members along which the TO reaches the outline; their width is the TO's depth from the outline
    for m in members:
        core, dt = m[3], m[4]
        reach = np.array([ring.distance(Point(p)) for p in core])
        m[3] = None
        if (reach <= dt + CHORD_GAP).mean() >= CHORD_SHARE:
            f, n = _inward(O, core[len(core) // 2])
            m[2] = float(np.median(reach + dt)); m[3] = (f, n)
        del m[4]
    # consecutive chords on one outline line: one chord, one width (length-weighted mean)
    ch = [i for i, m in enumerate(members) if m[3] is not None]
    grp = {i: i for i in ch}
    find = lambda i: i if grp[i] == i else find(grp[i])
    for i in ch:
        for j in ch:
            if i < j and set(members[i][:2]) & set(members[j][:2]):
                (fi, ni), (fj, nj) = members[i][3], members[j][3]
                if ni @ nj > 0.999 and abs((fj - fi) @ ni) < 0.1:
                    grp[find(j)] = find(i)
    for r in {find(i) for i in ch}:
        g = [i for i in ch if find(i) == r]
        ln = [np.linalg.norm(J[members[i][0]] - J[members[i][1]]) for i in g]
        w = float(np.average([members[i][2] for i in g], weights=ln))
        for i in g:
            members[i][2] = w; members[i][3] = members[r][3]
    W = _classes([m[2] for m in members], WIDTH_CLASSES[detail])
    for m, w in zip(members, W):
        m[2] = float(w)
    for m in members:                         # a chord's line: its centre line, half its width from the outline
        if m[3] is not None:
            f, n = m[3]
            lines[id(m)] = (f + n * m[2] / 2, np.array([n[1], -n[0]]))
    # joints: least-squares intersection of their members' lines, in the directions the lines determine; an end
    # joint slides onto its member's line
    for q in J:
        inc = [lines[id(m)] for m in members if q in m[:2]]
        A = np.zeros((2, 2)); r = np.zeros(2)
        for c0, d in inc:
            Pn = np.eye(2) - np.outer(d, d); A += Pn; r += Pn @ c0
        lam, E = np.linalg.eigh(A)
        for li, e in zip(lam, E.T):           # only along directions the lines fix (nearly parallel lines fix
            if li > COND:                     # nothing along them: the skeleton position is kept there)
                J[q] = J[q] + e * ((e @ r) - li * (e @ J[q])) / li
    for q in on_edge:                         # supports: where the member's line meets the outline
        m = next(m for m in members if q in m[:2])
        if m[3] is not None:
            continue
        c0, d = lines[id(m)]
        hit = LineString([c0 - 1e3 * d, c0 + 1e3 * d]).intersection(ring)
        pts = [np.array(g.coords[0]) for g in getattr(hit, 'geoms', [hit]) if not g.is_empty]
        if pts:
            J[q] = min(pts, key=lambda x: np.linalg.norm(x - J[q]))
    # free chord ends (one member): the strip runs on along the outline as far as the TO does (walked on the
    # projection at half the chord's depth)
    deg = {q: sum(q in m[:2] for m in members) for q in J}
    inproj = lambda x: (lambda ij: 0 <= ij[0] < proj.shape[0] and 0 <= ij[1] < proj.shape[1] and proj[ij])(
        tuple(np.round((x - o[:2]) / PITCH).astype(int)))
    for m in members:
        ext = [0.0, 0.0]
        if m[3] is not None:
            f, n = m[3]; tt = np.array([n[1], -n[0]])
            sa, sc = (J[m[0]] - f) @ tt, (J[m[1]] - f) @ tt
            u = 1.0 if sc >= sa else -1.0
            for i, (q, s0, dirn) in enumerate(((m[0], sa, -u), (m[1], sc, u))):
                if deg[q] != 1:
                    continue
                e = 0.0
                while inproj(f + (s0 + dirn * (e + PITCH)) * tt + n * m[2] / 2) and e < 1e3:
                    e += PITCH
                ext[i] = e
        if m[3] is not None:                  # an end within one width of the outline runs on to it
            for i, (q, s0, dirn) in enumerate(((m[0], sa, -u), (m[1], sc, u))):
                x = f + (s0 + dirn * ext[i]) * tt + n * m[2] / 2
                hit = LineString([x, x + dirn * tt * 1e3]).intersection(ring)
                pts = [np.array(g.coords[0]) for g in getattr(hit, 'geoms', [hit]) if not g.is_empty]
                if pts:
                    g = min(np.linalg.norm(p - x) for p in pts)
                    if g < max(m[2], 1.0):
                        ext[i] += g + 0.5
        m.append(tuple(ext))
    poly = section(J, members, O)
    log('%s: %d joints (%d on the outline), %d members, widths %s mm'
        % (label, len(J), len(on_edge), len(members), ', '.join('%.2f' % w for w in sorted(set(W)))))
    return dict(J=J, members=members, on_edge=on_edge, poly=poly, widths=sorted(set(W)))


def _join_deg2(nodes, edges):
    """Nodes left with two edges after pruning are dissolved (their edges joined)."""
    while True:
        deg = {}
        for a, c, _ in edges:
            deg[a] = deg.get(a, 0) + 1; deg[c] = deg.get(c, 0) + 1
        n2 = [n for n, d in deg.items() if d == 2 and not any(e[0] == e[1] == n for e in edges)]
        if not n2:
            return edges
        n = n2[0]
        e1, e2 = [e for e in edges if n in e[:2]]
        p1 = e1[2] if e1[1] == n else e1[2][::-1]
        p2 = e2[2] if e2[0] == n else e2[2][::-1]
        a = e1[0] if e1[1] == n else e1[1]; c = e2[1] if e2[0] == n else e2[0]
        edges = [e for e in edges if e is not e1 and e is not e2] + [(a, c, p1 + p2[1:])]


def _inward(O, x):
    """Foot of x on the outline of O and the unit inward normal there."""
    ring = O.exterior; s = ring.project(Point(x))
    f = np.array(ring.interpolate(s).coords[0])
    t = np.array(ring.interpolate(s + 0.01).coords[0]) - np.array(ring.interpolate(s - 0.01).coords[0])
    n = np.array([-t[1], t[0]]) / (np.linalg.norm(t) + 1e-12)
    if not O.contains(Point(f + 0.05 * n)):
        n = -n
    return f, n


def chord(A, C, w, edge, ext=(0.0, 0.0)):
    """Chord strip: from the foot of A to the foot of C on the outline line edge = (point, inward normal), w deep,
    overlapping its neighbours by 0.02 mm; ext: lengths by which the ends run on along the outline."""
    f, n = edge; t = np.array([n[1], -n[0]])
    sa, sc = (A - f) @ t, (C - f) @ t
    u = 1.0 if sc >= sa else -1.0
    sa, sc = sa - u * (ext[0] + 0.02), sc + u * (ext[1] + 0.02)
    fa, fc = f + sa * t, f + sc * t
    return Polygon([fa - 0.05 * n, fc - 0.05 * n, fc + w * n, fa + w * n]), (fa, fa + w * n), (fc, fc + w * n)


def section(J, members, O, log=None):
    """Idealized truss section: member rectangles between joints (chords flush with the outline), joint gussets,
    clipped to the outline O."""
    parts = []
    ends = {q: [] for q in J}
    for m in members:
        a, c, w, edge, ext = m
        A, C = J[a], J[c]; d = C - A; L = np.linalg.norm(d)
        if L < 1e-6:
            continue
        d /= L; n = np.array([-d[1], d[0]]) * w / 2
        if edge is not None:                  # chord: its strip; its end corners join the gussets
            strip, ea, ec = chord(A, C, w, edge, ext)
            parts.append(strip); ends[a] += list(ea); ends[c] += list(ec)
            continue
        parts.append(Polygon([A + n, C + n, C - n, A - n]))
        s = min(w, 0.45 * L)
        ends[a] += [A + d * s + n, A + d * s - n]
        ends[c] += [C - d * s + n, C - d * s - n]
    for q, pts in ends.items():
        if len(pts) >= 4:
            parts.append(MultiPoint([tuple(p) for p in pts] + [tuple(J[q])]).convex_hull)
    P = unary_union(parts).intersection(Polygon(O.exterior))
    return P if P.geom_type == 'Polygon' else max(P.geoms, key=lambda g: g.area)


def web_strips(webm, o, O, v_web, vs, log=print):
    """Webs: regions of the projection of the web's core (webm: its voxel columns); a region along the outline
    becomes a strip flush with it over the region's extent, its depth the region's mean section area over its
    length. Returns [(strip polygon, share of the core volume)]."""
    P2 = webm.any(2)
    lab, n = ndimage.label(P2, structure=np.ones((3, 3)))
    vol = ndimage.sum(webm.sum(2), lab, range(1, n + 1))
    out = []
    for i in range(1, n + 1):
        ij = np.argwhere(lab == i); X = o[:2] + PITCH * ij
        if len(ij) * PITCH ** 2 < WEB_MIN * O.area:
            continue
        reach = np.array([O.exterior.distance(Point(x)) for x in X[::max(1, len(X) // 400)]])
        if reach.min() > 2 * PITCH:
            log('  web region of %.0f mm^2 off the outline dropped' % (len(ij) * PITCH ** 2)); continue
        f, nrm = _inward(O, X[np.argmin([O.exterior.distance(Point(x)) for x in X[::max(1, len(X) // 400)]]) * max(1, len(X) // 400)])
        tt = np.array([nrm[1], -nrm[0]]); sv = (X - f) @ tt
        s0, s1 = sv.min() - PITCH / 2, sv.max() + PITCH / 2
        d = vol[i - 1] * vs / PITCH / webm.shape[2] / (s1 - s0)   # mean section area over the core, per length
        strip = Polygon([f + s0 * tt - 0.05 * nrm, f + s1 * tt - 0.05 * nrm, f + s1 * tt + d * nrm, f + s0 * tt + d * nrm])
        out.append((strip.intersection(Polygon(O.exterior)), vol[i - 1] / webm.sum()))
    return out
