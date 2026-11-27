"""Trusses of straight round struts (lattice arms, space frames, organic trusses).

Recognition from the mesh:
  1. local radius and centre-line points: rays from surface samples along the inward normal cross a round strut
     through its axis, so half the chord is the local radius and the sample moved inward by it is on the axis;
  2. struts: straight runs of centre points of one radius (RANSAC lines, grown over the points of that radius,
     split where the points stop for more than two radii); a run must be at least five radii long, its radius
     steady (spread below 8 %) and its surface well sampled;
  3. radii: measured per strut by rays from its axis (stations widened by a crossing strut or a joint dropped),
     rounded within the tolerance, not snapped to classes; a strut whose radius changes along it by more than
     10 % keeps a start and an end radius (a cone);
  4. nodes: strut ends whose axes meet (axis distance below half the smaller radius, ends within a few radii of
     the meeting point) form a node at the least-squares point of their axes; struts are cut or extended to their
     nodes; an end on another strut's axis (a T-joint) gets a node there; collinear struts through a node stay
     separate members (each node to node);
  5. the rest of the part (bushings, rings, pads: surface far from every strut) is fitted separately (fit_rest).
"""

import math
import time

import numpy as np

RUN_GAP = 2.0          # x radius: a gap in the centre points that ends a strut
MIN_LEN = 5.0          # x radius: shortest strut
MAX_SPREAD = 0.08      # radius spread within a strut
MIN_FILL = 0.35        # sampled share of a strut's surface
CLASS_GAP = 0.10       # relative radius gap between classes
JOINT_SCALE = 0.85    # joint ball radius / the largest ball inside the part at the joint (user: a little smaller)
FILLET_PIECES = 6     # joint-blend pieces segmented for the fillet radius (reported only)
TAPER = 0.10           # relative radius change along a strut that makes it a cone


def _intersector(mesh):
    """Embree ray intersector when embreex is installed (about 1000 times faster here, same distances to
    1e-13 mm on truss-connector), else trimesh's own."""
    try:
        from trimesh.ray.ray_pyembree import RayMeshIntersector
        return RayMeshIntersector(mesh)
    except Exception:
        from trimesh.ray.ray_triangle import RayMeshIntersector
        return RayMeshIntersector(mesh)


def _rays(mesh, n, seed=0):
    """Surface samples, normals and the chord length along the inward normal (inf where none)."""
    import trimesh
    P, fi = trimesh.sample.sample_surface(mesh, n, seed=seed)
    N = mesh.face_normals[fi]
    O = P - 1e-6 * float(np.linalg.norm(mesh.extents)) * N
    d = np.full(len(P), np.inf)
    rm = _intersector(mesh)
    for ch in np.array_split(np.arange(len(P)), max(1, len(P) // 1000)):
        loc, ri, _ = rm.intersects_location(O[ch], -N[ch], multiple_hits=False)
        d[ch[ri]] = np.linalg.norm(loc - O[ch][ri], axis=1)
    return P, N, d


def find_struts(C, rr, density, seed=0, log=print):
    """Straight runs of centre points C (local radius rr). density: surface samples per unit area."""
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    T = cKDTree(C)
    rem = np.ones(len(C), bool)
    segs = []
    fails = 0
    while fails < 400 and rem.sum() > 50:
        idx = np.flatnonzero(rem)
        i = int(rng.choice(idx))
        r_i = rr[i]
        nb = [j for j in T.query_ball_point(C[i], 8 * r_i) if rem[j] and abs(rr[j] - r_i) < 0.15 * r_i]
        if len(nb) < 15:
            fails += 1
            continue
        best = None
        for _ in range(30):
            j = nb[rng.integers(len(nb))]
            v = C[j] - C[i]
            L = np.linalg.norm(v)
            if L < 3 * r_i:
                continue
            v = v / L
            Q = C[nb] - C[i]
            n_in = (np.linalg.norm(Q - np.outer(Q @ v, v), axis=1) < 0.2 * r_i).sum()
            if best is None or n_in > best[0]:
                best = (n_in, v)
        if best is None or best[0] < 12:
            fails += 1
            continue
        v, mu, run = best[1], C[i], None
        for _ in range(3):
            Q = C[idx] - mu
            t = Q @ v
            dist = np.linalg.norm(Q - np.outer(t, v), axis=1)
            inl = idx[(dist < 0.25 * r_i) & (np.abs(rr[idx] - r_i) < 0.2 * r_i)]
            tt = (C[inl] - mu) @ v
            o = np.argsort(tt)
            ts, ins = tt[o], inl[o]
            if len(ts) == 0:
                run = None
                break
            k0 = min(int(np.searchsorted(ts, 0.0)), len(ts) - 1)
            lo = hi = k0
            gap = RUN_GAP * r_i
            while lo > 0 and ts[lo] - ts[lo - 1] < gap:
                lo -= 1
            while hi < len(ts) - 1 and ts[hi + 1] - ts[hi] < gap:
                hi += 1
            run = ins[lo:hi + 1]
            if len(run) < 10:
                break
            mu = C[run].mean(0)
            v = np.linalg.svd(C[run] - mu)[2][0]
            r_i = float(np.median(rr[run]))
        if run is None or len(run) < 25:
            fails += 1
            continue
        t = (C[run] - mu) @ v
        L = float(t.max() - t.min())
        rem[run] = False
        spread = float(np.std(rr[run])) / r_i
        fill = len(run) / (density * 2 * math.pi * r_i * max(L, 1e-9))
        if L < MIN_LEN * r_i or spread > MAX_SPREAD or fill < MIN_FILL:
            fails += 1
            continue
        fails = 0
        # radius along the strut: start and end values of a linear fit
        A = np.c_[np.ones(len(t)), t - t.min()]
        a0, a1 = np.linalg.lstsq(A, rr[run], rcond=None)[0]
        segs.append(dict(p0=mu + t.min() * v, p1=mu + t.max() * v, r=r_i, r0=float(a0), r1=float(a0 + a1 * L),
                         n=len(run)))
    return segs


def radius_classes(values, tol):
    """Groups of radii (sorted, split at relative gaps above CLASS_GAP); each value -> its class radius, the
    class median rounded to the coarsest digits within tol / 2."""
    from .wires import _snap
    v = np.sort(np.asarray(values, float))
    groups, cur = [], [v[0]]
    for x in v[1:]:
        if x - cur[-1] > CLASS_GAP * cur[-1]:
            groups.append(cur)
            cur = [x]
        else:
            cur.append(x)
    groups.append(cur)
    classes = [_snap(float(np.median(g)), 0.1 * tol) for g in groups]
    bounds = [(g[0], g[-1]) for g in groups]

    def of(x):
        k = min(range(len(bounds)), key=lambda i: 0 if bounds[i][0] <= x <= bounds[i][1]
                else min(abs(x - bounds[i][0]), abs(x - bounds[i][1])))
        return classes[k]
    return classes, of


def _line_dist(p, u, q, w):
    n = np.cross(u, w)
    if np.linalg.norm(n) < 1e-9:
        return float(np.linalg.norm(np.cross(q - p, u)))
    return float(abs((q - p) @ n) / np.linalg.norm(n))


def _ls_point(ends):
    A, b = np.zeros((3, 3)), np.zeros(3)
    for p, u in ends:
        M = np.eye(3) - np.outer(u, u)
        A += M
        b += M @ p
    return np.linalg.lstsq(A, b, rcond=None)[0]


def find_nodes(segs, log=print):
    """Nodes where strut axes meet; each strut's ends are moved to its nodes. Returns nodes (list of dict p,
    members) and the free ends (strut, end)."""
    E = []
    for k, s in enumerate(segs):
        v = s["p1"] - s["p0"]
        v = v / np.linalg.norm(v)
        E.append(dict(k=k, e=0, p=s["p0"], u=-v, r=s["r"]))
        E.append(dict(k=k, e=1, p=s["p1"], u=v, r=s["r"]))
    n = len(E)

    def reach(a, b):                              # an end lies within this of its node (blend + other radius)
        return 3 * (a["r"] + b["r"]) + 2.0

    def cluster(ids):
        """Groups of the ends ids whose axes meet ahead of the ends; an end that does not fit its group's
        common point is dropped (one at a time) and returned as left over."""
        parent = {i: i for i in ids}

        def root(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for x_, i in enumerate(ids):
            for j in ids[x_ + 1:]:
                a, b = E[i], E[j]
                R_ = reach(a, b)
                if a["k"] == b["k"] or np.linalg.norm(a["p"] - b["p"]) > 2 * R_:
                    continue
                if _line_dist(a["p"], a["u"], b["p"], b["u"]) > 0.5 * min(a["r"], b["r"]):
                    continue
                x = _ls_point([(a["p"], a["u"]), (b["p"], b["u"])])
                ea, eb = (x - a["p"]) @ a["u"], (x - b["p"]) @ b["u"]
                if -a["r"] < ea < R_ and -b["r"] < eb < R_:
                    parent[root(i)] = root(j)
        groups = {}
        for i in ids:
            groups.setdefault(root(i), []).append(i)
        found, left = [], []
        for g in groups.values():
            g = list(g)
            while len(g) >= 2:
                x = _ls_point([(E[i]["p"], E[i]["u"]) for i in g])
                res = [np.linalg.norm(np.cross(x - E[i]["p"], E[i]["u"])) / E[i]["r"] for i in g]
                w = int(np.argmax(res))
                if res[w] < 0.5:
                    break
                left.append(g.pop(w))
            if len(g) < 2:
                left += g
            else:
                found.append(dict(p=x, ends=[(E[i]["k"], E[i]["e"]) for i in g]))
        return found, left
    nodes, free = cluster(list(range(n)))
    for _ in range(2):                            # ends dropped from a group may form a node of their own
        more, free = cluster(free)
        nodes += more
    # T-joints: a free end ahead of which lies another strut's axis
    still = []
    for i in free:
        a = E[i]
        best = None
        for k, s in enumerate(segs):
            if k == a["k"]:
                continue
            v = s["p1"] - s["p0"]
            L = np.linalg.norm(v)
            v = v / L
            # closest point on s's axis to a's axis
            w0 = a["p"] - s["p0"]
            b_ = a["u"] @ v
            den = 1 - b_ * b_
            if den < 1e-6:
                continue
            ta = (b_ * (w0 @ v) - (w0 @ a["u"])) / den
            ts = (w0 @ v + b_ * ta)
            q = s["p0"] + ts * v
            pa = a["p"] + ta * a["u"]
            if -a["r"] < ta < reach(a, dict(r=s["r"])) and -s["r"] < ts < L + s["r"] and \
                    np.linalg.norm(q - pa) < 0.5 * a["r"]:
                if best is None or ta < best[0]:
                    best = (ta, q, k)
        if best is None:
            still.append((a["k"], a["e"]))
        else:
            nodes.append(dict(p=best[1], ends=[(a["k"], a["e"])], on=best[2]))
    # move strut ends to their nodes
    for nd in nodes:
        for k, e in nd["ends"]:
            segs[k]["p1" if e else "p0"] = nd["p"].copy()
    return nodes, still


def detect(mesh, tol, log=print, n_samples=60000, rays=None):
    """Struts, nodes and the rest of the part (tubes, knobs) with the joint fillet radius, or None."""
    t0 = time.time()
    P, N, d = rays if rays is not None else _rays(mesh, n_samples)
    n_samples = len(P)
    ok = np.isfinite(d) & (d < 0.08 * float(np.linalg.norm(mesh.extents)))
    C = P[ok] - 0.5 * d[ok][:, None] * N[ok]
    rr = 0.5 * d[ok]
    segs = find_struts(C, rr, n_samples / float(mesh.area), log=log)
    if len(segs) < 3:
        return None
    log(f"  struts: {len(segs)} straight struts from the centre lines ({time.time() - t0:.1f}s)")
    rmax = float(np.percentile([s["r"] for s in segs], 90))
    tubes, spheres, extra, rf = fit_rest(mesh, segs, tol, rmax, log=log)
    segs = segs + extra
    nodes, free = find_nodes(segs, log=lambda *a: None)       # ends at the joints, for the radius stations
    measure_radii(mesh, segs)
    # radii as measured (rounded within tol / 10), not snapped to classes (user decision)
    from .wires import _snap
    for s in segs:
        s["taper"] = abs(s["r1"] - s["r0"]) > TAPER * s["r"]
        s["r"] = _snap(s["r"], 0.1 * tol)
        if s["taper"]:
            s["r0"], s["r1"] = _snap(s["r0"], 0.1 * tol), _snap(s["r1"], 0.1 * tol)
            s["taper"] = s["r0"] != s["r1"]
    classes, _ = radius_classes([s["r"] for s in segs], tol)      # reported only
    log(f"  struts: {len(segs)} struts, radii {min(s['r'] for s in segs):g}..{max(s['r'] for s in segs):g} "
        f"({sum(s['taper'] for s in segs)} tapered)")
    log(f"  struts: {len(nodes)} nodes ({sum(1 for nd in nodes if 'on' in nd)} on a strut's side), "
        f"{len(free)} free ends")
    spheres = knobs(mesh, segs, free, spheres, rmax, log=log)
    free = attach(segs, free, tubes, spheres, nodes, log=log)
    trim_free_ends(mesh, segs, free, log=log)
    segs = prune_outside(mesh, segs, log=log)
    segs = extend_free_ends(mesh, segs, spheres, tubes, log=log)
    joints = joint_spheres(mesh, segs, spheres, tubes, log=log)
    return dict(segs=segs, nodes=nodes, free=free, classes=classes, tubes=tubes, spheres=spheres, fillet=rf,
                joints=joints)

def _seg_surface_distance(P, segs):
    D = np.full(len(P), np.inf)
    for s in segs:
        a, b = s["p0"], s["p1"]
        v = b - a
        L = np.linalg.norm(v)
        v = v / L
        t = np.clip((P - a) @ v, 0, L)
        q = a + np.outer(t, v)
        D = np.minimum(D, np.abs(np.linalg.norm(P - q, axis=1) - s["r"]))
    return D


def _segment_piece(args):
    from .segment import Segmenter
    c, tol = args
    return Segmenter(c, tol, kinds=("plane", "cylinder", "sphere", "cone", "torus"), log=lambda *a: None).run()


def _segment_pieces(pieces, tol):
    """Segmenter on each piece, in worker processes when more than one CPU is available (same result as
    serial: each piece is segmented on its own); serial if the pool cannot start."""
    import os
    n = min(len(pieces), os.cpu_count() or 1)
    if n > 1:
        try:
            from concurrent.futures import ProcessPoolExecutor
            order = sorted(range(len(pieces)), key=lambda i: -len(pieces[i].faces))   # largest first
            with ProcessPoolExecutor(n) as ex:
                res = dict(zip(order, ex.map(_segment_piece, [(pieces[i], tol) for i in order])))
            return [res[i] for i in range(len(pieces))]
        except Exception:
            pass
    return [_segment_piece((c, tol)) for c in pieces]


def fit_rest(mesh, segs, tol, rmax, log=print):
    """The surface far from every strut, split into connected pieces; each large piece is segmented into
    primitives (stl2cad's Segmenter). Returns tubes (coaxial outer/inner cylinders between two planes), spheres
    (knobs), extra struts (cylinders of strut size) and the joint fillet radius (tori minor radii)."""
    t0 = time.time()
    far = np.flatnonzero(_seg_surface_distance(mesh.triangles_center, segs) > max(2 * tol, 0.4))
    sub = mesh.submesh([far], append=True)
    comps = sorted(sub.split(only_watertight=False), key=lambda c: -c.area)
    tubes, spheres, extra, fillet_r = [], [], [], []
    rmin_ = 0.5 * min(float(s["r"]) for s in segs)
    # large pieces (tubes, knobs and the chords converging into them) are all segmented; small pieces are joint
    # blends: only the first FILLET_PIECES, for the fillet radius (tori). The pieces are segmented in parallel.
    todo, n_small = [], 0
    for c in comps:
        if c.area < 20 * rmax * rmax:
            break
        if c.area < 200 * rmax * rmax:
            n_small += 1
            if n_small > FILLET_PIECES:
                continue
        todo.append(c)
    for c, prims in zip(todo, _segment_pieces(todo, 0.3 * tol)):
        area = {id(p): float(c.area_faces[p.faces].sum()) for p in prims}
        cyl = [p for p in prims if p.kind == "cylinder"]
        big = sorted([p for p in cyl if p.p["r"] > 3 * rmax and area[id(p)] > 30 * rmax * rmax],
                     key=lambda p: -area[id(p)])
        # tube: the two largest coaxial cylinders, ends at the extreme planes normal to their axis
        if len(big) >= 2:
            a = big[0].p["a"] / np.linalg.norm(big[0].p["a"])
            co = [p for p in big[1:] if abs(abs(p.p["a"] @ a) - 1) < 1e-3
                  and np.linalg.norm(np.cross(p.p["c"] - big[0].p["c"], a)) < 0.5 * tol + 0.05]
            if co:
                ro, ri = sorted([big[0].p["r"], co[0].p["r"]], reverse=True)
                h = (c.vertices - big[0].p["c"]) @ a
                pl = [p for p in prims if p.kind == "plane" and abs(abs(p.p["n"] @ a) - 1) < 1e-2]
                hs = sorted(((p.p["d"] * np.sign(p.p["n"] @ a)) - big[0].p["c"] @ a) for p in pl)
                h0, h1 = (hs[0], hs[-1]) if len(hs) >= 2 else (h.min(), h.max())
                tubes.append(dict(c=big[0].p["c"] + 0.5 * (h0 + h1) * a, a=a, ro=float(ro), ri=float(ri),
                                  h=float(h1 - h0)))
        n_tubes = len(tubes)
        found_sphere = False
        for p in prims:
            if p.kind == "sphere" and 2 * rmax < p.p["r"] < 8 * rmax and \
                    area[id(p)] > 4 * math.pi * p.p["r"] ** 2 * 0.08:
                spheres.append(dict(c=p.p["c"], r=float(p.p["r"]), area=area[id(p)]))
                found_sphere = True
        for p in prims:                           # joint blends: tori in pieces that hold no tube or knob
            if p.kind == "torus" and 0 < p.p["r"] < 4 * rmax and area[id(p)] > 2.0 and len(tubes) == n_tubes \
                    and not found_sphere:
                fillet_r.append((float(p.p["r"]), area[id(p)]))
            # a missed strut: radius in the struts' range (a small cylinder is an edge blend), long enough
            if p.kind == "cylinder" and rmin_ <= p.p["r"] <= 1.5 * rmax and area[id(p)] > 15 * p.p["r"] ** 2:
                a = p.p["a"] / np.linalg.norm(p.p["a"])
                X = c.triangles_center[p.faces]
                t = (X - p.p["c"]) @ a
                if t.max() - t.min() > MIN_LEN * p.p["r"]:
                    p0, p1 = p.p["c"] + t.min() * a, p.p["c"] + t.max() * a
                    # not a strut found already (its axis inside an existing strut)
                    mid = 0.5 * (p0 + p1)
                    dup = any(np.linalg.norm(np.cross(mid - q["p0"], (q["p1"] - q["p0"]) /
                                                      np.linalg.norm(q["p1"] - q["p0"]))) < q["r"] and
                              -q["r"] < (mid - q["p0"]) @ (q["p1"] - q["p0"]) / np.linalg.norm(q["p1"] - q["p0"])
                              < np.linalg.norm(q["p1"] - q["p0"]) + q["r"] for q in segs + extra)
                    if not dup:
                        extra.append(dict(p0=p0, p1=p1, r=float(p.p["r"]), r0=float(p.p["r"]), r1=float(p.p["r"]),
                                          n=len(p.faces)))
    rf = None
    if fillet_r:
        v = np.array(fillet_r)
        o = np.argsort(v[:, 0])
        cw = np.cumsum(v[o, 1])
        rf = float(v[o, 0][np.searchsorted(cw, 0.5 * cw[-1])])      # area-weighted median
    # one knob may give several sphere fits: the one over the largest area is kept
    sp = []
    for s in sorted(spheres, key=lambda s: -s["area"]):
        if all(np.linalg.norm(s["c"] - q["c"]) > q["r"] for q in sp):
            sp.append(s)
    # one tube found in two pieces: coaxial, same radii -> one tube over both extents
    tb = []
    for t in tubes:
        for q in tb:
            if abs(abs(t["a"] @ q["a"]) - 1) < 1e-3 and np.linalg.norm(np.cross(t["c"] - q["c"], q["a"])) < 0.1 \
                    and abs(t["ro"] - q["ro"]) < 0.05 * q["ro"] and abs(t["ri"] - q["ri"]) < 0.05 * q["ri"] \
                    and abs((t["c"] - q["c"]) @ q["a"]) < 0.5 * (t["h"] + q["h"]):
                if t["h"] > q["h"]:
                    q.update(c=t["c"], h=t["h"], a=t["a"])
                break
        else:
            tb.append(dict(t))
    tubes = tb
    log(f"  struts: rest of the part: {len(tubes)} tube(s), {len(sp)} sphere(s), {len(extra)} more strut(s), "
        f"joint fillet radius {rf if rf is None else round(rf, 3)} ({time.time() - t0:.1f}s)")
    return tubes, sp, extra, rf


def _sphere_fit(X):
    A = np.c_[2 * X, np.ones(len(X))]
    b = (X ** 2).sum(1)
    sol = np.linalg.lstsq(A, b, rcond=None)[0]
    c = sol[:3]
    return c, float(np.sqrt(max(sol[3] + c @ c, 0.0)))


def knobs(mesh, segs, free, spheres, rmax, log=print, seed=0):
    """Spheres (knobs) where free strut ends gather: RANSAC sphere on the surface far from the struts near the
    ends; radius between 2 and 8 x rmax. Spheres found by fit_rest are kept; a new one must not repeat them."""
    rng = np.random.default_rng(seed)
    V = mesh.vertices
    far = _seg_surface_distance(V, segs) > 0.3
    Vf, Nf = V[far], mesh.vertex_normals[far]
    from scipy.spatial import cKDTree
    T = cKDTree(Vf)
    out = list(spheres)
    for k, e in free:
        p = segs[k]["p1" if e else "p0"]
        if any(np.linalg.norm(p - q["c"]) < q["r"] + 2 * rmax for q in out):
            continue
        idx = T.query_ball_point(p, 10 * rmax)
        if len(idx) < 30:
            continue
        X, NX = Vf[idx], Nf[idx]
        best = None
        for _ in range(300):
            c, r = _sphere_fit(X[rng.choice(len(X), 4, replace=False)])
            if not 2 * rmax < r < 4.5 * rmax:
                continue
            w = X - c
            dist = np.linalg.norm(w, axis=1)
            inl = (np.abs(dist - r) < 0.03 * r) & ((w * NX).sum(1) > 0.9 * dist)   # outward normals: a knob
            if best is None or inl.sum() > best[0]:
                best = (inl.sum(), inl)
        if best is None or best[0] < 30:
            continue
        c, r = _sphere_fit(X[best[1]])
        # enough of the sphere present (an eighth of its area by vertex count around it)
        if best[0] < 0.1 * len(X) or np.linalg.norm(p - c) > r + 3 * rmax or \
                any(np.linalg.norm(c - q["c"]) < q["r"] for q in out):
            continue
        out.append(dict(c=c, r=r, area=0.0))
        log(f"  struts: knob r {r:.3f} at ({c[0]:.2f}, {c[1]:.2f}, {c[2]:.2f}) from {best[0]} vertices")
    return out


def attach(segs, free, tubes, spheres, nodes, log=print):
    """Free strut ends into a knob (to its centre) or through a tube wall into its bore (along the strut; the
    bore cut trims the end). Returns the ends still free."""
    left = []
    for k, e in free:
        s = segs[k]
        p = s["p1" if e else "p0"]
        q = s["p0" if e else "p1"]
        u = (p - q) / np.linalg.norm(p - q)
        done = False
        for sp in spheres:
            if np.linalg.norm(p - sp["c"]) < sp["r"] + 3 * s["r"]:
                s["p1" if e else "p0"] = np.asarray(sp["c"], float).copy()
                done = True
                break
        if not done:
            for t in tubes:
                a = t["a"]
                # step along the strut until the axis distance reaches the wall's middle radius
                best, in_wall = None, None
                e1 = np.cross(u, [1.0, 0, 0] if abs(u[0]) < 0.9 else [0, 1.0, 0])
                e1 /= np.linalg.norm(e1)
                e2 = np.cross(u, e1)
                th = np.linspace(0, 2 * np.pi, 16, endpoint=False)
                for step in np.linspace(-2 * s["r"], 6 * t["ro"], 400):
                    x = p + step * u
                    w = x - t["c"]
                    h = w @ a
                    rad = np.linalg.norm(w - h * a)
                    # through the wall to inside the bore by the strut's radius: the bore cut trims the end
                    # flush (an end in the wall shows part of its end disc outside the tube)
                    if abs(h) < 0.5 * t["h"] and rad < t["ri"] - s["r"]:
                        best = x
                        break
                    if in_wall is None:          # else: the first end whose whole disc lies within the tube's
                        # outer cylinder (wall or bore; the bore cut trims what lies in the bore)
                        D = x + s["r"] * (np.outer(np.cos(th), e1) + np.outer(np.sin(th), e2)) - t["c"]
                        hd = D @ a
                        rd = np.linalg.norm(D - np.outer(hd, a), axis=1)
                        if (np.abs(hd) < 0.5 * t["h"]).all() and (rd < t["ro"]).all():
                            in_wall = x
                if best is None:
                    best = in_wall
                if best is not None and np.linalg.norm(best - p) < 4 * t["ro"]:
                    s["p1" if e else "p0"] = best
                    done = True
                    break
        if not done:
            left.append((k, e))
    log(f"  struts: {len(free) - len(left)} free end(s) attached to knobs or tubes, {len(left)} left")
    return left


def _parts(res):
    """CadQuery solids per stage: struts (cylinders, cones), joint spheres, knobs, tubes (outer - inner)."""
    import cadquery as cq
    V = lambda p: cq.Vector(*map(float, p))
    struts_ = []
    for s in res["segs"]:
        d = np.asarray(s["p1"], float) - np.asarray(s["p0"], float)
        L = float(np.linalg.norm(d))
        if L < 1e-6:
            continue
        if s.get("taper") and abs(s["r0"] - s["r1"]) > 1e-9:
            struts_.append(cq.Solid.makeCone(float(s["r0"]), float(s["r1"]), L, V(s["p0"]), V(d / L)))
        else:
            struts_.append(cq.Solid.makeCylinder(float(s["r"]), L, V(s["p0"]), V(d / L)))
    sph = lambda sp: cq.Solid.makeSphere(float(sp["r"]), V(sp["c"]), angleDegrees1=-90, angleDegrees2=90)
    tubes = []
    for t in res["tubes"]:
        base = V(np.asarray(t["c"]) - 0.5 * t["h"] * np.asarray(t["a"]))
        outer = cq.Solid.makeCylinder(float(t["ro"]), float(t["h"]), base, V(t["a"]))
        inner = cq.Solid.makeCylinder(float(t["ri"]), float(t["h"]), base, V(t["a"]))
        tubes.append(outer.cut(inner))
    return dict(struts=struts_, joints=[sph(x) for x in res.get("joints", [])],
                knobs=[sph(x) for x in res["spheres"]], tubes=tubes)


def _fuse(base, tools, fuzzy=1e-5):
    """One fuzzy Boolean fuse of the tools onto the base (base None: the first tool)."""
    import cadquery as cq
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCP.TopTools import TopTools_ListOfShape
    tools = list(tools)
    if base is None:
        base, tools = tools[0], tools[1:]
    if not tools:
        return base
    A, B = TopTools_ListOfShape(), TopTools_ListOfShape()
    A.Append(base.wrapped)
    for p in tools:
        B.Append(p.wrapped)
    op = BRepAlgoAPI_Fuse()
    op.SetArguments(A)
    op.SetTools(B)
    op.SetFuzzyValue(fuzzy)
    op.SetRunParallel(True)
    op.Build()
    return cq.Shape.cast(op.Shape())


def _bore(t):
    import cadquery as cq
    V = lambda p: cq.Vector(*map(float, p))
    base = V(np.asarray(t["c"]) - (0.5 * t["h"] + 1.0) * np.asarray(t["a"]))
    return cq.Solid.makeCylinder(float(t["ri"]), float(t["h"]) + 2.0, base, V(t["a"]))


def build(res, log=print, fuzzy=1e-5, bores=True):
    """Fused CadQuery solid of struts (cylinders, cones), joint spheres, knobs (spheres) and tubes; no
    fillets; the tube bores cut (struts reaching into a bore are cut back by it)."""
    P = _parts(res)
    shape = _fuse(None, P["struts"] + P["joints"] + P["knobs"] + P["tubes"], fuzzy)
    try:
        cleaned = shape.clean()             # merging faces can fail or leave an invalid solid: then not used
        if cleaned.isValid():
            shape = cleaned
    except Exception:
        pass
    nb = len(res["tubes"]) if bores is True else (0 if bores is False else int(bores))
    for t in res["tubes"][:nb]:
        shape = shape.cut(_bore(t))
    return shape


def concave_edges(shape):
    """Edges where the material meets at a re-entrant angle (strut to strut, strut to knob or tube). With n_i
    the outward normal of face i and t the edge's tangent as oriented in face 1, face 1 lies on the side
    d = n_1 x t; the edge is re-entrant when face 1 rises above face 2 (d . n_2 > 0)."""
    from OCP.TopExp import TopExp, TopExp_Explorer
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_REVERSED
    import cadquery as cq
    mp = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape.wrapped, TopAbs_EDGE, TopAbs_FACE, mp)
    out = []
    for k in range(1, mp.Extent() + 1):
        e0 = mp.FindKey(k)
        faces = []
        for f in mp.FindFromIndex(k):
            if not any(f.IsSame(g) for g in faces):
                faces.append(f)
        if len(faces) != 2:
            continue
        e = cq.Edge(e0)
        if e.Length() < 1e-6:
            continue
        p = e.positionAt(0.5)
        f1, f2 = faces
        # the edge as it occurs in face 1 (its orientation there)
        ex = TopExp_Explorer(f1, TopAbs_EDGE)
        sign = 1.0
        while ex.More():
            if ex.Current().IsSame(e0):
                sign = -1.0 if ex.Current().Orientation() == TopAbs_REVERSED else 1.0
                break
            ex.Next()
        t = e.tangentAt(0.5) * sign
        n1, n2 = cq.Face(f1).normalAt(p), cq.Face(f2).normalAt(p)
        d = n1.cross(t)
        if d.dot(n2) > 1e-6:
            out.append(e)
    return out


def fillet_joints(shape, r, log=print):
    """Fillet the concave edges with radius r; if OCC fails on all at once, edge groups by position (one joint
    at a time) with r, then 0.7 r, then 0.5 r; groups that fail at every radius stay sharp."""
    edges = concave_edges(shape)
    if not edges:
        return shape, 0, 0
    try:
        f = shape.fillet(r, edges)
        if f.isValid() and f.Solids():
            log(f"  struts: {len(edges)} joint edges filleted with r {r:g}")
            return f, len(edges), 0
    except Exception:
        pass
    # group edges by joint: edges whose midpoints lie within 4 r of each other
    mids = np.array([e.positionAt(0.5).toTuple() for e in edges])
    from scipy.cluster.hierarchy import fcluster, linkage
    labels = fcluster(linkage(mids, "single"), 4 * r, "distance") if len(mids) > 1 else np.ones(1, int)
    done = failed = 0
    for g in np.unique(labels):
        grp = [edges[i] for i in np.flatnonzero(labels == g)]
        ok = False
        for rr in (r, 0.7 * r, 0.5 * r):
            # the edges of this joint in the current shape (the earlier fillets changed the topology)
            cur = concave_edges(shape)
            gm = np.array([e.positionAt(0.5).toTuple() for e in grp])
            sel = [e for e in cur if np.min(np.linalg.norm(gm - np.array(e.positionAt(0.5).toTuple()), axis=1))
                   < 0.5 * r]
            if not sel:
                break
            try:
                f = shape.fillet(rr, sel)
                if f.isValid() and f.Solids() and len(f.Solids()) == len(shape.Solids()):
                    shape = f
                    done += len(sel)
                    ok = True
                    break
            except Exception:
                pass
        if not ok:
            failed += len(grp)
    log(f"  struts: joint fillets: {done} edges filleted, {failed} left sharp")
    return shape, done, failed


def measure_radii(mesh, segs, stations=15, dirs=12):
    """Radius of each strut measured on the mesh: rays from its axis outward (dirs directions at each of the
    stations over the middle 60 % of its length, where no joint blend is); the median per station, then a
    linear fit from start to end (r0, r1) and its middle value (r). Stations widened by a crossing strut or a
    joint are dropped: the fit is repeated without stations more than 4 % of the radius above it."""
    rm = _intersector(mesh)
    O, D, owner = [], [], []
    for k, s in enumerate(segs):
        a, b = np.asarray(s["p0"], float), np.asarray(s["p1"], float)
        v = b - a
        L = np.linalg.norm(v)
        v = v / L
        e1 = np.cross(v, [1.0, 0, 0] if abs(v[0]) < 0.9 else [0, 1.0, 0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(v, e1)
        for j, f in enumerate(np.linspace(0.2, 0.8, stations)):
            c = a + f * L * v
            for th in np.linspace(0, 2 * np.pi, dirs, endpoint=False):
                O.append(c)
                D.append(np.cos(th) * e1 + np.sin(th) * e2)
                owner.append((k, f))
    O, D = np.array(O), np.array(D)
    dist = np.full(len(O), np.inf)
    for ch in np.array_split(np.arange(len(O)), max(1, len(O) // 2000)):
        loc, ri, _ = rm.intersects_location(O[ch], D[ch], multiple_hits=False)
        dist[ch[ri]] = np.linalg.norm(loc - O[ch][ri], axis=1)
    owner = np.array(owner)
    for k, s in enumerate(segs):
        m_ = owner[:, 0] == k
        fs, rs = [], []
        for f in np.unique(owner[m_, 1]):
            dd = dist[m_ & (owner[:, 1] == f)]
            dd = dd[np.isfinite(dd) & (dd < 3 * s["r"])]
            if len(dd) >= dirs // 2:
                fs.append(f)
                rs.append(float(np.median(dd)))
        fs, rs = np.array(fs), np.array(rs)
        keep = rs <= 1.04 * np.median(rs) if len(rs) else rs > 0
        for _ in range(5):
            if keep.sum() < 3:
                break
            A = np.c_[np.ones(keep.sum()), fs[keep]]
            c0, c1 = np.linalg.lstsq(A, rs[keep], rcond=None)[0]
            new = rs <= c0 + c1 * fs + 0.04 * (c0 + 0.5 * c1)
            if (new == keep).all():
                break
            keep = new
        if keep.sum() >= 3:
            s["r0"], s["r1"], s["r"] = float(c0), float(c0 + c1), float(c0 + 0.5 * c1)


def _stage_volumes(res):
    """Expected part volumes after each SolidWorks stage (CadQuery): struts, joint spheres, knobs, each tube,
    each bore. Each stage is fused onto the previous stage's solid (one Boolean per stage)."""
    P = _parts(res)
    shape = _fuse(None, P["struts"])
    out = [("struts", float(shape.Volume()))]
    if P["joints"]:
        shape = _fuse(shape, P["joints"])
        out.append(("joints", float(shape.Volume())))
    shape = _fuse(shape, P["knobs"])
    out.append(("knobs", float(shape.Volume())))
    for k, tb in enumerate(P["tubes"]):
        shape = _fuse(shape, [tb])
        out.append((f"tube {k + 1}", float(shape.Volume())))
    for k, t in enumerate(res["tubes"]):
        shape = shape.cut(_bore(t))
        out.append((f"bore {k + 1}", float(shape.Volume())))
    out.append(("bores", out[-1][1]))
    return out


def solidworks_code(res, name):
    """SolidWorks script (pywin32) for the strut model without joint fillets: per cylindrical strut a 3D sketch
    line and a circular-profile sweep; per tapered strut planes normal to its line at the ends, a circle on
    each, a loft along the line; per knob a half circle revolved about its diameter; per tube two circles
    extruded; the tube bores cut. Calls checked with sw_probe_struts.py and sw_probe_helix.py (SolidWorks 2025);
    every stage's volume is checked against stl2cad's CadQuery model."""
    import json
    stages = _stage_volumes(res)
    cyl = [dict(p0=[float(v) for v in s["p0"]], p1=[float(v) for v in s["p1"]], r=float(s["r"]))
           for s in res["segs"] if not (s.get("taper") and s["r0"] != s["r1"])]
    cones = [dict(p0=[float(v) for v in s["p0"]], p1=[float(v) for v in s["p1"]], r0=float(s["r0"]),
                  r1=float(s["r1"])) for s in res["segs"] if s.get("taper") and s["r0"] != s["r1"]]
    knobs = [dict(c=[float(v) for v in k["c"]], r=float(k["r"])) for k in res["spheres"]]
    joints = [dict(c=[float(v) for v in k["c"]], r=float(k["r"])) for k in res.get("joints", [])]
    tubes = []
    for t in res["tubes"]:
        a = np.asarray(t["a"], float)
        ax = int(np.argmax(np.abs(a)))
        if np.linalg.norm(a - np.sign(a[ax]) * np.eye(3)[ax]) > 1e-3:
            raise ValueError("struts: a tube axis is not along a coordinate axis")
        tubes.append(dict(c=[float(v) for v in t["c"]], axis=ax, ro=float(t["ro"]), ri=float(t["ri"]),
                          h=float(t["h"])))
    data = dict(struts=cyl, cones=cones, joints=joints, knobs=knobs, tubes=tubes, stages=stages,
                classes=[float(c) for c in res["classes"]])
    head = (f'"""{name}: truss of {len(cyl) + len(cones)} round struts ({len(cones)} tapered; radii '
            f'{min(s["r"] for s in res["segs"]):g}..{max(s["r"] for s in res["segs"]):g}), {len(joints)} joint sphere(s), {len(knobs)} knob(s), '
            f'{len(tubes)} tube(s).\n'
            "SolidWorks script (pywin32), generated by stl2cad (struts). Joint fillets are not made (open item).\n"
            "Stages: struts (3D sketch line + circular-profile sweep; tapered: loft between circles on planes\n"
            "normal to the line), joint spheres at strut ends (revolve), knobs (revolve), tubes (extrude), bores (cut); each stage's volume is checked\n"
            'against the CadQuery model."""\n')
    return head + "DATA_JSON = r'''" + json.dumps(data) + "'''\n" + _SW_STRUTS


_SW_STRUTS = r'''import json
import math
import sys
import time

import pythoncom
import win32com.client
from win32com.client import VARIANT

NOTHING = VARIANT(pythoncom.VT_DISPATCH, None)
MM = 0.001
DATA = json.loads(DATA_JSON)
FM_SWEEP = 17                    # swFeatureNameID_e.swFmSweep (SolidWorks 2025)
BASE = {0: "Right Plane", 1: "Top Plane", 2: "Front Plane"}     # the plane normal to x, y, z
REL_TOL = 1e-3                   # stage volumes against the CadQuery model


def call0(obj, name):
    r = getattr(obj, name)
    if isinstance(r, win32com.client.CDispatch) or not callable(r):
        return r
    return r()


def to_sk(T, p):
    """Model point (mm) -> sketch coordinates (m), by the sketch's ModelToSketchTransform."""
    x = [c * MM for c in p]
    return tuple((sum(x[i] * T[3 * i + j] for i in range(3)) + T[9 + j]) * T[12] for j in range(3))


class SW:
    def __init__(self):
        self.app = win32com.client.Dispatch("SldWorks.Application")
        self.app.Visible = True
        try:
            self.app.SetUserPreferenceToggle(10, False)
        except Exception:
            pass
        self.model = call0(self.app, "NewPart")
        if self.model is None:
            self.model = self.app.NewDocument(self.app.GetUserPreferenceStringValue(8), 0, 0, 0)
        self.fm = self.model.FeatureManager
        self.quiet(True)

    def quiet(self, on):
        for obj, attr in ((lambda: self.app, "CommandInProgress"), (lambda: self.fm, "EnableFeatureTree"),
                          (lambda: call0(self.model, "ActiveView"), "EnableGraphicsUpdate")):
            try:
                setattr(obj(), attr, on if attr == "CommandInProgress" else not on)
            except Exception:
                pass

    def last(self):
        return self.model.FeatureByPositionReverse(0)

    def volume(self):
        self.model.ForceRebuild3(False)
        mp = call0(self.model.Extension, "CreateMassProperty")
        return 0.0 if mp is None else mp.Volume / MM ** 3

    def delete(self, *feats):
        for f in feats:
            if f is not None:
                self.model.ClearSelection2(True)
                f.Select2(False, 0)
                self.model.Extension.DeleteSelection2(1)

    # ---- sketches ----
    def line3d(self, p0, p1):
        """A 3D sketch with one line (mm); the first sketch after NewPart may not open and the first line may
        not be made (sw_probe_struts): both are retried."""
        sm = self.model.SketchManager
        self.model.ClearSelection2(True)
        for _ in range(3):                      # an odd number of toggles leaves the sketch open
            sm.Insert3DSketch(True)
            if call0(sm, "ActiveSketch") is not None or call0(self.model, "GetActiveSketch2") is not None:
                break
        sm.AddToDB = True
        ln = None
        for _ in range(4):
            ln = sm.CreateLine(*(v * MM for v in p0), *(v * MM for v in p1))
            if ln is not None:
                break
        sm.AddToDB = False
        sm.Insert3DSketch(True)
        if ln is None:
            raise RuntimeError(f"no 3D sketch line {p0} -> {p1}")
        return self.last(), ln

    def open_sketch(self, plane):
        sm = self.model.SketchManager
        for _ in range(4):
            self.model.ClearSelection2(True)
            plane.Select2(False, 0)
            sm.InsertSketch(True)
            time.sleep(0.1)
            sk = call0(sm, "ActiveSketch") or call0(self.model, "GetActiveSketch2")
            if sk is not None:
                return sm, sk, list(sk.ModelToSketchTransform.ArrayData)
        raise RuntimeError(f"no sketch opened on {plane.Name}")

    def plane_at(self, axis, point, name):
        """The base plane normal to the axis (0 x, 1 y, 2 z), or a plane offset from it through the point (mm);
        the side is checked through a sketch's transform."""
        base = self.model.FeatureByName(BASE[axis])
        d = point[axis]
        if abs(d) < 1e-9:
            return base
        for flip in (0, 256):
            self.model.ClearSelection2(True)
            base.Select2(False, 0)
            p = self.fm.InsertRefPlane(8 | flip, abs(d) * MM, 0, 0, 0, 0)
            if p is None:
                continue
            f = self.last()
            sm, sk, T = self.open_sketch(f)
            sm.InsertSketch(True)
            s_ = self.last()
            if call0(s_, "GetTypeName2") == "ProfileFeature":
                self.delete(s_)
            if abs(to_sk(T, point)[2]) < 1e-7:
                f.Name = name
                return f
            self.delete(f)
        raise RuntimeError(f"could not make {name}")

    # ---- features ----
    def overshoot(self, f, s):
        """How far (mm) the feature's faces run along the strut's own cylinder past its ends, beyond the
        joint allowance of two radii (tessellated faces: points at distance r from the axis line); SolidWorks
        sometimes makes a circular-profile sweep longer than its path (sw_probe_strutlines: strut 20 swept
        11 mm past its start)."""
        p0 = [float(v) for v in s["p0"]]
        d = [float(b) - a for a, b in zip(p0, s["p1"])]
        L = math.sqrt(sum(x * x for x in d))
        u = [x / L for x in d]
        r = s["r"]
        worst = 0.0
        for fc in call0(f, "GetFaces") or ():
            tri = fc.GetTessTriangles(True)
            if not tri:
                continue
            for i in range(0, len(tri) - 2, 3):
                w = [tri[i + j] / MM - p0[j] for j in range(3)]
                t = sum(a * b for a, b in zip(w, u))
                rad = math.sqrt(max(sum(a * a for a in w) - t * t, 0.0))
                if abs(rad - r) < 0.02 * r:
                    worst = max(worst, -t - 2 * r, t - L - 2 * r)
        return worst

    def strut(self, s, k):
        """Circular-profile sweep along a 3D sketch line; a sweep that runs past its line is deleted and the
        strut is made as a loft between equal circles on planes normal to the line at its ends."""
        sk, ln = self.line3d(s["p0"], s["p1"])
        sk.Name = f"Strut {k} line"
        self.model.ClearSelection2(True)
        sk.Select2(True, 4)
        d = self.fm.CreateDefinition(FM_SWEEP)
        d.CircularProfile = True
        d.CircularProfileDiameter = 2 * s["r"] * MM
        try:
            d.Merge = True
        except Exception:
            pass
        f = self.fm.CreateFeature(d)
        if f is not None:
            try:
                over = self.overshoot(f, s)
            except Exception:
                over = 0.0
            if over <= 0.01:
                f.Name = f"Strut {k}"
                return
            print(f"  strut {k}: the sweep runs {over:.2f} mm past its line; made as a loft")
        self.delete(f, sk)
        if not self.cone(dict(p0=s["p0"], p1=s["p1"], r0=s["r"], r1=s["r"]), k):
            raise RuntimeError(f"strut {k}: no sweep or loft")

    def cone(self, c, k):
        """Tapered strut: planes normal to its line at the ends (perpendicular + coincident), a circle on each,
        a loft from one to the other along the line (then without the line); if neither is made, everything is
        made again from the other end."""
        for cc in (c, dict(p0=c["p1"], p1=c["p0"], r0=c["r1"], r1=c["r0"])):
            made = self.cone_once(cc, k)
            if made:
                return True
        raise RuntimeError(f"strut {k}: no loft")

    def cone_once(self, c, k):
        sk, ln = self.line3d(c["p0"], c["p1"])
        sk.Name = f"Strut {k} line"
        secs, planes = [], []
        for end, pt, ctr, r in ((0, call0(ln, "GetStartPoint2"), c["p0"], c["r0"]),
                                (1, call0(ln, "GetEndPoint2"), c["p1"], c["r1"])):
            self.model.ClearSelection2(True)
            sd0 = call0(self.model.SelectionManager, "CreateSelectData")
            sd0.Mark = 0
            sd1 = call0(self.model.SelectionManager, "CreateSelectData")
            sd1.Mark = 1
            ln.Select4(False, sd0)
            pt.Select4(True, sd1)
            if self.fm.InsertRefPlane(2, 0, 4, 0, 0, 0) is None:
                raise RuntimeError(f"strut {k}: no plane at end {end}")
            pl = self.last()
            pl.Name = f"Strut {k} end {end}"
            planes.append(pl)
            sm, s_, T = self.open_sketch(pl)
            u = to_sk(T, ctr)
            sm.AddToDB = True
            sm.CreateCircleByRadius(u[0], u[1], 0.0, r * MM)
            sm.AddToDB = False
            sm.InsertSketch(True)
            secs.append(self.last())
        v0 = self.volume()
        for label, centerline in (("along its line", True), ("between the end circles", False)):
            self.model.ClearSelection2(True)
            secs[0].Select2(False, 1)
            secs[1].Select2(True, 1)
            if centerline:
                sk.Select2(True, 4)
            f = self.fm.InsertProtrusionBlend2(False, True, False, 1.0, 0, 0, 1.0, 1.0, True, True, False, 0.0, 0.0,
                                               0, True, True, True, 0)
            if f is not None and self.volume() > v0:
                f.Name = f"Strut {k}"
                if not centerline:
                    print(f"  strut {k}: loft {label}")
                return True
            if f is None:
                print(f"  strut {k}: no loft {label}")
            self.delete(f)
        self.delete(*secs, *planes, sk)
        return False

    def knob(self, kb, k, label="Knob"):
        """Sphere: on the plane parallel to the Front Plane through its centre, a half circle closed by its
        diameter and a centerline on the diameter, revolved 360 degrees."""
        c, r = kb["c"], kb["r"]
        pl = self.plane_at(2, c, f"{label} {k} plane")
        sm, sk, T = self.open_sketch(pl)
        u = to_sk(T, c)
        sm.AddToDB = True
        sm.CreateArc(u[0], u[1], 0.0, u[0], u[1] - r * MM, 0.0, u[0], u[1] + r * MM, 0.0, 1)
        sm.CreateLine(u[0], u[1] - r * MM, 0.0, u[0], u[1] + r * MM, 0.0)
        cl = sm.CreateCenterLine(u[0], u[1] - r * MM, 0.0, u[0], u[1] + r * MM, 0.0)
        sm.AddToDB = False
        sm.InsertSketch(True)
        prof = self.last()
        prof.Name = f"{label} {k} sketch"
        v0 = self.volume()
        for mark in (None, 4, 16):
            self.model.ClearSelection2(True)
            prof.Select2(False, 0)
            if mark is not None and cl is not None:
                sd = call0(self.model.SelectionManager, "CreateSelectData")
                sd.Mark = mark
                cl.Select4(True, sd)
            f = self.fm.FeatureRevolve2(True, True, False, False, False, False, 0, 0, 2 * math.pi, 0.0,
                                        False, False, 0.0, 0.0, 0, 0.0, 0.0, True, True, True)
            if f is not None and self.volume() > v0:
                f.Name = f"{label} {k}"
                return
            self.delete(f)
        raise RuntimeError(f"{label.lower()} {k}: no revolve")

    def extrude_circles(self, t, radii, start, depth, cut, name, expected):
        """Circles (radii) about the tube axis on a plane normal to it at axial position start (mm from the
        tube centre), extruded (or cut) depth along the axis; the direction giving the expected volume is kept."""
        ax = t["axis"]
        p = list(t["c"])
        p[ax] += start
        pl = self.plane_at(ax, p, f"{name} plane")
        sm, sk, T = self.open_sketch(pl)
        u = to_sk(T, p)
        sm.AddToDB = True
        for r in radii:
            sm.CreateCircleByRadius(u[0], u[1], 0.0, r * MM)
        sm.AddToDB = False
        sm.InsertSketch(True)
        prof = self.last()
        prof.Name = f"{name} sketch"
        for rev in (False, True):
            self.model.ClearSelection2(True)
            prof.Select2(False, 0)
            if cut:
                f = self.fm.FeatureCut4(True, False, rev, 0, 0, depth * MM, 0.0, False, False, False, False,
                                        0.0, 0.0, False, False, False, False, False, True, True,
                                        False, False, False, 0, 0.0, False, False)
            else:
                f = self.fm.FeatureExtrusion3(True, False, rev, 0, 0, depth * MM, 0.0, False, False, False, False,
                                              0.0, 0.0, False, False, False, False, True, True, True, 0, 0.0, False)
            if f is None:
                continue
            v = self.volume()
            if abs(v - expected) <= REL_TOL * expected:
                f.Name = name
                return v
            print(f"  {name}: reverse={rev}: volume {v:.3f} mm^3, expected {expected:.3f}; other direction")
            self.delete(f)
        raise RuntimeError(f"{name}: no direction gives the expected volume")


def hide_construction(model):
    """Hide the reference planes, axes and sketches the script made (BlankRefGeom, BlankSketch), so that only
    the solid shows. They are hidden, not suppressed: suppressing a plane would suppress what is built on it."""
    n = {"planes": 0, "axes": 0, "sketches": 0}
    base = {"Front Plane", "Top Plane", "Right Plane"}
    f = call0(model, "FirstFeature")
    feats = []
    while f is not None and len(feats) < 100000:
        feats.append(f)
        f = call0(f, "GetNextFeature")
    for f in feats:
        kind = call0(f, "GetTypeName2")
        try:
            if kind in ("RefPlane", "RefAxis") and f.Name not in base:
                model.ClearSelection2(True)
                f.Select2(False, 0)
                model.BlankRefGeom()
                n["planes" if kind == "RefPlane" else "axes"] += 1
            elif kind in ("ProfileFeature", "3DProfileFeature"):
                model.ClearSelection2(True)
                f.Select2(False, 0)
                model.BlankSketch()
                n["sketches"] += 1
        except Exception:
            pass
        # sketches absorbed by a feature (a sweep's path, a loft's sections) are its sub-features
        sub = call0(f, "GetFirstSubFeature")
        guard = 0
        while sub is not None and guard < 1000:
            guard += 1
            try:
                if call0(sub, "GetTypeName2") in ("ProfileFeature", "3DProfileFeature"):
                    model.ClearSelection2(True)
                    sub.Select2(False, 0)
                    model.BlankSketch()
                    n["sketches"] += 1
            except Exception:
                pass
            sub = call0(sub, "GetNextSubFeature")
    model.ClearSelection2(True)
    print(f"hidden: {n['planes']} plane(s), {n['axes']} axis(es), {n['sketches']} sketch(es)")


def check(sw, stage, expected):
    v = sw.volume()
    ok = abs(v - expected) <= REL_TOL * expected
    print(f"{stage}: volume {v:.3f} mm^3 (expected {expected:.3f}, {100 * (v / expected - 1):+.3f} %)"
          + ("" if ok else "  <-- differs"))
    return v


def build():
    sw = SW()
    st = dict(DATA["stages"])
    n = len(DATA["struts"]) + len(DATA["cones"])
    print(f"connected to SolidWorks; {n} struts, {len(DATA['knobs'])} knobs, {len(DATA['tubes'])} tubes")
    failed = []

    def attempt(label, fn, *args):
        """A failed feature is reported and left out; the build goes on (the volume checks show the gap)."""
        try:
            fn(*args)
        except Exception as ex:
            failed.append(label)
            print(f"  {label}: FAILED ({type(ex).__name__}: {ex}); left out")
    k = 0
    for s in DATA["struts"]:
        k += 1
        attempt(f"strut {k}", sw.strut, s, k)
    for c in DATA["cones"]:
        k += 1
        attempt(f"strut {k} (tapered)", sw.cone, c, k)
    check(sw, "struts", st["struts"])
    for j, kb in enumerate(DATA.get("joints", []), 1):
        attempt(f"joint {j}", sw.knob, kb, j, "Joint")
    if DATA.get("joints"):
        check(sw, "joints", st["joints"])
    for j, kb in enumerate(DATA["knobs"], 1):
        attempt(f"knob {j}", sw.knob, kb, j)
    if DATA["knobs"]:
        check(sw, "knobs", st["knobs"])
    for j, t in enumerate(DATA["tubes"], 1):
        # with a feature left out the stage volumes no longer match: the tube's own volume is checked instead
        if failed:
            exp = sw.volume() + st[f"tube {j}"] - (st[f"tube {j - 1}"] if j > 1 else
                                                    (st["knobs"] if DATA["knobs"] else
                                                     st["joints"] if DATA.get("joints") else st["struts"]))
        else:
            exp = st[f"tube {j}"]
        attempt(f"tube {j}", sw.extrude_circles, t, (t["ro"], t["ri"]), -0.5 * t["h"], t["h"], False, f"Tube {j}", exp)
    prev = st[f"tube {len(DATA['tubes'])}"] if DATA["tubes"] else None
    for j, t in enumerate(DATA["tubes"], 1):
        if abs(st[f"bore {j}"] - prev) <= 1e-6 * prev:
            print(f"Bore {j}: no strut reaches into the bore; no cut")
            continue
        sw.extrude_circles(t, (t["ri"],), -0.5 * t["h"] - 1.0, t["h"] + 2.0, True, f"Bore {j}", st[f"bore {j}"])
        prev = st[f"bore {j}"]
    check(sw, "part", st["bores"])
    if failed:
        print(f"left out: {', '.join(failed)}")
    try:
        hide_construction(sw.model)
    except Exception as ex:
        print(f"construction geometry not hidden: {type(ex).__name__}: {ex}")
    sw.quiet(False)
    try:
        sw.model.GraphicsRedraw2()
        sw.app.SetUserPreferenceToggle(10, True)
    except Exception:
        pass
    return sw


if __name__ == "__main__":
    sw = build()
    if len(sys.argv) > 1:
        sw.model.SaveAs3(sys.argv[1], 0, 0)
'''


def cadquery_code(res, name):
    """CadQuery script of the strut model (no joint fillets): one fused solid, exported as STEP."""
    import json
    data = dict(segs=[dict(p0=[float(v) for v in s["p0"]], p1=[float(v) for v in s["p1"]], r=float(s["r"]),
                           r0=float(s.get("r0", s["r"])), r1=float(s.get("r1", s["r"])),
                           taper=bool(s.get("taper") and s["r0"] != s["r1"])) for s in res["segs"]],
                spheres=[dict(c=[float(v) for v in k["c"]], r=float(k["r"]))
                         for k in list(res.get("joints", [])) + list(res["spheres"])],
                tubes=[dict(c=[float(v) for v in t["c"]], a=[float(v) for v in t["a"]], ro=float(t["ro"]),
                            ri=float(t["ri"]), h=float(t["h"])) for t in res["tubes"]])
    return (f'"""{name}: truss of {len(data["segs"])} round struts, {len(data["spheres"])} knob(s), '
            f'{len(data["tubes"])} tube(s); radii {min(x["r"] for x in data["segs"]):g}..{max(x["r"] for x in data["segs"]):g}.\n'
            'Generated by stl2cad (struts). Joint fillets are not modelled (OCC cannot fillet these joints)."""\n'
            "import json\n\nimport cadquery as cq\nfrom OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse\n"
            "from OCP.TopTools import TopTools_ListOfShape\n\n"
            "DATA = json.loads(r'''" + json.dumps(data) + "''')\n" + r'''

def V(p):
    return cq.Vector(*p)


def build():
    parts = []
    for s in DATA["segs"]:
        d = V(s["p1"]) - V(s["p0"])
        L = d.Length
        if s["taper"]:
            parts.append(cq.Solid.makeCone(s["r0"], s["r1"], L, V(s["p0"]), d.normalized()))
        else:
            parts.append(cq.Solid.makeCylinder(s["r"], L, V(s["p0"]), d.normalized()))
    for k in DATA["spheres"]:
        parts.append(cq.Solid.makeSphere(k["r"], V(k["c"]), angleDegrees1=-90, angleDegrees2=90))
    for t in DATA["tubes"]:
        a = V(t["a"])
        base = V(t["c"]) - a * (0.5 * t["h"])
        parts.append(cq.Solid.makeCylinder(t["ro"], t["h"], base, a).cut(cq.Solid.makeCylinder(t["ri"], t["h"], base, a)))
    A, B = TopTools_ListOfShape(), TopTools_ListOfShape()
    A.Append(parts[0].wrapped)
    for p in parts[1:]:
        B.Append(p.wrapped)
    op = BRepAlgoAPI_Fuse()
    op.SetArguments(A)
    op.SetTools(B)
    op.SetFuzzyValue(1e-5)
    op.Build()
    shape = cq.Shape.cast(op.Shape())
    try:
        cleaned = shape.clean()             # merging faces can fail or leave an invalid solid: then not used
        if cleaned.isValid():
            shape = cleaned
    except Exception:
        pass
    for t in DATA["tubes"]:
        a = V(t["a"])
        shape = shape.cut(cq.Solid.makeCylinder(t["ri"], t["h"] + 2.0, V(t["c"]) - a * (0.5 * t["h"] + 1.0), a))
    return shape


result = build()

if __name__ == "__main__":
    import sys
    out = sys.argv[1] if len(sys.argv) > 1 else "''' + name + r'''.step"
    cq.exporters.export(result, out)
    print("valid", result.isValid(), "solids", len(result.Solids()), "volume", round(result.Volume(), 3))
''')


def _end_radius(s, e):
    return float(s["r1" if e else "r0"]) if s.get("taper") else float(s["r"])


def _seg_point(q, s):
    """Nearest point on a strut's axis segment, its parameter and distance."""
    a, b = np.asarray(s["p0"], float), np.asarray(s["p1"], float)
    ab = b - a
    t = float(np.clip((q - a) @ ab / (ab @ ab), 0.0, 1.0))
    x = a + t * ab
    return x, t, float(np.linalg.norm(q - x))


def _free_ends(segs, spheres, tubes, eps=1e-3):
    """Ends that meet nothing: no other end, knob, tube or strut side at the end point."""
    out = []
    for k, s in enumerate(segs):
        for e in (0, 1):
            p = np.asarray(s["p1" if e else "p0"], float)
            hit = any(np.linalg.norm(p - np.asarray(sp["c"], float)) < sp["r"] for sp in spheres)
            for t in tubes:
                w = p - np.asarray(t["c"], float)
                h = w @ np.asarray(t["a"], float)
                hit |= abs(h) < 0.5 * t["h"] and np.linalg.norm(w - h * np.asarray(t["a"], float)) < t["ro"]
            for j, o in enumerate(segs):
                if hit:
                    break
                if j != k and _seg_point(p, o)[2] < max(eps, 0.5 * float(o["r"])):
                    hit = True
            if not hit:
                out.append((k, e))
    return out


def extend_free_ends(mesh, segs, spheres, tubes, log=print):
    """A free end whose axis runs on inside the part into a knob or onto another strut's axis is extended
    there (to the knob's centre, or the nearest point of the other axis): the strut was cut short where its
    centre line was lost near the joint. The axis must stay inside the part by a quarter radius all the way."""
    import trimesh
    pq = trimesh.proximity.ProximityQuery(mesh)
    n = 0
    for k, e in _free_ends(segs, spheres, tubes):
        s = segs[k]
        p = np.asarray(s["p1" if e else "p0"], float)
        q = np.asarray(s["p0" if e else "p1"], float)
        L = float(np.linalg.norm(p - q))
        u = (p - q) / L
        r = _end_radius(s, e)
        steps = np.arange(0.25 * r, max(15 * r, 0.5 * L), 0.25 * r)
        sd = pq.signed_distance(p + np.outer(steps, u))
        target = None
        for lam, d in zip(steps, sd):
            if d < 0.25 * r:
                break
            x = p + lam * u
            for sp in spheres:
                if np.linalg.norm(x - np.asarray(sp["c"], float)) < sp["r"]:
                    target = np.asarray(sp["c"], float).copy()
                    break
            if target is not None:
                break
            for j, o in enumerate(segs):
                if j == k:
                    continue
                y, _, dist = _seg_point(x, o)
                if dist < 0.5 * float(o["r"]):
                    target = y
                    break
            if target is not None:
                break
        if target is None:
            continue
        # the whole way to the target inside the part
        L2 = float(np.linalg.norm(target - p))
        ts = np.linspace(0, L2, max(4, int(L2 / (0.25 * r))))
        if (pq.signed_distance(p + np.outer(ts, (target - p) / max(L2, 1e-12))) < 0.25 * r).any():
            continue
        s["p1" if e else "p0"] = target
        n += 1
    if n:
        log(f"  struts: {n} free end(s) extended along the strut to a knob or another strut")
    return segs


def joint_spheres(mesh, segs, spheres, tubes, log=print, eps=1e-3):
    """A sphere at every strut end that meets another strut: the part's joint ball, JOINT_SCALE times the
    largest sphere centred at the joint inside the part, and at least the end radius (the end discs of struts meeting at a joint
    otherwise show as blades beyond the other struts). Ends in knobs or tubes, free ends, and spheres that a strut through the point already
    contains or that leave the part are skipped."""
    import trimesh
    pq = trimesh.proximity.ProximityQuery(mesh)
    free = set(_free_ends(segs, spheres, tubes))
    pts = []
    for k, s in enumerate(segs):
        for e in (0, 1):
            if (k, e) in free:
                continue
            p = np.asarray(s["p1" if e else "p0"], float)
            if any(np.linalg.norm(p - np.asarray(sp["c"], float)) < sp["r"] for sp in spheres):
                continue
            inside_tube = False
            for t in tubes:
                w = p - np.asarray(t["c"], float)
                h = w @ np.asarray(t["a"], float)
                inside_tube |= abs(h) < 0.5 * t["h"] + _end_radius(s, e) and \
                    np.linalg.norm(w - h * np.asarray(t["a"], float)) < t["ro"] + _end_radius(s, e)
            if inside_tube:
                continue
            for c in pts:
                if np.linalg.norm(c["c"] - p) < eps:
                    c["r"] = max(c["r"], _end_radius(s, e))
                    break
            else:
                pts.append(dict(c=p.copy(), r=_end_radius(s, e)))
    # the part's joint ball: JOINT_SCALE times the largest sphere centred at the joint that stays inside the
    # part (the distance from the joint to the surface), rounded down to 0.01 mm, at least the end radius
    if pts:
        dist = pq.signed_distance(np.array([c["c"] for c in pts]))
        for c, d_ in zip(pts, dist):
            c["r"] = max(c["r"], math.floor(100 * JOINT_SCALE * float(d_)) / 100)
    out, n_in, n_out = [], 0, 0
    fib = np.arange(64) + 0.5
    phi, th = np.arccos(1 - 2 * fib / 64), np.pi * (1 + 5 ** 0.5) * fib
    unit = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    for c in pts:
        covered = False
        for o in segs:
            y, t, dist = _seg_point(c["c"], o)
            L = float(np.linalg.norm(np.asarray(o["p1"], float) - np.asarray(o["p0"], float)))
            ro = float(o["r"]) if not o.get("taper") else min(float(o["r0"]), float(o["r1"]))
            if dist < eps and ro >= c["r"] - 1e-9 and c["r"] / L < t < 1 - c["r"] / L:
                covered = True
                break
        if covered:
            n_in += 1
            continue
        if (pq.signed_distance(c["c"] + c["r"] * unit) < -0.1 * c["r"]).any():
            n_out += 1
            continue
        out.append(c)
    log(f"  struts: {len(out)} joint sphere(s) at strut ends ({n_in} inside a strut, {n_out} leaving the part)")
    return out


def prune_outside(mesh, segs, log=print):
    """Struts whose axis leaves the mesh (a strut grown past a knob, an end pushed too far) are cut where the
    axis is outside for more than two radii; inside pieces shorter than MIN_LEN radii are dropped."""
    import trimesh
    pq = trimesh.proximity.ProximityQuery(mesh)
    out_segs, n_cut = [], 0
    for s in segs:
        p0, p1 = np.asarray(s["p0"], float), np.asarray(s["p1"], float)
        L = float(np.linalg.norm(p1 - p0))
        v = (p1 - p0) / L
        ts = np.linspace(0, L, max(10, int(L / (0.25 * s["r"]))))
        inside = pq.signed_distance(p0 + np.outer(ts, v)) > -0.1 * s["r"]
        if inside.all():
            out_segs.append(s)
            continue
        # runs of inside samples; gaps (outside runs) shorter than 2 r are bridged
        runs, start = [], None
        for i, ok in enumerate(inside):
            if ok and start is None:
                start = i
            if not ok and start is not None:
                runs.append([start, i - 1])
                start = None
        if start is not None:
            runs.append([start, len(ts) - 1])
        merged = []
        for r_ in runs:
            if merged and ts[r_[0]] - ts[merged[-1][1]] < 2 * s["r"]:
                merged[-1][1] = r_[1]
            else:
                merged.append(r_)
        if len(merged) == 1 and merged[0][0] == 0 and merged[0][1] == len(ts) - 1:
            out_segs.append(s)
            continue
        n_cut += 1
        for a, b in merged:
            t0, t1 = ts[a], ts[b]
            if t1 - t0 < MIN_LEN * s["r"]:
                continue
            q = dict(s)
            q["p0"], q["p1"] = p0 + t0 * v, p0 + t1 * v
            if s.get("taper"):
                q["r0"] = s["r0"] + (s["r1"] - s["r0"]) * t0 / L
                q["r1"] = s["r0"] + (s["r1"] - s["r0"]) * t1 / L
            out_segs.append(q)
    if n_cut:
        log(f"  struts: {n_cut} strut(s) cut where their axis leaves the part")
    return out_segs


def trim_free_ends(mesh, segs, free, log=print):
    """A free strut end (no node, knob or tube) is pulled back along the strut until its end disc lies inside
    the mesh (every disc point within 0.05 r of the surface or inside): an end left in mid-air otherwise pokes
    out of the part."""
    import trimesh
    pq = trimesh.proximity.ProximityQuery(mesh)
    n = 0
    for k, e in free:
        s = segs[k]
        p = np.asarray(s["p1" if e else "p0"], float)
        q = np.asarray(s["p0" if e else "p1"], float)
        L = float(np.linalg.norm(p - q))
        u = (p - q) / L
        e1 = np.cross(u, [1.0, 0, 0] if abs(u[0]) < 0.9 else [0, 1.0, 0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(u, e1)
        th = np.linspace(0, 2 * np.pi, 16, endpoint=False)
        r = s["r1" if e else "r0"] if s.get("taper") else s["r"]
        best = p
        for back in np.linspace(0, 0.5 * L, 60):
            c = p - back * u
            D = c + r * (np.outer(np.cos(th), e1) + np.outer(np.sin(th), e2))
            if (pq.signed_distance(np.vstack([D, c[None]])) > -0.05 * r).all():
                best = c
                break
        else:
            continue
        if np.linalg.norm(best - p) > 1e-9:
            s["p1" if e else "p0"] = best
            n += 1
    if n:
        log(f"  struts: {n} free end(s) pulled back inside the part")
