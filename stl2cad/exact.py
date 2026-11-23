"""Exact reconstruction for CAD-like parts (prismatic, turned, filleted).

Pipeline: segment the mesh into analytic surfaces -> snap to design intent ->
(optionally) symmetrize -> build a B-rep by cell arrangement:
  1. every fitted surface becomes a tool face, extended a little past its region,
  2. a box around the part is split by all tool faces (OCC general fuse splitter),
  3. each cell is classified inside/outside the STL (majority vote, symmetric
     vote over the symmetry group when symmetry is imposed),
  4. inside cells are glued and coplanar/co-cylindrical faces are unified.
Because the true boundary lies on the union of the tool faces, every cell is
either fully inside or fully outside the part.
"""

import time

import numpy as np
from cadquery import Shape, Vector

from . import brep_helpers as BH
from .core import Result
from .primitives import perp_basis, unit
from .segment import (Segmenter, cones_to_tori, refit_odd_regions, symmetric_consensus, decompose_revolutions, spheres_to_cylinders, facets_to_cylinders, reassign_faces, simplify_degenerate, snap_directions, snap_grid, snap_tangency,
                      summarize, symmetrize, unify_parallel)
from .symmetry import group_matrices

TWO_PI = 2 * np.pi


def tool_spec(pr, margin):
    """Plain-data description of the tool surface of a primitive: its region plus a margin."""
    k, p, X = pr.kind, pr.p, getattr(pr, "extent_pts", pr.pts)
    margin = margin + getattr(pr, "extra_margin", 0.0)      # neighbours of removed blends reach the sharp edge
    if k == "plane":
        n = p["n"]
        u, v = perp_basis(n)
        c = X.mean(0)
        o = c - (c @ n - p["d"]) * n
        U, V = (X - o) @ u, (X - o) @ v
        return dict(kind="plane", n=n, d=float(p["d"]), o=o, u=u,
                    rect=(U.min() - margin, U.max() + margin, V.min() - margin, V.max() + margin))
    if k == "cylinder":
        h = (X - p["c"]) @ p["a"]
        h_reg = (pr.pts - p["c"]) @ p["a"]
        return dict(kind="cylinder", c=p["c"], a=p["a"], r=float(p["r"]), h=(h.min() - margin, h.max() + margin),
                    hole=bool(getattr(pr, "is_hole", False)), h_region=(h_reg.min(), h_reg.max()))
    if k == "cone":
        h = (X - p["v"]) @ p["a"]
        return dict(kind="cone", v=p["v"], a=p["a"], t=float(p["t"]), h=(max(h.min() - margin, 0.0), h.max() + margin))
    if k == "sphere":
        return dict(kind="sphere", c=p["c"], r=float(p["r"]))
    if k == "torus":
        if p["r"] >= p["R"]:
            return None
        return dict(kind="torus", c=p["c"], a=p["a"], R=float(p["R"]), r=float(p["r"]))
    if k == "extrusion":
        from scipy.interpolate import splev
        closed = bool(p.get("closed", False))
        s_ = np.linspace(0, 1, 80, endpoint=not closed)
        Q = np.array(splev(s_, p["tck"])).T
        if not closed:
            t0, t1 = unit(Q[1] - Q[0]), unit(Q[-1] - Q[-2])
            Q = np.vstack([Q[0] - t0 * margin, Q[1:-1], Q[-1] + t1 * margin])
        a, o, u, w = p["a"], p["o"], p["u"], p["w"]
        h0 = p["h0"] - margin
        base = o + (h0 - o @ a) * a
        pts = np.array([base + q[0] * u + q[1] * w for q in Q])
        return dict(kind="extrusion", points=pts, vec=a * ((p["h1"] + margin) - h0), closed=closed)
    if k == "revolution":
        from scipy.interpolate import splev
        c, a = p["c"], p["a"]
        Q = np.array(splev(np.linspace(0, 1, 80), p["tck"])).T
        t0, t1 = unit(Q[1] - Q[0]), unit(Q[-1] - Q[-2])
        start = Q[0] - t0 * margin if Q[0][1] > margin else Q[0]
        end = Q[-1] + t1 * margin if Q[-1][1] > margin else Q[-1]
        Q = np.vstack([start, Q[1:-1], end])
        Q[:, 1] = np.maximum(Q[:, 1], 1e-6)
        e, _ = perp_basis(a)
        pts = np.array([c + q[0] * a + q[1] * e for q in Q])
        return dict(kind="revolution", points=pts, c=c, a=a)
    if k == "facets":
        return dict(kind="facets", tris=np.asarray(p["tris"]))
    raise ValueError(k)


def make_tool(sp):
    k = sp["kind"]
    if k == "plane":
        return BH.plane_face(sp["n"], sp["o"], sp["u"], sp["rect"])
    if k == "cylinder":
        return BH.cylinder_face(sp["c"], sp["a"], sp["r"], sp["h"])
    if k == "cone":
        return BH.cone_face(sp["v"], sp["a"], sp["t"], sp["h"])
    if k == "sphere":
        return BH.sphere_face(sp["c"], sp["r"])
    if k == "torus":
        return BH.torus_face(sp["c"], sp["a"], sp["R"], sp["r"])
    if k == "extrusion":
        return BH.extrusion_face(sp["points"], sp["vec"], sp["closed"])
    if k == "profile":
        return BH.profile_face(sp["edges"], sp["vec"])
    if k == "revolution":
        return BH.revolution_face(sp["points"], sp["c"], sp["a"])
    if k == "facets":
        return BH.facet_faces(sp["tris"])
    raise ValueError(k)


COVER_PLANES = True


def tangent_chain_shells(prims, pairs, margin, log=print):
    """Build one extruded shell per chain of mutually tangent planes and cylinders.

    Planes tangent to cylinders (all parallel to a common axis) form a 2D G1 profile of
    lines and arcs; extruding that profile gives a gap-free tool surface, which a set of
    separately extended faces cannot guarantee at tangent contacts.
    Returns (shells, member_indices).
    """
    import networkx as nx

    # one tangency graph per cylinder-axis direction: a profile shell needs a common extrusion axis
    groups = []
    for (ci, pj) in pairs:
        a_ = prims[ci].p["a"]
        for g in groups:
            if abs(g[0] @ a_) > 0.9999:
                g[1].append((ci, pj))
                break
        else:
            groups.append((a_, [(ci, pj)]))
    comps = []
    for _, edges in groups:
        Gg = nx.Graph()
        Gg.add_edges_from(edges)
        comps += [(Gg, c) for c in nx.connected_components(Gg)]
    shells, members = [], set()
    for G, comp in comps:
        comp = sorted(comp)
        cyls = [i for i in comp if prims[i].kind == "cylinder"]
        if not cyls:
            continue
        a = prims[cyls[0]].p["a"]
        u, w = perp_basis(a)
        to2 = lambda X: np.c_[X @ u, X @ w] if X.ndim == 2 else np.array([X @ u, X @ w])
        # extents include the regions of removed fillets (extent_pts), so the shell reaches the
        # sharp edges the fillets will be re-applied on
        H = np.concatenate([getattr(prims[i], "extent_pts", prims[i].pts) @ a for i in comp])
        xm = max(getattr(prims[i], "extra_margin", 0.0) for i in comp)
        h0, h1 = H.min() - margin - xm, H.max() + margin + xm
        # order the chain: path or cycle in the bipartite tangency graph
        sub = G.subgraph(comp)
        if any(d > 2 for _, d in sub.degree()):
            continue
        ends = [n for n, d in sub.degree() if d == 1]
        order = list(nx.dfs_preorder_nodes(sub, ends[0] if ends else comp[0]))
        closed = not ends

        def tangent_point(ci, pj):
            c, n, d = prims[ci].p["c"], prims[pj].p["n"], prims[pj].p["d"]
            return to2(c - (n @ c - d) * n)

        def plane_far_end(pj, ci):
            # far end of a chain-end plane region along the line, away from its tangent point
            t = tangent_point(ci, pj)
            n2 = to2(prims[pj].p["n"])
            dvec = np.array([-n2[1], n2[0]])
            Q = to2(getattr(prims[pj], "extent_pts", prims[pj].pts))
            s = (Q - t) @ dvec
            sgn = 1 if s.max() > -s.min() else -1
            return t + dvec * sgn * (np.abs(s).max() + margin)

        def arc_points(ci, p_start, p_end_or_none):
            c2, r = to2(prims[ci].p["c"]), prims[ci].p["r"]
            Q = to2(prims[ci].pts) - c2
            mid_dir = np.mean(Q / np.linalg.norm(Q, axis=1, keepdims=True), 0)
            mid_dir /= np.linalg.norm(mid_dir)
            if p_end_or_none is None:
                # chain-end arc: from the tangent point across the region's angular span (found from the
                # largest gap between point angles, so spans over 180 deg work), plus a margin
                th = np.sort(np.arctan2(Q[:, 1], Q[:, 0]))
                gaps = np.diff(np.r_[th, th[0] + 2 * np.pi])
                k = int(np.argmax(gaps))
                t_start = th[(k + 1) % len(th)]
                span = 2 * np.pi - gaps[k]
                a0 = np.arctan2(*(p_start - c2)[::-1])
                d_start = abs(np.angle(np.exp(1j * (a0 - t_start))))
                d_end = abs(np.angle(np.exp(1j * (a0 - (t_start + span)))))
                ext = margin / r
                if d_start <= d_end:          # tangent point at the start: run counter-clockwise
                    a1 = a0 + span + ext
                else:                         # tangent point at the end: run clockwise
                    a1 = a0 - span - ext
                am = 0.5 * (a0 + a1)
                return p_start, c2 + r * np.array([np.cos(am), np.sin(am)]), c2 + r * np.array([np.cos(a1), np.sin(a1)])
            return p_start, c2 + r * mid_dir, p_end_or_none

        pts_edges = []          # list of ("line", P, Q) / ("arc", P, M, Q) in 2D
        covered = set()
        n = len(order)
        idx = list(range(n)) + ([0] if closed else [])
        for k in range(len(idx) - 1 if closed else n):
            i = order[idx[k]]
            prev = order[idx[k - 1]] if (k > 0 or closed) else None
            nxt = order[idx[k + 1]] if k + 1 < len(idx) else None
            if closed and k == 0:
                prev = order[-1]
            if prims[i].kind == "plane":
                A = tangent_point(prev, i) if prev is not None else plane_far_end(i, nxt)
                B = tangent_point(nxt, i) if nxt is not None else plane_far_end(i, prev)
                pts_edges.append(("line", A, B))
                # the line segment carries this plane in the shell; if it spans the whole region,
                # the plane's own face would coincide with the shell face and break the splitter
                L = np.linalg.norm(B - A)
                if L > 1e-12:
                    s = (to2(prims[i].pts) - A) @ ((B - A) / L)
                    if s.min() > -margin and s.max() < L + margin:
                        covered.add(i)
            else:
                A = tangent_point(i, prev) if prev is not None else None
                B = tangent_point(i, nxt) if nxt is not None else None
                if A is None:
                    P, M, Q = arc_points(i, B, None)
                    pts_edges.append(("arc", Q, M, P))
                elif B is None:
                    pts_edges.append(("arc", *arc_points(i, A, None)))
                else:
                    pts_edges.append(("arc", *arc_points(i, A, B)))
        to3 = lambda q: q[0] * u + q[1] * w + h0 * a
        edges = []
        for e in pts_edges:
            if e[0] == "line" and np.linalg.norm(e[1] - e[2]) < 1e-9:
                continue
            edges.append((e[0], *[to3(q) for q in e[1:]]))
        spec_ = dict(kind="profile", edges=edges, vec=a * (h1 - h0),
                     members=[(prims[i].kind, prims[i].p.get("r", prims[i].p.get("d"))) for i in order])
        try:
            make_tool(spec_)
        except Exception as ex:
            log(f"  tangent chain skipped ({ex})")
            continue
        shells.append(spec_)
        members.update(i for i in comp if prims[i].kind == "cylinder")
        members.update(covered if COVER_PLANES else ())
    if shells:
        log(f"  tangent chains: {len(shells)} profile shells replace {len(members)} surfaces")
    # members: cylinders of every chain, plus chain planes whose region the shell fully covers
    return shells, members


def junction_separators(mesh, prims, margin, angle_deg=20.0):
    """Triangles of thin ruled strips along every near-tangent junction between two regions.

    At a tangent (G1) junction the two extended tool faces touch instead of crossing, and a
    numerical gap of the order of the fit error lets inside and outside cells merge. A strip
    normal to both surfaces along the junction crosses both, which closes that channel; where it
    passes through material or air it only splits a cell, and inside cells are re-fused later."""
    import networkx as nx
    owner = np.full(len(mesh.faces), -1)
    for i, p in enumerate(prims):
        owner[p.faces] = i
    adj = mesh.face_adjacency
    a, b = owner[adj[:, 0]], owner[adj[:, 1]]
    sel = (a != b) & (a >= 0) & (b >= 0) & (mesh.face_adjacency_angles < np.radians(angle_deg))
    if not sel.any():
        return []
    E = mesh.face_adjacency_edges[sel]
    key = np.minimum(a[sel], b[sel]) * (len(prims) + 1) + np.maximum(a[sel], b[sel])
    V, VN = mesh.vertices, mesh.vertex_normals
    tris = []
    for k in np.unique(key):
        G = nx.Graph()
        G.add_edges_from(map(tuple, E[key == k]))
        for comp in nx.connected_components(G):
            H = G.subgraph(comp)
            ends = [v for v in H if H.degree(v) == 1]
            order = list(nx.dfs_preorder_nodes(H, ends[0] if ends else next(iter(H))))
            P = V[order]
            N = VN[order]
            if len(P) < 2:
                continue
            t0 = P[0] - P[1]
            t1 = P[-1] - P[-2]
            P = np.vstack([P[0] + margin * t0 / np.linalg.norm(t0), P, P[-1] + margin * t1 / np.linalg.norm(t1)])
            N = np.vstack([N[0], N, N[-1]])
            for i in range(len(P) - 1):
                q0, q1 = P[i] - margin * N[i], P[i] + margin * N[i]
                r0, r1 = P[i + 1] - margin * N[i + 1], P[i + 1] + margin * N[i + 1]
                tris.append((q0, q1, r1))
                tris.append((q0, r1, r0))
    return tris


def _inside_points(cell, rng, k=5, n_try=2000):
    """Up to k points strictly inside a cell, found on a coarse tessellation of the cell."""
    import trimesh
    bb = cell.BoundingBox()
    size = max(bb.xlen, bb.ylen, bb.zlen)
    try:
        verts, tris = cell.tessellate(1e-3 * size + 1e-9, 0.5)
        cm = trimesh.Trimesh(np.array([[p.x, p.y, p.z] for p in verts]), np.array(tris), process=True)
    except Exception:
        return np.zeros((0, 3))
    if len(cm.faces) < 4:
        return np.zeros((0, 3))
    lo, hi = cm.bounds
    if len(cm.faces) > 20000:            # large cells: fewer ray tests keep memory bounded
        n_try = 400
    Q = lo + rng.random((n_try, 3)) * (hi - lo)
    try:
        ins = cm.contains(Q)
    except Exception:
        ins = np.zeros(len(Q), bool)
    P = Q[ins]
    if len(P) >= 1:
        # prefer points far from the cell boundary (approximate distance to surface samples),
        # confirmed by the exact OCC classifier
        from scipy.spatial import cKDTree
        S, _ = trimesh.sample.sample_surface(cm, 3000, seed=0)
        d, _ = cKDTree(S).query(P[:400])
        P = P[:400][np.argsort(-d)]
    else:
        eps = 1e-3 * size                       # thin cells: step inward from face centroids
        P = cm.triangles_center - eps * cm.face_normals
    if cm.is_watertight and len(P) and len(Q[ins]):
        return P[:k]                     # watertight cell mesh: its inside test is reliable
    out = []
    for q in P[:4 * k]:                  # stop as soon as enough points are confirmed
        if cell.isInside(Vector(*q), 1e-9):
            out.append(q)
            if len(out) >= k:
                break
    tries = 0
    lo = np.array([bb.xmin, bb.ymin, bb.zmin])
    hi = np.array([bb.xmax, bb.ymax, bb.zmax])
    while len(out) < 1 and tries < 20:
        q = lo + rng.random(3) * (hi - lo)
        tries += 1
        if cell.isInside(Vector(*q), 1e-9):
            out.append(q)
    return np.array(out).reshape(-1, 3)


def shape_fix(shape):
    from OCP.ShapeFix import ShapeFix_Shape
    try:
        fx = ShapeFix_Shape(shape.wrapped)
        fx.Perform()
        fixed = Shape.cast(fx.Shape())
    except Exception:
        return shape
    ok = fixed.isValid() and fixed.Solids() and abs(fixed.Volume() - shape.Volume()) < 1e-3 * abs(shape.Volume())
    return fixed if ok else shape


def _flat(tools):
    out = []
    for t in tools:
        out += t if isinstance(t, list) else [t]
    return out


def _split_dropping(lo, hi, specs, tools, fuzzy, chunk=8):
    """Add tools in chunks; a chunk that breaks the split is retried tool by tool and the
    offending tools are dropped. Typical culprits are duplicate fits of one surface, whose
    faces nearly coincide."""
    keep, dropped, cells = [], [], []
    k = 0
    while k < len(tools):
        batch = list(range(k, min(k + chunk, len(tools))))
        c, ok = BH.split_box(lo, hi, _flat([tools[i] for i in keep + batch]), fuzzy, return_ok=True)
        if ok:
            keep += batch
            cells = c
        else:
            for i in batch:
                c, ok = BH.split_box(lo, hi, _flat([tools[j] for j in keep + [i]]), fuzzy, return_ok=True)
                if ok:
                    keep.append(i)
                    cells = c
                else:
                    dropped.append(specs[i])
        k += chunk
    if not keep:
        cells = BH.split_box(lo, hi, [], fuzzy)
    return [specs[i] for i in keep], [tools[i] for i in keep], cells, dropped


def build_arrangement(mesh, prims, margin, fuzzy, G=None, log=print, tangent_pairs=(), separators=False):
    """Returns (solid, info); info holds the tool specs, box and one interior point per kept cell,
    which is everything a generated script needs to rebuild the solid."""
    t0 = time.time()
    shells, skip = tangent_chain_shells(prims, tangent_pairs, margin, log) if tangent_pairs else ([], set())
    specs = list(shells)
    for k, pr in enumerate(prims):
        if k in skip:                            # covered by a chain shell (planes only if fully covered)
            continue
        sp_ = tool_spec(pr, margin)
        if sp_ is not None:
            specs.append(sp_)
    sep = junction_separators(mesh, prims, margin) if separators else []
    if sep:
        specs.append(dict(kind="facets", tris=np.array(sep), separator=True))
        log(f"  junction separators: {len(sep)} strip triangles at near-tangent junctions")
    tools = [make_tool(sp_) for sp_ in specs]
    lo, hi = mesh.bounds
    pad = 4 * margin
    box_lo, box_hi = lo - pad, hi + pad
    cells, ok = BH.split_box(box_lo, box_hi, _flat(tools), fuzzy, return_ok=True)
    if not ok:
        specs, tools, cells, dropped = _split_dropping(box_lo, box_hi, specs, tools, fuzzy)
        if dropped:
            log(f"  splitter failed; dropped {len(dropped)} tool(s) that break it: "
                + ", ".join(sp_['kind'] for sp_ in dropped))
    nt = len(_flat(tools))
    log(f"  arrangement: {nt} tool faces -> {len(cells)} cells ({time.time() - t0:.1f}s)")

    t1 = time.time()
    rng = np.random.default_rng(0)
    G = G or [np.eye(4)]
    pts, owner, per_cell = [], [], []
    for i, cell in enumerate(cells):
        P = _inside_points(cell, rng)
        per_cell.append(P)
        for T in G:                      # G[0] is the identity
            pts.append(P @ T[:3, :3].T + T[:3, 3])
            owner += [i] * len(P)
    pts, owner = np.vstack(pts), np.array(owner)
    if mesh.is_watertight:
        inside = mesh.contains(pts)
    else:
        from trimesh.proximity import signed_distance
        inside = signed_distance(mesh, pts) > 0
    frac = np.bincount(owner, weights=inside, minlength=len(cells)) / np.maximum(
        np.bincount(owner, minlength=len(cells)), 1)
    # ties follow the identity element (the unmapped original)
    ident_mask = np.zeros(len(pts), bool)
    k = 0
    for i, P in enumerate(per_cell):
        ident_mask[k:k + len(P)] = True
        k += len(P) * len(G)
    tie_in = np.bincount(owner[ident_mask], weights=inside[ident_mask], minlength=len(cells)) / np.maximum(
        np.bincount(owner[ident_mask], minlength=len(cells)), 1)
    inside_cells = [c for c, f, t in zip(cells, frac, tie_in) if f > 0.5 or (abs(f - 0.5) < 1e-9 and t > 0.5)]
    log(f"  classified cells: {len(inside_cells)} inside ({time.time() - t1:.1f}s)")
    if not inside_cells:
        raise RuntimeError("no cell classified as inside")
    t2 = time.time()
    solid = BH.fuse_cells(inside_cells)
    if not solid.isValid():
        solid = shape_fix(solid)
    sols = solid.Solids()
    if len(sols) > 1:                    # drop zero-volume slivers from the assembly
        big = [x for x in sols if x.Volume() > 1e-6 * max(y.Volume() for y in sols)]
        if len(big) < len(sols):
            from cadquery import func as cf_
            solid = big[0] if len(big) == 1 else cf_.compound(big)
            sols = big
    log(f"  fused and unified: {len(sols)} solid(s), {len(solid.Faces())} faces, "
        f"valid={solid.isValid()} ({time.time() - t2:.1f}s)")
    keep_pts = [per_cell[i][0] for i, (f, t) in enumerate(zip(frac, tie_in))
                if (f > 0.5 or (abs(f - 0.5) < 1e-9 and t > 0.5)) and len(per_cell[i])]
    info = dict(specs=specs, box=(box_lo, box_hi), points=np.array(keep_pts), fuzzy=fuzzy)
    return solid, info


def mark_holes(prims, mesh):
    """Flag full concave cylinders (drilled holes); the script generator writes them as hole features."""
    for p in prims:
        p.is_hole = False
        if p.kind != "cylinder" or len(p.faces) < 6:
            continue
        a, c = p.p["a"], p.p["c"]
        fc = mesh.triangles_center[p.faces]
        w = fc - c
        rad = w - np.outer(w @ a, a)
        rad /= np.linalg.norm(rad, axis=1, keepdims=True)
        concave = np.mean(np.einsum("ij,ij->i", rad, mesh.face_normals[p.faces]) < 0) > 0.9
        u, v = perp_basis(a)
        ang = np.arctan2(rad @ v, rad @ u)
        coverage = len(np.unique(np.floor((ang + np.pi) / (2 * np.pi) * 36))) / 36
        p.is_hole = bool(concave and coverage > 0.9)


def reconstruct(mesh, opts, spec=None, log=print):
    """opts: dict with tol, kinds, allow_extrusion, allow_facets, snap_deg, angle_step, tangency,
    grid, margin, fuzzy."""
    tol = opts["tol"]
    log(f"exact mode: tolerance {tol:.4f} mm, surfaces {', '.join(opts['kinds'])}"
        f"{', spline extrusions' if opts.get('allow_extrusion') else ''}")
    if opts.get("prism25", True):
        from . import prism25
        res = prism25.reconstruct(mesh, tol, log=log)
        if res is not None:
            return res
        from . import revolve
        res = revolve.reconstruct(mesh, tol, log=log)    # a body of revolution: one revolved profile
        if res is not None:
            return res
        res = revolve.reconstruct_with_cuts(mesh, tol, log=log)   # with straight cuts across it
        if res is not None:
            return res
    t = time.time()
    seg = Segmenter(mesh, tol, kinds=opts["kinds"], allow_extrusion=opts.get("allow_extrusion", False),
                    allow_facets=opts.get("allow_facets", True), log=log)
    prims = seg.run()
    s, analytic = summarize(prims, mesh)
    log(f"  segmentation ({time.time() - t:.1f}s): {s}")
    if not opts.get("allow_facets", True) and any(p.kind == "facets" for p in prims):
        n = sum(p.kind == "facets" for p in prims)
        log(f"  warning: {n} region(s) could not be fitted with the allowed surfaces; "
            f"they are left as planar facets")
    from .symmetry import principal_frame
    axes = spec.axes if spec is not None else principal_frame(mesh)[1]
    prims = simplify_degenerate(prims, tol, log=log)
    prims = facets_to_cylinders(prims, mesh, tol, log=log)
    prims = spheres_to_cylinders(prims, mesh, tol, log=log)
    prims = decompose_revolutions(prims, mesh, tol, log=log)
    if "torus" in opts["kinds"]:
        prims = cones_to_tori(prims, mesh, tol, log=log)
    hint = opts.get("symmetry_hint")
    Gh = None
    if hint is not None:
        from .symmetry import group_matrices
        Gh = [T for T in group_matrices(hint) if not np.allclose(T, np.eye(4))]
    for _ in range(2):                                # consensus and refit help each other
        prims = symmetric_consensus(prims, tol, Gh, log=log)
        prims = refit_odd_regions(prims, mesh, tol, log=log, group=Gh)
        if any(p.kind == "torus" for p in prims):
            from .segment import _grow_tori
            prims = _grow_tori(prims, mesh, tol)      # seam faces go to the torus their normals follow
    prims = symmetric_consensus(prims, tol, Gh, log=log)
    if opts.get("snap_deg", 1.0) > 0:
        prims = snap_directions(prims, axes, opts.get("snap_deg", 1.0), opts.get("angle_step", 15.0), log=log, tol=tol)
    prims = unify_parallel(prims, tol, log=log)
    pairs = []
    if opts.get("tangency", True):
        prims, pairs = snap_tangency(prims, mesh, tol, log=log)
        from .segment import drop_seam_bridges
        prims, pairs = drop_seam_bridges(prims, mesh, pairs, log=log)
    if opts.get("grid", 0):
        prims = snap_grid(prims, mesh.bounds[0], axes, opts["grid"], log=log)
    G = None
    impose = opts.get("impose_symmetry", True)
    if spec is not None and len(group_matrices(spec)) > 1:
        Gs = group_matrices(spec)
        prims = symmetrize(prims, Gs, tol, float(np.linalg.norm(mesh.extents)), log=log, complete=impose)
        G = Gs if impose else None
        pairs = []                       # indices change under symmetrization; rebuild tangency
        prims = reassign_faces(prims, mesh, tol)
        if opts.get("tangency", True):
            prims, pairs = snap_tangency(prims, mesh, tol, log=log)
    mark_holes(prims, mesh)
    blends = []
    if opts.get("fillets", True) and any(p.kind in ("cylinder", "torus") for p in prims):
        from .fillets import detect_blends, remove_blends
        # plane-tangent cylinders stay with the tangent-chain builder
        chain_cyl = {i for pr_ in pairs for i in pr_ if prims[i].kind == "cylinder"}
        blends = detect_blends(prims, mesh, tol, exclude=chain_cyl, log=lambda *a: None)
        if any(b.kind == "torus" for b in blends):
            from .fillets import network_cylinders
            net = network_cylinders(prims, mesh, blends, tol, chain_cyl)
            if net:
                chain_cyl -= net
                pairs = [pr_ for pr_ in pairs if not (set(pr_) & net)]
        blends = detect_blends(prims, mesh, tol, exclude=chain_cyl, log=log)
        if blends:
            prims, pairs, blends = remove_blends(prims, blends, pairs)
    margin = opts.get("margin") or max(5 * tol, 0.004 * float(np.linalg.norm(mesh.extents)))
    fuzzies = [opts["fuzzy"]] if opts.get("fuzzy") else [0.0]
    solid, best, info = None, None, None
    for fz in fuzzies:
        try:
            cand, cinfo = build_arrangement(mesh, prims, margin, fz, G=G, log=log, tangent_pairs=pairs,
                                             separators=opts.get("separators", False))
        except Exception as e:
            log(f"  arrangement with fuzzy {fz:.2e} failed: {e}")
            continue
        vol_err = abs(cand.Volume() / mesh.volume - 1) if cand.Solids() else 1.0
        good = cand.isValid() and len(cand.Solids()) == mesh.body_count and vol_err < 0.01
        if best is None or vol_err < best:
            solid, best, info = cand, vol_err, cinfo
        if good:
            break
        log(f"  result check: valid={cand.isValid()}, solids={len(cand.Solids())} "
            f"(mesh bodies {mesh.body_count}), volume error {100 * vol_err:.2f} %")
    if solid is None:
        raise RuntimeError("exact reconstruction failed")
    if blends:
        from .fillets import apply_fillets
        t1 = time.time()
        try:
            frec = {}
            solid, n_done = apply_fillets(solid, prims, blends, tol, log=log, record=frec)
            if info is not None:
                info["fillets"] = frec.get("steps", [])
                info["tori"] = frec.get("tori", [])
        except Exception as e:
            log(f"  fillets failed: {e}")
        best = abs(solid.Volume() / mesh.volume - 1) if solid.Solids() else 1.0
        log(f"  after fillets: valid={solid.isValid()}, volume error {100 * best:.2f} % ({time.time() - t1:.1f}s)")
    res = Result(shape=solid, mode="exact")
    res.extras = dict(prims=prims, analytic_fraction=analytic, summary=s, volume_error=best, brep=info,
                      valid=solid.isValid(), passed=bool(best is not None and best < 0.01))
    return res
