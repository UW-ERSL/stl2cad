"""Relations and dimensions that fully define a closed sketch outline, without redundancy.

Unknowns: the corners of the outline and the centres of its arcs (x, y each). Every relation or
dimension is an equation in them. A candidate is kept only if it raises the rank of the system's
Jacobian (evaluated at the actual geometry), i.e. if it removes a degree of freedom; the outline is
fully defined when the rank equals the number of unknowns. Candidates, in order: horizontal and
vertical relations, tangency at line-arc corners, the corner on the origin (coincident) or its
distances from the origin, arc radii, lengths of horizontal and vertical edges, horizontal and
vertical distances between the ends of slanted edges, and as a last resort corner distances from
the origin.
"""
import math

import numpy as np


def _centre(p, q, r):
    ax, ay = p
    bx, by = q
    cx, cy = r
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    return np.array([ux, uy])


def _tangent(line_dir, corner, centre):
    rad = corner - centre
    c = abs(line_dir @ rad) / (np.linalg.norm(line_dir) * np.linalg.norm(rad) + 1e-30)
    return c < math.sin(math.radians(0.5))


# singular values below RANK_TOL (relative) count as zero: numerical noise in the constraint
# Jacobians sits near 1e-7 (finite differences, nearly symmetric data), real constraints at 1e-3 and up
RANK_TOL = 1e-5


def plan_outline(edges, tag="", linear=None, rotational=None, pitch=None, copies=None):
    """edges: ("line", a, b) / ("arc", a, mid, b) in (a, b) sketch coordinates (mm).
    linear: (p, N) when edges[:p N] are N copies of the first p edges made by a sketch linear
    pattern: the copies then follow the period (each corner moved by the period's own end-to-start
    offset), and no operation picks a copy, which the script did not draw.
    rotational: (p, N, centre) when the edges are N copies of the first p turned about centre by a
    sketch circular pattern; handled the same way.
    copies: (p, [(start, k), ...]) when edges start .. start + p - 1 are made by a sketch linear
    pattern of the first p edges (the seed), moved k times the seed's own end-to-start offset;
    the pattern's spacing and direction fix that offset, and no operation picks a copy.
    pitch: p when the first edges are drawn copies of a period of p edges (tied by Equal and
    Parallel relations): the pitch, from corner 0 to corner p, is dimensioned first, so that a
    change of a feature's size keeps the pitch.
    Returns (ops, fully_defined); ops as executed by the SolidWorks script."""
    n = len(edges)
    p_, N_ = linear if linear is not None else (rotational[:2] if rotational is not None else (0, 0))
    made = set(range(p_, p_ * N_))                                      # edges made by a pattern
    if copies is not None:
        made = {(c + j) % n for c, _ in copies[1] for j in range(copies[0])}
    ok_e = lambda i: i % n not in made                                  # edge drawn by the script
    ok_c = lambda i: i % n not in made or (i - 1) % n not in made       # corner on a drawn edge
    V = np.array([np.asarray(e[1], float) for e in edges])            # corner i = start of edge i
    arcs = [i for i, e in enumerate(edges) if e[0] == "arc"]
    C = {i: _centre(*[np.asarray(q, float) for q in edges[i][1:]]) for i in arcs}
    nv = 2 * n + 2 * len(arcs)
    x0 = np.concatenate([V.ravel()] + [C[i] for i in arcs])
    cidx = {i: 2 * n + 2 * k for k, i in enumerate(arcs)}

    def corner(x, i):
        return x[2 * (i % n): 2 * (i % n) + 2]

    def centre(x, i):
        return x[cidx[i]: cidx[i] + 2]

    def jac(fs):
        rows = []
        h = 1e-4                # central differences; 1e-6 left rounding noise near the rank tolerance
        for f in fs:
            g = np.zeros(nv)
            for k in range(nv):
                xp, xm = x0.copy(), x0.copy()
                xp[k] += h
                xm[k] -= h
                g[k] = (f(xp) - f(xm)) / (2 * h)
            rows.append(g)
        return np.array(rows)

    rows = np.zeros((0, nv))
    rank = 0
    ops = []

    def try_add(fs, op):
        nonlocal rows, rank
        R = np.vstack([rows, jac(fs)])
        r = np.linalg.matrix_rank(R, tol=RANK_TOL * max(1.0, np.abs(R).max()))
        if r > rank:
            rows, rank = R, r
            ops.append(op)
            return True
        return False

    mid = lambda e: tuple(((np.asarray(e[1]) + np.asarray(e[2])) / 2).tolist()) if e[0] == "line" else tuple(map(float, e[2]))
    P = lambda v: (float(v[0]), float(v[1]))
    kinds = []
    for e in edges:
        if e[0] == "arc":
            kinds.append("A")
            continue
        d = np.asarray(e[2]) - np.asarray(e[1])
        kinds.append("V" if abs(d[0]) < 1e-9 else "H" if abs(d[1]) < 1e-9 else "S")
    # arcs: start and end on one circle (part of the sketch entity itself, not an operation)
    for i in arcs:
        f = lambda x, i=i: np.sum((corner(x, i) - centre(x, i)) ** 2) - np.sum((corner(x, i + 1) - centre(x, i)) ** 2)
        R = np.vstack([rows, jac([f])])
        rows, rank = R, np.linalg.matrix_rank(R, tol=RANK_TOL * max(1.0, np.abs(R).max()))
    # copies of a sketch pattern follow the period
    if linear is not None:
        fs = [lambda x, m=m, j=j, a=a: (corner(x, m * p_ + j) - corner(x, j) - m * (corner(x, p_) - corner(x, 0)))[a]
              for m in range(1, N_) for j in range(p_) for a in (0, 1)]
        R = np.vstack([rows, jac(fs)])
        rows, rank = R, np.linalg.matrix_rank(R, tol=RANK_TOL * max(1.0, np.abs(R).max()))
    if copies is not None:          # pattern copies, and the pattern's spacing and direction
        q_ = copies[0]
        fs = [lambda x, c=c, k=k, j=j, a=a: (corner(x, c + j) - corner(x, j) - k * (corner(x, q_) - corner(x, 0)))[a]
              for c, k in copies[1] for j in range(q_) for a in (0, 1)]
        fs += [lambda x, a=a: (corner(x, q_) - corner(x, 0))[a] for a in (0, 1)]
        R = np.vstack([rows, jac(fs)])
        rows, rank = R, np.linalg.matrix_rank(R, tol=RANK_TOL * max(1.0, np.abs(R).max()))
    if rotational is not None:
        c0 = np.asarray(rotational[2], float)

        def turned(x, j, m):
            t = 2 * np.pi * m / N_
            q = corner(x, j) - c0
            return c0 + np.array([np.cos(t) * q[0] - np.sin(t) * q[1], np.sin(t) * q[0] + np.cos(t) * q[1]])
        fs = [lambda x, m=m, j=j, a=a: (corner(x, m * p_ + j) - turned(x, j, m))[a]
              for m in range(1, N_) for j in range(p_) for a in (0, 1)]
        R = np.vstack([rows, jac(fs)])
        rows, rank = R, np.linalg.matrix_rank(R, tol=RANK_TOL * max(1.0, np.abs(R).max()))
    # mirror symmetry about the sketch axes through the origin: a construction centerline per
    # axis, mirrored corners tied with Symmetric relations, corners on the axis put on it
    ext = float(np.abs(V).max()) if n else 1.0
    for ax in (0, 1):
        M = np.array([1.0, 1.0])
        M[ax] = -1.0
        part = {}
        for i in range(n):
            d = np.linalg.norm(V - V[i] * M, axis=1)
            k = int(np.argmin(d))
            if d[k] > 1e-6 * max(1.0, ext):
                break
            part[i] = k
        else:
            mids = [np.asarray(edges[i][2], float) for i in arcs]
            if not all(any(np.linalg.norm(m * M - q) < 1e-6 * max(1.0, ext) for q in mids) for m in mids):
                continue
            Lc = 1.3 * ext + 1.0
            a_ = np.zeros(2)
            a_[1 - ax] = Lc                                  # along the axis line
            cl = (tuple((-a_).tolist()), tuple(a_.tolist()))
            pick = P(0.9 * a_)
            ops.append(("centerline", ax, cl, f"{tag}centerline {'a' if ax == 0 else 'b'} = 0"))
            for i, k in sorted(part.items()):
                if not (ok_c(i) and ok_c(k)) or k < i:
                    continue
                if k == i:
                    try_add([lambda x, i=i, ax=ax: corner(x, i)[ax]],
                            ("rel", "sgCOINCIDENT", [("SKETCHPOINT", P(V[i])), ("SKETCHSEGMENT", pick)],
                             f"{tag}corner {P(V[i])} on the centerline"))
                else:
                    try_add([lambda x, i=i, k=k, M=M, c=c_: (corner(x, k) - corner(x, i) * M)[c] for c_ in (0, 1)],
                            ("rel", "sgSYMMETRIC", [("SKETCHPOINT", P(V[i])), ("SKETCHPOINT", P(V[k])),
                                                    ("SKETCHSEGMENT", pick)],
                             f"{tag}corners {P(V[i])} and {P(V[k])} symmetric"))
    # relations
    for i, k in enumerate(kinds):
        if not ok_e(i):
            continue
        if k == "H":
            try_add([lambda x, i=i: corner(x, i + 1)[1] - corner(x, i)[1]],
                    ("rel", "sgHORIZONTAL2D", [("SKETCHSEGMENT", mid(edges[i]))], f"{tag}horizontal edge at {P(mid(edges[i]))}"))
        elif k == "V":
            try_add([lambda x, i=i: corner(x, i + 1)[0] - corner(x, i)[0]],
                    ("rel", "sgVERTICAL2D", [("SKETCHSEGMENT", mid(edges[i]))], f"{tag}vertical edge at {P(mid(edges[i]))}"))
    for i in range(n):                                          # tangency at corner i
        p, q = (i - 1) % n, i
        if not (ok_e(p) and ok_e(q)):
            continue
        kp, kq = kinds[p], kinds[q]
        if {kp, kq} != {"A"} and "A" in (kp, kq) and kp != kq:
            li, ai = (p, q) if kp != "A" else (q, p)
            other = (li + 1) if li == p else li + 1             # far end of the line from corner i
            far = V[li % n] if li == p else V[(li + 1) % n]
            d = V[i] - far
            if _tangent(d, V[i], C[ai]):
                if li == p:
                    f = lambda x, i=i, li=li, ai=ai: (corner(x, i) - corner(x, li)) @ (corner(x, i) - centre(x, ai))
                else:
                    f = lambda x, i=i, li=li, ai=ai: (corner(x, i) - corner(x, li + 1)) @ (corner(x, i) - centre(x, ai))
                try_add([f], ("rel", "sgTANGENT", [("SKETCHSEGMENT", mid(edges[li % n])), ("SKETCHSEGMENT", mid(edges[ai]))],
                              f"{tag}tangency at {P(V[i])}"))
    # equal radii and equal lengths: one dimension drives each group (Equal relations to the
    # group's first member, which is dimensioned below)
    def group(vals):
        gs = []
        for i, v in vals:
            for g in gs:
                if abs(g[0][1] - v) < 1e-6 * max(1.0, abs(v)):
                    g.append((i, v))
                    break
            else:
                gs.append([(i, v)])
        return [g for g in gs if len(g) > 1]
    radius = lambda x, i: np.linalg.norm(corner(x, i) - centre(x, i))
    length = lambda x, i: np.linalg.norm(corner(x, i + 1) - corner(x, i))
    # each member is tied to the previous one (a chain), not all to the first, so the relation
    # icons spread over the sketch instead of piling up on one edge
    for g in group([(i, float(np.linalg.norm(V[i] - C[i]))) for i in arcs if ok_e(i)]):
        for (i0, _), (i, _) in zip(g[:-1], g[1:]):
            try_add([lambda x, i=i, i0=i0: radius(x, i) - radius(x, i0)],
                    ("rel", "sgSAMELENGTH", [("SKETCHSEGMENT", mid(edges[i0])), ("SKETCHSEGMENT", mid(edges[i]))],
                     f"{tag}equal radius, arcs at {P(mid(edges[i0]))} and {P(mid(edges[i]))}"))
    for g in group([(i, float(np.linalg.norm(V[(i + 1) % n] - V[i]))) for i in range(n)
                    if kinds[i] != "A" and ok_e(i)]):
        for (i0, _), (i, _) in zip(g[:-1], g[1:]):
            try_add([lambda x, i=i, i0=i0: length(x, i) - length(x, i0)],
                    ("rel", "sgSAMELENGTH", [("SKETCHSEGMENT", mid(edges[i0])), ("SKETCHSEGMENT", mid(edges[i]))],
                     f"{tag}equal length, edges at {P(mid(edges[i0]))} and {P(mid(edges[i]))}"))
    # the pitch of drawn copies, before any length: changing a tooth then changes the gaps
    if pitch:
        a, b = V[0], V[pitch % n]
        d = b - a
        kind = "H" if abs(d[1]) < 1e-9 else ("V" if abs(d[0]) < 1e-9 else "any")
        text = ((a[0] + b[0]) / 2, min(a[1], b[1]) - 20) if kind != "V" else (min(a[0], b[0]) - 20, (a[1] + b[1]) / 2)
        fs = ([lambda x, ax=(0 if kind == "H" else 1): corner(x, pitch)[ax] - corner(x, 0)[ax]] if kind != "any"
              else [lambda x: np.linalg.norm(corner(x, pitch) - corner(x, 0))])
        pitch_op = ("dim", kind, [("SKETCHPOINT", P(a)), ("SKETCHPOINT", P(b))], P(text), f"{tag}pitch {P(a)} to {P(b)}")
    # slanted lines of one direction: Parallel to the first of their group (with Equal lengths
    # this ties repeated features, such as the teeth of a rack, to the first one)
    dirs = []
    for i in range(n):
        if kinds[i] != "S" or not ok_e(i):
            continue
        d = V[(i + 1) % n] - V[i]
        d = d / np.linalg.norm(d)
        for g in dirs:
            if abs(g[0][0] * d[1] - g[0][1] * d[0]) < 1e-6:
                g[1].append(i)
                break
        else:
            dirs.append((d, [i]))
    for _, idx in dirs:
        for i0, i in zip(idx[:-1], idx[1:]):
            try_add([lambda x, i=i, i0=i0: ((corner(x, i + 1) - corner(x, i))[0] * (corner(x, i0 + 1) - corner(x, i0))[1]
                                            - (corner(x, i + 1) - corner(x, i))[1] * (corner(x, i0 + 1) - corner(x, i0))[0])],
                    ("rel", "sgPARALLEL", [("SKETCHSEGMENT", mid(edges[i0])), ("SKETCHSEGMENT", mid(edges[i]))],
                     f"{tag}parallel, edges at {P(mid(edges[i0]))} and {P(mid(edges[i]))}"))
    if pitch:
        try_add(fs, pitch_op)
    # position: corner nearest the origin
    j = min((i for i in range(n) if ok_c(i)), key=lambda i: float(np.linalg.norm(V[i])))
    if np.linalg.norm(V[j]) < 1e-9:
        try_add([lambda x, j=j: corner(x, j)[0], lambda x, j=j: corner(x, j)[1]],
                ("rel", "sgCOINCIDENT", [("SKETCHPOINT", P(V[j])), ("ORIGIN", None)], f"{tag}corner on the origin"))
    else:
        for axis, kind in ((0, "H"), (1, "V")):
            if abs(V[j][axis]) > 1e-9:
                text = (V[j][0] / 2, V[j][1] - 8) if kind == "H" else (V[j][0] - 8, V[j][1] / 2)
                try_add([lambda x, j=j, axis=axis: corner(x, j)[axis]],
                        ("dim", kind, [("SKETCHPOINT", P(V[j])), ("ORIGIN", None)], text,
                         f"{tag}{'x' if axis == 0 else 'y'} of corner {P(V[j])}"))
    area = sum(V[i][0] * V[(i + 1) % n][1] - V[(i + 1) % n][0] * V[i][1] for i in range(n))
    out = 1.0 if area > 0 else -1.0
    # arc radii
    for i in [i for i in arcs if ok_e(i)]:
        c, m = C[i], np.asarray(edges[i][2], float)
        try_add([lambda x, i=i: np.linalg.norm(corner(x, i) - centre(x, i))],
                ("dim", "any", [("SKETCHSEGMENT", P(m))], P(m + (m - c)), f"{tag}radius of arc at {P(m)}"))
    stagger = {"H": 0, "V": 0}                  # successive extent dimensions step outwards

    def slanted_extents(which):
        for i in which:
            if kinds[i] == "S" and ok_e(i):
                a, b = V[i], V[(i + 1) % n]
                for axis, kind in ((0, "H"), (1, "V")):
                    off = 6 + 3 * (stagger[kind] % 4)
                    text = ((a[0] + b[0]) / 2, min(a[1], b[1]) - off) if kind == "H" else (min(a[0], b[0]) - off, (a[1] + b[1]) / 2)
                    stagger[kind] += 1
                    try_add([lambda x, i=i, axis=axis: corner(x, i + 1)[axis] - corner(x, i)[axis]],
                            ("dim", kind, [("SKETCHPOINT", P(a)), ("SKETCHPOINT", P(b))], P(text),
                             f"{tag}{'horizontal' if kind == 'H' else 'vertical'} extent of slanted edge at {P(mid(edges[i]))}"))
    if pitch or copies:          # the period's flanks first: a feature's shape, not the gap, is dimensioned
        slanted_extents(range(pitch or copies[0]))
    # lengths of horizontal and vertical edges
    for i, k in enumerate(kinds):
        if k in "HV" and ok_e(i):
            a, b = V[i], V[(i + 1) % n]
            d = b - a
            L = float(np.linalg.norm(d))
            m = np.asarray(mid(edges[i]))
            t = m + 6 * out * np.array([d[1], -d[0]]) / L
            axis = 0 if k == "H" else 1
            try_add([lambda x, i=i, axis=axis: corner(x, i + 1)[axis] - corner(x, i)[axis]],
                    ("dim", "any", [("SKETCHSEGMENT", P(m))], P(t), f"{tag}length of edge at {P(m)}"))
    # slanted edges: horizontal and vertical distance between their ends
    slanted_extents(range(n))
    # last resort: corner positions from the origin
    for i in [i for i in range(n) if ok_c(i)]:
        for axis, kind in ((0, "H"), (1, "V")):
            if rank >= nv:
                break
            if abs(V[i][axis]) < 1e-9:
                continue
            text = (V[i][0] / 2, V[i][1] - 8) if kind == "H" else (V[i][0] - 8, V[i][1] / 2)
            try_add([lambda x, i=i, axis=axis: corner(x, i)[axis]],
                    ("dim", kind, [("SKETCHPOINT", P(V[i])), ("ORIGIN", None)], P(text),
                     f"{tag}{'x' if axis == 0 else 'y'} of corner {P(V[i])}"))
    return ops, rank >= nv
