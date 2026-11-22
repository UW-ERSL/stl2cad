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


def plan_outline(edges, tag=""):
    """edges: ("line", a, b) / ("arc", a, mid, b) in (a, b) sketch coordinates (mm).
    Returns (ops, fully_defined); ops as executed by the SolidWorks script."""
    n = len(edges)
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
        h = 1e-6
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
        r = np.linalg.matrix_rank(R, tol=1e-7 * max(1.0, np.abs(R).max()))
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
        rows, rank = R, np.linalg.matrix_rank(R, tol=1e-7 * max(1.0, np.abs(R).max()))
    # relations
    for i, k in enumerate(kinds):
        if k == "H":
            try_add([lambda x, i=i: corner(x, i + 1)[1] - corner(x, i)[1]],
                    ("rel", "sgHORIZONTAL2D", [("SKETCHSEGMENT", mid(edges[i]))], f"{tag}horizontal edge at {P(mid(edges[i]))}"))
        elif k == "V":
            try_add([lambda x, i=i: corner(x, i + 1)[0] - corner(x, i)[0]],
                    ("rel", "sgVERTICAL2D", [("SKETCHSEGMENT", mid(edges[i]))], f"{tag}vertical edge at {P(mid(edges[i]))}"))
    for i in range(n):                                          # tangency at corner i
        p, q = i - 1, i
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
    # position: corner nearest the origin
    j = int(np.argmin(np.linalg.norm(V, axis=1)))
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
    for i in arcs:
        c, m = C[i], np.asarray(edges[i][2], float)
        try_add([lambda x, i=i: np.linalg.norm(corner(x, i) - centre(x, i))],
                ("dim", "any", [("SKETCHSEGMENT", P(m))], P(m + (m - c)), f"{tag}radius of arc at {P(m)}"))
    # lengths of horizontal and vertical edges
    for i, k in enumerate(kinds):
        if k in "HV":
            a, b = V[i], V[(i + 1) % n]
            d = b - a
            L = float(np.linalg.norm(d))
            m = np.asarray(mid(edges[i]))
            t = m + 6 * out * np.array([d[1], -d[0]]) / L
            axis = 0 if k == "H" else 1
            try_add([lambda x, i=i, axis=axis: corner(x, i + 1)[axis] - corner(x, i)[axis]],
                    ("dim", "any", [("SKETCHSEGMENT", P(m))], P(t), f"{tag}length of edge at {P(m)}"))
    # slanted edges: horizontal and vertical distance between their ends
    for i, k in enumerate(kinds):
        if k == "S":
            a, b = V[i], V[(i + 1) % n]
            for axis, kind in ((0, "H"), (1, "V")):
                text = ((a[0] + b[0]) / 2, min(a[1], b[1]) - 6) if kind == "H" else (min(a[0], b[0]) - 6, (a[1] + b[1]) / 2)
                try_add([lambda x, i=i, axis=axis: corner(x, i + 1)[axis] - corner(x, i)[axis]],
                        ("dim", kind, [("SKETCHPOINT", P(a)), ("SKETCHPOINT", P(b))], P(text),
                         f"{tag}{'horizontal' if kind == 'H' else 'vertical'} extent of slanted edge at {P(mid(edges[i]))}"))
    # last resort: corner positions from the origin
    for i in range(n):
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
