"""Write a CadQuery script for a TO result: the original features, then one plate per detected plate.

plate = original part cut to the plate's slab [lo, hi] along its axis, minus its cuts;
the result is the union of the plates. Each cut is sketched in the plate's (u, v) plane.
"""
import math

HELPERS = '''
# ---- topology-optimized part: plates -----------------------------------------
# Each plate is the original part cut to a slab [lo, hi] along its axis, minus cuts
# sketched in the plate plane (u, v) = the other two axes in increasing order.
# A cut is a closed polygon whose vertices carry tangent fillet radii (x, y, r)
# (r = 0: sharp), or an arc notch closed by a polyline outside the original outline.
original = part


def to3(axis, u, v, w):
    """Point (u, v) of a plate plane at level w along axis, in model coordinates."""
    p = [0.0, 0.0, 0.0]
    a, b = [i for i in range(3) if i != axis]
    p[a], p[b], p[axis] = u, v, w
    return Vector(*p)


def rounded_polygon(verts, axis, w):
    """Closed polygon with a tangent fillet of radius r at each vertex (r = 0: sharp)."""
    V = [Vector(x, y, 0) for x, y, _ in verts]
    R = [r for _, _, r in verts]
    n = len(V)
    ends = []
    for k in range(n):
        a, b, c = V[k - 1], V[k], V[(k + 1) % n]
        u = (b - a).normalized()
        q = (c - b).normalized()
        if R[k] <= 0:
            ends.append((b, b, None))
            continue
        th = math.acos(max(-1.0, min(1.0, u.dot(q))))
        t = R[k] * math.tan(th / 2)
        bis = (q - u).normalized()
        centre = b + bis * (R[k] / math.cos(th / 2))
        ends.append((b - u * t, b + q * t, centre - bis * R[k]))
    P = lambda p: to3(axis, p.x, p.y, w)
    edges = []
    for k in range(n):
        t1, t2, mid = ends[k]
        if mid is not None:
            edges.append(Edge.makeThreePointArc(P(t1), P(mid), P(t2)))
        nxt = ends[(k + 1) % n][0]
        if (nxt - t2).Length > 1e-6:
            edges.append(cf.segment(P(t2), P(nxt)))
    return cf.wire(*edges)


def arc_notch(cx, cy, r, a0, a1, path, axis, w):
    pt = lambda a: to3(axis, cx + r * math.cos(a), cy + r * math.sin(a), w)
    edges = [Edge.makeThreePointArc(pt(a0), pt((a0 + a1) / 2), pt(a1))]
    poly = [pt(a1)] + [to3(axis, x, y, w) for x, y in path] + [pt(a0)]
    edges += [cf.segment(p, q) for p, q in zip(poly[:-1], poly[1:])]
    return cf.wire(*edges)


def slab(axis, lo, hi, box):
    """Box spanning the original part, cut to [lo, hi] along axis."""
    (x0, y0, z0), (x1, y1, z1) = box
    lo_, hi_ = [x0 - 1, y0 - 1, z0 - 1], [x1 + 1, y1 + 1, z1 + 1]
    lo_[axis], hi_[axis] = lo, hi
    return Solid.makeBox(hi_[0] - lo_[0], hi_[1] - lo_[1], hi_[2] - lo_[2], Vector(*lo_))


def through(wire, axis, lo, hi):
    d = [0.0, 0.0, 0.0]
    d[axis] = hi - lo + 2
    return cf.extrude(cf.face(wire), d)


def concave_edges(solid, angle=0.01):
    """Re-entrant edges, from the two face normals at each edge (OCC BRepOffset_Analyse:
    convex, concave or tangent from the sign of (n1 x n2) . t along the edge)."""
    an = BRepOffset_Analyse(solid.wrapped, angle)
    return [e for e in solid.Edges()
            if any(i.Type() == ChFiDS_TypeOfConcavity.ChFiDS_Concave for i in an.Type(e.wrapped))]


def fillet_reentrant(solid, r):
    """Fillet every re-entrant edge with radius r (all at once; if that fails, one edge at a time)."""
    edges = concave_edges(solid)
    if not edges:
        return solid, 0, 0
    try:
        return solid.fillet(r, edges), len(edges), 0
    except Exception:
        pass
    done = skipped = 0
    for e0 in edges:
        mid = Vertex.makeVertex(*e0.positionAt(0.5).toTuple())
        cand = [e for e in concave_edges(solid) if e.distance(mid) < 1e-3]
        if not cand:
            continue
        try:
            solid = solid.fillet(r, cand[:1]); done += 1
        except Exception:
            skipped += 1
    return solid, done, skipped


def plate(axis, lo, hi, cuts, box):
    body = original.intersect(slab(axis, lo, hi, box))
    for cut in cuts:
        if cut[0] == "poly":
            body = body - through(rounded_polygon(cut[1], axis, lo - 1), axis, lo, hi)
        else:
            (cx, cy, r, a0, a1), path = cut[1], cut[2]
            body = body - through(arc_notch(cx, cy, r, a0, a1, path, axis, lo - 1), axis, lo, hi)
    return body

'''


def fmt_cut(c):
    if c[0] == 'poly':
        rows = ''.join('                (%9.4f, %9.4f, %7.3f),\n' % v for v in c[1])
        return '            ("poly", [\n%s            ]),\n' % rows
    rows = ''.join('                (%9.4f, %9.4f),\n' % p for p in c[2])
    return '            ("arc", (%.4f, %.4f, %.4f, %.6f, %.6f), [\n%s            ]),\n' % (*c[1], rows)


def write(path, base_script, name, plates, box, r_fillet=1.0):
    src = open(base_script).read()
    head, tail = src.split('part = part.clean()')
    head = head.replace('from cadquery import Edge, Vector', 'from cadquery import Edge, Solid, Vector, Vertex\nfrom OCP.BRepOffset import BRepOffset_Analyse\nfrom OCP.ChFiDS import ChFiDS_TypeOfConcavity')
    doc = head.split('"""')[1]
    head = head.replace(doc, '%s: topology-optimized part as plates of the original part.\n\n'
                        'Original features first, then %d plate(s); units mm.\nRun:  python %s.py [out.step]\n'
                        % (name, len(plates), name), 1)
    body = HELPERS + 'BOX = (%s, %s)\n\nPLATES = [\n' % (tuple(round(float(x), 4) for x in box[0]),
                                                      tuple(round(float(x), 4) for x in box[1]))
    for p in plates:
        body += '    dict(  # normal to %s, %.1f mm thick, %d cut(s)\n' % ('xyz'[p['axis']], p['hi'] - p['lo'], len(p['cuts']))
        body += '        axis=%d, lo=%.4f, hi=%.4f,\n        cuts=[\n' % (p['axis'], p['lo'], p['hi'])
        body += ''.join(fmt_cut(c) for c in p['cuts'])
        body += '        ],\n    ),\n'
    body += ']\n\npart = None\nfor p in PLATES:\n    b = plate(p["axis"], p["lo"], p["hi"], p["cuts"], BOX)\n' \
            '    part = b if part is None else part + b\n' \
            '\n# ---- no re-entrant (concave) edges: fillet them ---------------------------------\n' \
            'R_FILLET = %.2f    # fillet radius on re-entrant edges (mm)\n' \
            'part, n_fillet, n_skip = fillet_reentrant(part.clean(), R_FILLET)\n' \
            'print("re-entrant edges filleted: %%d%%s" %% (n_fillet, ", skipped: %%d" %% n_skip if n_skip else ""))\n' % r_fillet
    tail = tail.replace('orig_exact', name)
    open(path, 'w').write(head + body + '\npart = part.clean()' + tail)
