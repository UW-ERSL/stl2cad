"""Arms of the bike hook as smooth curved members with elliptical sections.

Each arm runs from its boss (the frame's anchor inside the boss) to a point inside the trunk (the frame joint where
the arm enters it). Its path is a natural cubic spline (chord-length parameters, as SolidWorks' CreateSpline2 with
natural ends) through the two end points and N_FREE interior points, which are fitted by least squares to the
centroids of the mesh sections normal to the path; the sections are measured again on the fitted path (ITER
passes). Each section piece (the one round the path) gives an ellipse of the same area and principal directions
(semi-axes 2 sqrt(lambda_i), scaled to the area). The loft sections are ellipses at the two ends and at the interior
stations where the linear change along the arc length departs from the measured axes by more than SEC_TOL.
Stations are used only outside the bosses and above the trunk (TRUNK_TOP), and only before the section grows where
two arms merge (MERGE); past the last station the path keeps the last measured section into the trunk.
"""
import math, sys, os
import numpy as np, trimesh, cadquery as cq
from cadquery import Vector
from scipy.optimize import least_squares
from shapely.geometry import Point
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'to'))
import frame3d

N_FREE, ITER, N_ST = 2, 3, 40
SPACING = 7.0      # mm: spline points (and sections) at most this far apart along the arm, so that a loft
                   # through the sections follows the spline in SolidWorks and in CadQuery alike
TRUNK_TOP, BOSS_CLEAR, SEC_TOL, MERGE, ROUND = -26.0, 2.5, 0.08, 1.6, 1.15


def arm_chains(R, bosses, r_boss):
    """(start, end) of each arm: the bars from a hub; a bar ending at a free junction above the trunk is continued
    by the bar from that junction into the trunk."""
    J = np.array(R['joints']); B = R['bars']
    out = []
    for b in B:
        ks = b['joints']
        for kk in (ks, ks[::-1]):
            if R['kind'][kk[0]] == 'hub' and J[kk[0]][2] > -1 and J[kk[-1]][2] < 0:
                pts = [J[k] for k in kk]
                end = J[kk[-1]]
                if end[2] > -28:                               # a junction above the trunk: continue downward
                    for c in B:
                        if c is b:
                            continue
                        cs = c['joints']
                        for cc in (cs, cs[::-1]):
                            if np.allclose(J[cc[0]], end) and J[cc[-1]][2] < end[2] - 5:
                                pts += [J[k] for k in cc[1:]]
                out.append(np.array(pts))
                break
    return out


def spline(P):
    return frame3d.natural_spline(np.asarray(P, float))


def sample(P, n=400):
    f, d = spline(P)
    t = np.linspace(d[0], d[-1], n)
    X = np.array([np.asarray(f(x), float).reshape(-1)[:3] for x in t])
    T = np.array([np.asarray(f(x, 1), float).reshape(-1)[:3] for x in t])
    T /= np.linalg.norm(T, axis=1)[:, None]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(X, axis=0), axis=1))]
    return X, T, s


def frame_vecs(t):
    a = np.eye(3)[int(np.argmin(np.abs(t)))]
    u = np.cross(t, a); u /= np.linalg.norm(u)
    return u, np.cross(t, u)


def section_at(mesh, p, t):
    """Ellipse of the section piece round p in the plane normal to t: (centroid 3D, a, b, major unit vector,
    area), or None."""
    s = mesh.section(plane_origin=p, plane_normal=t)
    if s is None:
        return None
    u, v = frame_vecs(t)
    V = s.vertices - p
    s2 = trimesh.path.Path2D(entities=[e.copy() for e in s.entities], vertices=np.c_[V @ u, V @ v])
    polys = s2.polygons_full
    if not polys:
        return None
    o = Point(0, 0)
    pg = min(polys, key=lambda g: g.distance(o))
    if pg.distance(o) > 3.0:
        return None
    A = pg.area
    c = np.array(pg.centroid.coords[0])
    # second moments of the region (triangulation-free: sample points inside)
    xs = np.array(pg.exterior.coords)
    x0, y0, x1, y1 = pg.bounds
    rng = np.random.default_rng(0)
    q = np.c_[rng.uniform(x0, x1, 4000), rng.uniform(y0, y1, 4000)]
    from shapely import contains_xy
    q = q[contains_xy(pg, q[:, 0], q[:, 1])]
    if len(q) < 50:
        return None
    C = np.cov((q - c).T)
    w, E = np.linalg.eigh(C)
    a, b = 2 * math.sqrt(max(w[1], 1e-9)), 2 * math.sqrt(max(w[0], 1e-9))
    k = math.sqrt(A / (math.pi * a * b)); a, b = a * k, b * k
    major = E[0, 1] * u + E[1, 1] * v
    return p + c[0] * u + c[1] * v, a, b, major / np.linalg.norm(major), A


def stations(mesh, P, bosses, r_boss):
    X, T, s = sample(P)
    out = []
    for si in np.linspace(0, s[-1], N_ST + 2)[1:-1]:
        i = int(np.searchsorted(s, si))
        p, t = X[i], T[i]
        if p[2] < TRUNK_TOP:
            continue
        if any(math.hypot(p[0] - x, p[1] - y) < r_boss + BOSS_CLEAR and p[2] > -BOSS_CLEAR for x, y in bosses):
            continue
        sec = section_at(mesh, p, t)
        if sec is not None:
            out.append((si, p, t) + sec)
    # where two arms merge (or an arm meets the trunk) the section piece grows: stations from the first one whose
    # area exceeds MERGE x the median of the first half are dropped
    if len(out) >= 4:
        ref = float(np.median([x[7] for x in out[:max(2, len(out) // 2)]]))
        for k, x in enumerate(out):
            if x[7] > MERGE * ref or max(x[4], x[5]) > MERGE * math.sqrt(ref / math.pi) * 1.5:
                out = out[:k]
                break
    return out


def fit_path(P0, cents):
    """Interior points of the natural spline from P0[0] to P0[-1] fitted to the centroids."""
    a, b = P0[0], P0[-1]
    X, T, s = sample(P0)
    init = [X[int(np.searchsorted(s, s[-1] * (k + 1) / (N_FREE + 1)))] for k in range(N_FREE)]

    def res(z):
        Q = np.vstack([a, z.reshape(-1, 3), b])
        Y, _, _ = sample(Q, 300)
        d = np.linalg.norm(cents[:, None, :] - Y[None], axis=2).min(axis=1)
        return d
    r = least_squares(res, np.ravel(init), x_scale=5.0, max_nfev=200)
    return np.vstack([a, r.x.reshape(-1, 3), b]), float(np.sqrt(np.mean(r.fun ** 2)))


def fit_arm(mesh, P, bosses, r_boss):
    Q = np.vstack([P[0], P[-1]]) if len(P) < 2 else P
    Q = P
    rms = None
    for _ in range(ITER):
        st = stations(mesh, Q, bosses, r_boss)
        if len(st) < 3:
            break
        cents = np.array([x[3] for x in st])
        Q, rms = fit_path(Q, cents)
    st = stations(mesh, Q, bosses, r_boss)
    return Q, st, rms


def loft_sections(Q, st):
    """Sections (arc length, a, b, major vector) at the ends and where the axes depart from a linear change."""
    X, T, s = sample(Q)
    S = np.array([x[0] for x in st]); A = np.array([x[4] for x in st]); B = np.array([x[5] for x in st])
    M = [x[6] for x in st]
    Tn = [x[2] for x in st]
    for k in range(1, len(M)):
        # a nearly round section has no clear major direction: it takes the previous one (projected on its plane),
        # so the loft does not twist; otherwise the sign closer to the previous one
        if A[k] < ROUND * B[k]:
            m = M[k - 1] - (M[k - 1] @ Tn[k]) * Tn[k]
            M[k] = m / np.linalg.norm(m)
        elif M[k] @ M[k - 1] < 0:
            M[k] = -M[k]
    keep = {0, len(S) - 1}
    stack = [(0, len(S) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        lin = lambda arr, k: arr[i] + (arr[j] - arr[i]) * (S[k] - S[i]) / (S[j] - S[i])
        dev = [(max(abs(A[k] - lin(A, k)) / A[k], abs(B[k] - lin(B, k)) / B[k]), k) for k in range(i + 1, j)]
        e, k = max(dev)
        if e > SEC_TOL:
            keep.add(k); stack += [(i, k), (k, j)]
    ks = sorted(keep)
    secs = [(float(S[k]), float(A[k]), float(B[k]), M[k]) for k in ks]
    # the ends of the path take the nearest measured section
    secs[0] = (0.0,) + secs[0][1:]
    secs[-1] = (float(s[-1]),) + secs[-1][1:]
    return secs


def ellipse_wire(p, t, a, b, major):
    major = major - (major @ t) * t; major /= np.linalg.norm(major)
    b = min(a, b)                       # never swap the axes: the ellipse's start point stays on +major, so a loft
    e = cq.Edge.makeEllipse(a, b, Vector(*p), Vector(*t), Vector(*major))   # through the sections does not twist
    return cq.Wire.assembleEdges([e])


def arm_solid(Q, secs):
    """Pieces swept along the spline between consecutive sections (SolidWorks' loft with the spline as centerline:
    the section changes along the arc length)."""
    f, d = spline(Q)
    X, T, s = sample(Q, 600)
    pieces = []
    for (s0, a0, b0, m0), (s1, a1, b1, m1) in zip(secs[:-1], secs[1:]):
        i0, i1 = int(np.searchsorted(s, s0 - 1e-9)), int(np.searchsorted(s, s1 - 1e-9))
        i1 = min(i1, len(s) - 1)
        pts = [Vector(*X[i]) for i in range(i0, i1 + 1, max(1, (i1 - i0) // 40))]
        if pts[-1].sub(Vector(*X[i1])).Length > 1e-9:
            pts.append(Vector(*X[i1]))
        path = cq.Wire.assembleEdges([cq.Edge.makeSpline(pts)])
        w0 = ellipse_wire(X[i0], T[i0], a0, b0, m0)
        w1 = ellipse_wire(X[i1], T[i1], a1, b1, m1)
        pieces.append(cq.Solid.sweep_multi([w0, w1], path, True, False))
    return pieces


def fit_all(mesh, R, bosses, r_boss, log=print):
    arms = []
    for P in arm_chains(R, bosses, r_boss):
        Q, st, rms = fit_arm(mesh, P, bosses, r_boss)
        secs = loft_sections(Q, st)
        log('arm from %s to %s: %d stations, centroid rms %.2f mm, %d sections, axes %s' % (
            np.round(Q[0], 1), np.round(Q[-1], 1), len(st), rms or -1, len(secs),
            ', '.join('%.1fx%.1f' % (2 * a, 2 * b) for _, a, b, _ in secs)))
        arms.append(dict(points=Q, secs=secs))
    return arms


def path_points(Q, secs, min_interior=2, spacing=None):
    """Interpolation points of the arm's spline as written (SolidWorks CreateSpline2): the ends, a point at each
    interior section, and more points at even arc length up to min_interior interior points; the fitted spline
    sampled there. Returns (points, [(point index, a, b, major)])."""
    X, T, s = sample(Q, 800)
    S = [sec[0] for sec in secs]
    extra = []
    n_int = len(S) - 2
    if n_int < min_interior and not spacing:
        for k in range(1, min_interior - n_int + 1):
            extra.append(s[-1] * k / (min_interior - n_int + 1))
    if spacing:                                   # points at most `spacing` apart along the arc (dense sections)
        for i in range(len(S) - 1):
            k = int(math.ceil((S[i + 1] - S[i]) / spacing))
            extra += [S[i] + (S[i + 1] - S[i]) * j / k for j in range(1, k)]
    allS = sorted(set([round(v, 6) for v in S + extra]))
    pts = np.array([[np.interp(v, s, X[:, k]) for k in range(3)] for v in allS])
    pts[0], pts[-1] = Q[0], Q[-1]
    out = []
    for sec in secs:
        i = allS.index(round(sec[0], 6))
        out.append((i, sec[1], sec[2], np.asarray(sec[3], float)))
    return pts, out


def every_point(P, sec):
    """A section at every spline point: the half-axes interpolated linearly in arc length between the fitted
    sections, the major direction interpolated and projected on the point's normal plane (so a loft through the
    sections follows the spline)."""
    f, d = spline(P)
    t, sl = frame3d.arc_lengths(f, d)
    sk = np.interp(d, t, sl)
    have = {i: (a, b, np.asarray(m, float)) for i, a, b, m in sec}
    idx = sorted(have)
    out = []
    for i in range(len(P)):
        tan = np.asarray(f(d[i], 1), float).reshape(-1)[:3]; tan /= np.linalg.norm(tan)
        if i in have:
            a, b, m = have[i]
        else:
            j0 = max(k for k in idx if k < i); j1 = min(k for k in idx if k > i)
            w = (sk[i] - sk[j0]) / (sk[j1] - sk[j0])
            a = (1 - w) * have[j0][0] + w * have[j1][0]; b = (1 - w) * have[j0][1] + w * have[j1][1]
            m1 = have[j1][2] if have[j1][2] @ have[j0][2] >= 0 else -have[j1][2]
            m = (1 - w) * have[j0][2] + w * m1
        m = m - (m @ tan) * tan; m /= np.linalg.norm(m)
        out.append((i, a, min(a, b), m))
    return P, out


def arm_loft(P, sec):
    """The arm in CadQuery: ruled lofts between consecutive ellipses at its spline points (at most SPACING apart),
    fused. OCC's smooth ThruSections through all of them was 4-6 % off on two arms (with its compatibility
    matching or without), and sweeps between them failed on short pieces; SolidWorks lofts the same sections with
    the spline as centerline (smooth), within about 1 % of this."""
    f, d = spline(P)
    ws = []
    for i, a, b, m in sec:
        tan = np.asarray(f(d[i], 1), float).reshape(-1)[:3]; tan /= np.linalg.norm(tan)
        ws.append(ellipse_wire(P[i], tan, a, b, np.asarray(m, float)))
    s = cq.Solid.makeLoft(ws[:2], ruled=True)
    for k in range(1, len(ws) - 1):
        s = s.fuse(cq.Solid.makeLoft(ws[k:k + 2], ruled=True), glue=True)
    return s.clean()


def arm_solid_pts(P, sec):
    """The arm on the spline through P (as SolidWorks has it), sections at point indices."""
    f, d = spline(P)
    t, sl = frame3d.arc_lengths(f, d)
    sk = np.interp(d, t, sl)
    secs = [(float(sk[i]), a, b, m) for i, a, b, m in sec]
    return arm_solid(P, secs)


def fit_from_frame(mesh, R, bosses, r_boss, detail=5, z_trunk=-28.0, z_keep=-30.0, log=print):
    """Arms from the frame: members rising above z_keep (none with both ends on hubs), joints below z_trunk end
    a bar; the bars from each hub fitted (fit_all) and written as spline points and sections."""
    J = np.array(R['joints'])
    keep = [m for m in R['members'] if max(J[m['i']][2], J[m['j']][2]) > z_keep
            and not (R['kind'][m['i']] == 'hub' and R['kind'][m['j']] == 'hub')]
    R2 = dict(R); R2['members'] = keep
    R2['kind'] = [('hub' if (k == 'hub' or J[i][2] < z_trunk) else k) for i, k in enumerate(R['kind'])]
    R2['bars'] = frame3d.make_bars(R2, frame3d.BAR_TOL[detail], frame3d.DETAIL[detail])
    out = []
    for a in fit_all(mesh, R2, bosses, r_boss, log=log):
        P, sec = every_point(*path_points(a['points'], a['secs'], spacing=SPACING))
        out.append(dict(points=[tuple(round(float(v), 2) for v in p) for p in P],
                        sec=[(i, round(a_, 2), round(b_, 2), tuple(round(float(v), 4) for v in m)) for i, a_, b_, m in sec]))
    return out
