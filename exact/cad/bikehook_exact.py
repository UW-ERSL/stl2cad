"""Bike hook (bikehook.stl), written by stl2cad.

Preserve regions (exact): six bosses D15 x 5 with D5.6 holes; the J cradle, one sweep of a
tapered section (8 mm at the centre, 4 mm at the sides) along lines and R15 arcs, with a rib under the bottom.
Connector: a ruled loft through 11 horizontal sections (12 corners each, R2 tangent fillets) and 6 curved
arms (natural cubic splines; elliptical sections at spline points at most 7 mm apart; here ruled lofts between
them, in SolidWorks one loft along the spline). Cut last: above the boss tops, the tyre space inside the J, the holes.

python bikehook_exact.py [out.step]
"""
import math, sys
import numpy as np
import cadquery as cq
from cadquery import Vector


BOSSES = [(-10, 25), (-10, -25), (34, 45), (34, -45), (78, 25), (78, -25)]
R_BOSS, R_HOLE, H_BOSS = 7.5, 2.8, 5.0

# J path in the xz plane: the inner (tyre-side) face of the cradle at its centre
R_J, C_L, C_R, Z_FLAT = 15.0, (-5.0, -63.0), (5.0, -63.0), -78.0
Z_LEG, Z_LIP = -46.0, -61.0           # ends of the two straight legs
# section: u across the cradle (global y), v outward from the path
W, T, E = 25.0, 8.0, 2.0              # half width, centre thickness, drop of each face at the side
SLOPE = 0.2437                        # taper slope (13.7 deg)
R_BLEND, R_SIDE = 30.0, 1.4           # flat-to-taper blend, side corners
END_R = 1.9                           # end fillets (3 at the centre in the mesh, 1.75 at |y| = 22: a constant 1.9 fits the 4 mm edge)
# rib under the bottom: |y| <= RIB_W, flat at RIB_Z, its ends arcs
RIB_W, RIB_Z, RIB_F = 13.0, -90.0, 1.5


def path(d=0.0):
    """The J path; d > 0 moves it d towards the tyre (for an inflated cradle) and extends both ends by d."""
    a = math.radians(45); r = R_J - d
    return (cq.Workplane('XZ')
            .moveTo(C_R[0] + r, Z_LEG + d).lineTo(C_R[0] + r, C_R[1])
            .threePointArc((C_R[0] + r * math.cos(a), C_R[1] - r * math.sin(a)), (C_R[0], C_R[1] - r))
            .lineTo(C_L[0], C_L[1] - r)
            .threePointArc((C_L[0] - r * math.cos(a), C_L[1] - r * math.sin(a)), (C_L[0] - r, C_L[1]))
            .lineTo(C_L[0] - r, Z_LIP + d))



def fillet_polygon(pts, radii):
    """Closed polygon with a fillet of radii[i] at corner i: list of ('L', p) and ('A', mid, end) steps from the
    first tangent point."""
    n = len(pts); P = [np.array(p, float) for p in pts]
    tang = []
    for i in range(n):
        a, b, c = P[i - 1], P[i], P[(i + 1) % n]
        u, w = (a - b) / np.linalg.norm(a - b), (c - b) / np.linalg.norm(c - b)
        phi = math.acos(np.clip(u @ w, -1, 1))            # interior angle
        r = radii[i]
        t = r / math.tan(phi / 2)
        t1, t2 = b + u * t, b + w * t
        bis = (u + w) / np.linalg.norm(u + w)
        ctr = b + bis * r / math.sin(phi / 2)
        m = ctr + (b - ctr) / np.linalg.norm(b - ctr) * r
        tang.append((t1, m, t2))
    return tang


def section_face(pl, d=0.0):
    """Section on plane pl (local x = u along y, local y = v outward), corners filleted; d inflates it."""
    u0 = W - E / SLOPE
    Wd, T0, T1 = W + d, 0.0, T + 2 * d
    pts = [(-u0, T0), (u0, T0), (Wd, E), (Wd, T1 - E), (u0, T1), (-u0, T1), (-Wd, T1 - E), (-Wd, E)]
    radii = [R_BLEND, R_BLEND, R_SIDE, R_SIDE, R_BLEND, R_BLEND, R_SIDE, R_SIDE]
    tg = fillet_polygon(pts, radii)
    w = cq.Workplane(pl).moveTo(*tg[0][2])
    for i in range(1, len(pts)):
        w = w.lineTo(*tg[i][0]).threePointArc(tuple(tg[i][1]), tuple(tg[i][2]))
    w = w.lineTo(*tg[0][0]).threePointArc(tuple(tg[0][1]), tuple(tg[0][2])).close()
    return w


def cradle(d=0.0):
    P = path(d)
    pl = cq.Plane(origin=(C_R[0] + R_J - d, 0.0, Z_LEG + d), xDir=(0, 1, 0), normal=(0, 0, -1))   # local y = +x
    return section_face(pl, d).sweep(P, transition='right')


def rib(d=0.0):
    """Rib under the bottom, |y| <= RIB_W: flat at RIB_Z between x = -4 and 4, then arcs up (centre (-4, -70)
    R 20 on the lip side, until it enters the cradle wall; centre (4, -63) R 27 on the leg side, up to x = 20,
    where the organic trunk takes over). Closed inside the cradle wall."""
    CL, RL = np.array([-4.0, RIB_Z + 20.0]), 20.0 + d
    CR, RR = np.array([4.0, RIB_Z + 27.0]), 27.0 + d
    R_MID = R_J + T / 2
    pt = lambda c, r, a: c + r * np.array([math.cos(a), math.sin(a)])
    # lip side: from -90 deg towards -180 deg and beyond, until the arc is R_MID from the J's lip centre
    aL = next(a for a in np.radians(np.arange(-90, -270, -0.05)) if np.linalg.norm(pt(CL, RL, a) - C_L) <= R_MID)
    # leg side: from -90 deg upwards until x = 20
    aR = -math.pi / 2 + math.asin((20.0 - CR[0]) / RR)
    pl, pr = pt(CL, RL, aL), pt(CR, RR, aR)
    bl, br = (CL[0], RIB_Z - d), (CR[0], RIB_Z - d)
    qr = np.array(C_R) + R_MID * (pr - C_R) / np.linalg.norm(pr - C_R)      # into the wall
    jr, jl = (C_R[0], C_R[1] - R_MID), (C_L[0], C_L[1] - R_MID)
    def arcmid(c, r, p, q):
        a0, a1 = math.atan2(p[1] - c[1], p[0] - c[0]), math.atan2(q[1] - c[1], q[0] - c[0])
        if a1 - a0 > math.pi: a1 -= 2 * math.pi
        if a0 - a1 > math.pi: a1 += 2 * math.pi
        return tuple(pt(c, r, (a0 + a1) / 2))
    w = (cq.Workplane('XZ').moveTo(*pl).threePointArc(arcmid(CL, RL, pl, bl), bl).lineTo(*br)
         .threePointArc(arcmid(CR, RR, br, pr), tuple(pr)).lineTo(*qr)
         .threePointArc(arcmid(np.array(C_R), R_MID, qr, jr), jr).lineTo(*jl)
         .threePointArc(arcmid(np.array(C_L), R_MID, jl, pl), tuple(pl)).close())
    return w.extrude(RIB_W + d, both=True)


def cradle_with_rib():
    body = cradle().union(rib())
    es = [e for f in body.faces().vals() if f.geomType() == 'PLANE' and abs(abs(f.Center().y) - RIB_W) < 1e-6
          and f.Center().z < -80 for e in f.Edges()
          if (e.geomType() == 'LINE' and abs(e.Center().z - RIB_Z) < 1e-3) or (e.geomType() == 'CIRCLE' and e.Center().z < -85.5)]
    # bottom edges only: OCC refuses the concave base fillet (the rib's lip-side arc meets the cradle at a
    # grazing angle)
    body = body.newObject([body.val()]).newObject(es).fillet(RIB_F)
    return body


def end_fillets(body):
    for z, x in ((Z_LEG, 24.0), (Z_LIP, -24.0)):
        f = body.faces(cq.selectors.NearestToPointSelector((x, 0, z))).val()
        es = f.Edges()
        body = body.newObject([body.val()]).newObject(es).fillet(END_R)
    return body


def bosses():
    s = None
    for x, y in BOSSES:
        b = cq.Workplane('XY').center(x, y).circle(R_BOSS).extrude(H_BOSS)
        s = b if s is None else s.union(b)
    return s




# ---- connector ----
R_F = 2.0
SECTIONS = [  # (z, corners (x, y))
    (-28.0, [(16.0, 9.1), (17.0, -17.7), (21.1, -22.7), (25.7, -21.5), (20.0, -3.6), (20.1, 9.0), (38.0, 4.6), (45.8, 9.4), (41.9, 11.9), (23.9, 11.8), (25.0, 25.0), (19.1, 24.8)]),
    (-34.0, [(18.1, 17.3), (17.7, -8.0), (18.5, -20.0), (23.2, -20.9), (24.3, -7.6), (23.2, 2.8), (32.9, 1.1), (44.2, 5.5), (40.9, 9.1), (23.9, 8.1), (24.4, 22.9), (19.5, 23.2)]),
    (-40.0, [(21.8, 22.2), (19.4, 21.1), (19.5, -19.8), (22.9, -20.0), (25.4, -2.1), (33.9, 2.3), (43.9, 0.3), (46.0, 3.1), (43.1, 7.2), (28.2, 6.1), (24.6, 10.5), (24.0, 20.3)]),
    (-46.0, [(22.9, 22.3), (20.2, 20.5), (20.2, -20.3), (24.8, -19.9), (25.1, -7.8), (27.7, 0.7), (37.7, 0.9), (44.5, -2.4), (46.8, 1.2), (43.9, 5.8), (27.1, 4.6), (25.3, 18.6)]),
    (-52.0, [(22.0, 24.9), (20.0, 14.8), (20.0, -14.8), (26.0, -25.0), (28.0, -1.4), (33.8, 2.4), (57.0, -16.0), (58.3, -9.0), (47.2, -8.2), (46.0, 4.7), (29.0, 4.7), (26.0, 25.0)]),
    (-58.0, [(20.0, 17.0), (20.0, -17.4), (25.9, -25.1), (29.4, -7.2), (27.3, -5.6), (31.1, 6.5), (51.9, -16.4), (58.5, -10.2), (47.4, -6.0), (45.2, 4.0), (30.7, 5.0), (26.0, 24.8)]),
    (-64.0, [(20.0, 14.8), (20.0, -14.8), (25.8, -25.0), (29.5, -11.8), (41.6, -9.4), (48.0, -14.9), (55.9, -13.6), (54.3, -7.4), (47.5, -4.4), (43.1, 4.3), (32.8, 5.0), (22.2, 26.5)]),
    (-70.0, [(18.3, 14.8), (18.3, -17.0), (20.7, -25.1), (24.5, -24.7), (26.7, -18.3), (36.5, -13.3), (50.8, -15.6), (51.4, -6.7), (34.8, 7.7), (26.8, 16.3), (24.6, 25.1), (20.7, 24.7)]),
    (-76.0, [(16.0, 25.0), (16.2, -24.7), (20.9, -24.8), (28.7, -16.8), (43.6, -13.7), (45.5, -12.0), (46.6, -9.1), (45.9, -5.8), (42.5, -0.7), (38.8, 6.2), (23.9, 17.2), (21.5, 25.0)]),
    (-81.0, [(16.0, 22.3), (16.0, -23.8), (18.3, -19.4), (19.5, -13.5), (30.0, -12.6), (33.1, -13.4), (37.7, -12.3), (40.7, -8.8), (33.3, 8.2), (26.0, 12.9), (19.8, 13.6), (16.8, 22.7)]),
    (-85.0, [(16.0, 13.0), (16.0, -12.5), (16.9, -13.0), (19.1, -12.4), (19.7, -10.8), (20.9, -10.3), (21.3, -9.4), (21.5, 7.0), (21.2, 9.0), (19.7, 10.9), (19.3, 12.2), (17.4, 13.0)]),
]

ARMS = [  # spline points (natural cubic, chord-length parameters); sections (point index, half-axes a, b,
          # major-axis direction)
    dict(points=[(76.2, -20.1, 2.2), (72.8, -18.74, -3.01), (69.74, -17.55, -8.46), (67.43, -16.73, -14.33), (65.91, -16.29, -20.73), (64.84, -16.05, -27.23), (63.73, -15.78, -33.72), (62.11, -15.25, -40.09), (59.8, -14.36, -46.2), (57.01, -13.24, -52.06), (54.0, -12.0, -57.8)],
         sec=[(0, 6.35, 3.34, (-0.6618, -0.7058, 0.2527)), (1, np.float64(5.96), np.float64(3.67), (-0.6783, -0.6903, 0.252)), (2, np.float64(5.56), np.float64(4.0), (-0.7132, -0.6641, 0.2241)), (3, 5.16, 4.33, (-0.7554, -0.6333, 0.1682)), (4, np.float64(5.12), np.float64(4.36), (-0.7671, -0.6313, 0.1143)), (5, np.float64(5.07), np.float64(4.39), (-0.7692, -0.6314, 0.0987)), (6, np.float64(5.02), np.float64(4.43), (-0.7669, -0.6308, 0.1182)), (7, np.float64(4.98), np.float64(4.46), (-0.7536, -0.6361, 0.1658)), (8, np.float64(4.93), np.float64(4.49), (-0.7322, -0.6488, 0.2071)), (9, np.float64(4.89), np.float64(4.52), (-0.7165, -0.6597, 0.2267)), (10, 4.84, 4.55, (-0.7097, -0.6652, 0.2318))]),
    dict(points=[(-5.6, 22.1, 2.2), (-1.91, 19.08, -0.29), (1.72, 16.35, -3.15), (5.18, 14.22, -6.65), (8.4, 12.75, -10.69), (11.85, 11.72, -15.76), (15.05, 11.28, -21.08), (18.07, 11.28, -26.51), (21.0, 11.5, -32.0)],
         sec=[(0, 3.73, 3.32, (-0.5625, -0.0234, -0.8265)), (1, np.float64(3.7), np.float64(3.34), (-0.5885, -0.018, -0.8083)), (2, np.float64(3.66), np.float64(3.36), (-0.6534, 0.0267, -0.7565)), (3, np.float64(3.63), np.float64(3.38), (-0.7301, 0.0785, -0.6788)), (4, 3.6, 3.4, (-0.7858, 0.114, -0.6078)), (5, np.float64(3.91), np.float64(3.37), (-0.8457, -0.1198, -0.52)), (6, np.float64(4.22), np.float64(3.34), (-0.802, -0.3958, -0.4473)), (7, np.float64(4.52), np.float64(3.31), (-0.6726, -0.6344, -0.3809)), (8, 4.83, 3.28, (-0.5162, -0.798, -0.3109))]),
    dict(points=[(32.1, 40.0, 2.2), (30.02, 36.09, -1.33), (27.95, 32.32, -5.03), (25.9, 28.85, -9.01), (24.53, 26.6, -11.98), (23.3, 24.41, -15.05), (22.01, 21.4, -19.32), (21.29, 18.21, -23.59), (21.03, 14.89, -27.81), (21.0, 11.5, -32.0)],
         sec=[(0, 6.52, 4.89, (-0.8728, 0.4871, -0.0299)), (1, np.float64(6.32), np.float64(5.11), (-0.8556, 0.5141, -0.0596)), (2, np.float64(6.11), np.float64(5.32), (-0.8287, 0.5543, -0.0769)), (3, 5.91, 5.54, (-0.7972, 0.5972, -0.088)), (4, np.float64(6.25), np.float64(5.3), (-0.3386, 0.8242, -0.454)), (5, 6.58, 5.07, (0.2644, 0.7417, -0.6165)), (6, np.float64(6.93), np.float64(4.99), (0.3061, 0.7364, -0.6033)), (7, np.float64(7.28), np.float64(4.92), (0.3445, 0.7258, -0.5954)), (8, np.float64(7.64), np.float64(4.84), (0.3824, 0.7161, -0.5839)), (9, 7.99, 4.77, (0.4199, 0.7048, -0.5718))]),
    dict(points=[(31.5, -40.4, 2.2), (29.62, -36.41, -1.19), (27.78, -32.59, -4.79), (26.04, -29.14, -8.79), (24.44, -25.94, -13.33), (23.14, -22.88, -18.05), (22.19, -19.5, -23.03), (21.74, -15.73, -27.8), (21.68, -11.68, -32.36), (21.8, -7.5, -36.8)],
         sec=[(0, 5.06, 4.19, (-0.8777, -0.4737, -0.0727)), (1, np.float64(4.92), np.float64(4.3), (-0.8764, -0.4766, -0.0689)), (2, np.float64(4.78), np.float64(4.42), (-0.8706, -0.4885, -0.0576)), (3, 4.63, 4.53, (-0.8648, -0.5001, -0.0454)), (4, np.float64(4.93), np.float64(4.28), (-0.9226, 0.1173, 0.3675)), (5, 5.24, 4.02, (0.474, -0.6832, -0.5555)), (6, np.float64(5.44), np.float64(3.95), (0.4805, -0.6752, -0.5597)), (7, np.float64(5.64), np.float64(3.89), (0.4876, -0.6548, -0.5774)), (8, np.float64(5.85), np.float64(3.82), (0.4968, -0.6408, -0.5853)), (9, 6.05, 3.75, (0.5048, -0.6344, -0.5854))]),
    dict(points=[(-5.7, -22.0, 2.2), (-0.81, -18.7, -1.27), (3.87, -15.96, -5.43), (8.1, -14.03, -10.45), (11.83, -12.58, -16.0), (14.69, -11.44, -21.05), (17.24, -10.21, -26.23), (19.58, -8.88, -31.5), (21.8, -7.5, -36.8)],
         sec=[(0, 3.54, 3.51, (0.6581, -0.7121, 0.2445)), (1, np.float64(3.62), np.float64(3.51), (0.6512, -0.7177, 0.2467)), (2, np.float64(3.71), np.float64(3.51), (0.64, -0.7242, 0.2568)), (3, np.float64(3.8), np.float64(3.5), (0.6414, -0.7227, 0.2576)), (4, 3.89, 3.5, (0.6567, -0.7167, 0.2345)), (5, np.float64(4.05), np.float64(3.47), (0.6085, -0.7808, 0.142)), (6, np.float64(4.2), np.float64(3.44), (0.5459, -0.8365, 0.0475)), (7, np.float64(4.36), np.float64(3.41), (0.4655, -0.8846, -0.029)), (8, 4.51, 3.38, (0.3745, -0.9232, -0.0858))]),
    dict(points=[(74.0, 21.6, 2.2), (70.39, 18.34, -2.65), (66.54, 15.31, -7.46), (62.19, 12.83, -12.16), (57.18, 11.26, -16.56), (51.59, 10.67, -20.5), (45.85, 10.19, -24.24), (40.3, 8.93, -28.08), (35.06, 6.9, -32.03), (30.0, 4.5, -36.0)],
         sec=[(0, 3.52, 3.46, (0.3721, -0.8715, 0.3193)), (1, np.float64(3.52), np.float64(3.46), (0.3585, -0.8816, 0.307)), (2, np.float64(3.52), np.float64(3.46), (0.321, -0.9097, 0.2635)), (3, np.float64(3.53), np.float64(3.45), (0.2435, -0.9521, 0.1849)), (4, np.float64(3.53), np.float64(3.45), (0.116, -0.9892, 0.09)), (5, np.float64(3.53), np.float64(3.45), (0.0315, -0.9985, 0.0447)), (6, np.float64(3.54), np.float64(3.44), (0.09, -0.9932, 0.0735)), (7, np.float64(3.54), np.float64(3.44), (0.2034, -0.9689, 0.1406)), (8, np.float64(3.54), np.float64(3.44), (0.2738, -0.9429, 0.1896)), (9, 3.55, 3.43, (0.2958, -0.9335, 0.203))]),
]

def natural_spline(P):
    """The natural cubic spline through points P with chord-length parameters (SolidWorks' CreateSpline2 with natural
    ends: a loft along it matched SolidWorks to 0.03 % in sw_probe_bars). Returns (spline, knots)."""
    from scipy.interpolate import CubicSpline
    P = np.asarray(P, float)
    d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    if len(P) == 2:
        return (lambda t, nu=0: P[0] + (P[1] - P[0]) * (t / d[1]) if nu == 0 else np.broadcast_to((P[1] - P[0]) / d[1], (np.size(t), 3)) if np.ndim(t) else (P[1] - P[0]) / d[1]), d
    cs = CubicSpline(d, P, bc_type='natural')
    return cs, d


def arc_lengths(f, d, n=2000):
    t = np.linspace(0.0, d[-1], n)
    X = np.array([f(x) for x in t]) if not hasattr(f, 'derivative') else f(t)
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(X, axis=0), axis=1))]
    return t, s


def spline(P):
    return natural_spline(np.asarray(P, float))


def sample(P, n=400):
    f, d = spline(P)
    t = np.linspace(d[0], d[-1], n)
    X = np.array([np.asarray(f(x), float).reshape(-1)[:3] for x in t])
    T = np.array([np.asarray(f(x, 1), float).reshape(-1)[:3] for x in t])
    T /= np.linalg.norm(T, axis=1)[:, None]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(X, axis=0), axis=1))]
    return X, T, s


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


def radii(V, r):
    """Fillet radius per corner, reduced so each fillet's tangent length is at most 0.45 x its sides."""
    n = len(V); out = []
    for i in range(n):
        a, b, c = V[i - 1], V[i], V[(i + 1) % n]
        u, w = (a - b) / np.linalg.norm(a - b), (c - b) / np.linalg.norm(c - b)
        phi = math.acos(np.clip(u @ w, -1, 1))
        tmax = 0.45 * min(np.linalg.norm(a - b), np.linalg.norm(c - b))
        out.append(min(r, tmax * math.tan(phi / 2)))
    return out


def section_wire(V, z, r=R_F):
    tg = fillet_polygon(V, radii(V, r))
    P3 = lambda p: cq.Vector(float(p[0]), float(p[1]), float(z))
    es = []
    for i in range(len(V)):
        t1, m, t2 = tg[i]
        es.append(cq.Edge.makeThreePointArc(P3(t1), P3(m), P3(t2)))
        n1 = tg[(i + 1) % len(V)][0]
        if np.linalg.norm(n1 - t2) > 1e-6:
            es.append(cq.Edge.makeLine(P3(t2), P3(n1)))
    return cq.Wire.assembleEdges(es)



def trunk():
    return cq.Solid.makeLoft([section_wire(np.array(V), z) for z, V in SECTIONS], ruled=True)


def arm(a):
    """The arm: ruled lofts between the ellipses at its spline points (SolidWorks lofts the same sections with the
    spline as centerline, smooth)."""
    sec = [(i, a_, b_, np.array(m)) for i, a_, b_, m in a['sec']]
    return arm_loft(np.array(a['points'], float), sec)


def keepout():
    """The tyre's space inside the J."""
    a = math.radians(45)
    w = (cq.Workplane('XZ').moveTo(C_R[0] + R_J, Z_LEG + 20).lineTo(C_R[0] + R_J, C_R[1])
         .threePointArc((C_R[0] + R_J * math.cos(a), C_R[1] - R_J * math.sin(a)), (C_R[0], Z_FLAT))
         .lineTo(C_L[0], Z_FLAT)
         .threePointArc((C_L[0] - R_J * math.cos(a), C_L[1] - R_J * math.sin(a)), (C_L[0] - R_J, C_L[1]))
         .lineTo(C_L[0] - R_J, Z_LEG + 20).close())
    return w.extrude(40, both=True).val()


def build():
    part = end_fillets(cradle_with_rib()).val()
    for x, y in BOSSES:
        part = part.fuse(cq.Solid.makeCylinder(R_BOSS, H_BOSS, cq.Vector(x, y, 0)))
    core = trunk()                      # arms fused to the trunk first, then the cradle and bosses
    for a in ARMS:
        core = core.fuse(arm(a))
    part = part.fuse(core)
    part = part.cut(cq.Solid.makeBox(300, 300, 50, cq.Vector(-150, -150, H_BOSS)))
    part = part.cut(keepout())
    for x, y in BOSSES:
        part = part.cut(cq.Solid.makeCylinder(R_HOLE, 60, cq.Vector(x, y, -30)))
    return part.clean()


if __name__ == '__main__':
    part = build()
    print('valid %s, solids %d, faces %d, volume %.1f mm^3' % (part.isValid(), len(part.Solids()), len(part.Faces()), part.Volume()))
    cq.exporters.export(cq.Workplane().add(part), sys.argv[1] if len(sys.argv) > 1 else 'bikehook_exact.step')
