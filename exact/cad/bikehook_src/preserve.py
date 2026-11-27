"""Bike hook preserve regions (exact): six bosses and the J cradle (sweep of a tapered section) with its rib."""
import math
import numpy as np
import cadquery as cq

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


if __name__ == '__main__':
    c = end_fillets(cradle_with_rib())
    print('cradle', c.val().isValid(), round(c.val().Volume(), 1))
    cq.exporters.export(c, 'cradle.step')


def inflated(d):
    """The cradle grown by about d on every face (no fillets): the region the organic fit leaves out."""
    return cradle(d).union(rib(d))
