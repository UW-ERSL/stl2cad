"""stent_exact: 20 wires (10 right-handed, 10 left-handed), wire diameter 0.12, helix radius 1.5, pitch 16.
One solid per wire: OCC does not fuse tubes of equal radius whose centre lines meet; SolidWorks
does (Combine or merged sweeps)."""
import math
import cadquery as cq

AXIS_POINT = (1.566, 1.566, 0.0)
AXIS = (0.0, 0.0, 1.0)
E1 = (1.0, 0.0, 0.0)      # angle 0 direction
R = 1.5             # helix radius
WIRE_R = 0.06        # wire radius
CUTS = (None, 16.03)   # flat cuts normal to the axis (axial coordinate)
WIRES = [  # (pitch, right-handed, angle at the start z0 (deg, from E1), z0, z1)
    (16.0, True, 0.0, 0.03, 16.15),
    (16.0, True, 36.0, 0.03, 16.15),
    (16.0, True, 72.0, 0.03, 16.15),
    (16.0, True, 108.0, 0.03, 16.15),
    (16.0, True, 144.0, 0.03, 16.15),
    (16.0, True, 180.0, 0.03, 16.15),
    (16.0, True, 216.0, 0.03, 16.15),
    (16.0, True, 252.0, 0.03, 16.15),
    (16.0, True, 288.0, 0.03, 16.15),
    (16.0, True, 324.0, 0.03, 16.15),
    (16.0, False, 0.0, 0.03, 16.15),
    (16.0, False, 36.0, 0.03, 16.15),
    (16.0, False, 72.0, 0.03, 16.15),
    (16.0, False, 108.0, 0.03, 16.15),
    (16.0, False, 144.0, 0.03, 16.15),
    (16.0, False, 180.0, 0.03, 16.15),
    (16.0, False, 216.0, 0.03, 16.15),
    (16.0, False, 252.0, 0.03, 16.15),
    (16.0, False, 288.0, 0.03, 16.15),
    (16.0, False, 324.0, 0.03, 16.15),
]


def frame():
    a = cq.Vector(*AXIS).normalized()
    e1 = cq.Vector(*E1).normalized()
    return a, e1, a.cross(e1)


def wire_solid(pitch, right, start_deg, z0, z1):
    a, e1, e2 = frame()
    phi0 = math.radians(start_deg)                   # angle at the start of the wire
    h = cq.Wire.makeHelix(pitch, z1 - z0, R, lefthand=not right)   # local: starts at (R, 0, 0), up z
    # local frame -> part: x along the start angle, z along the axis, origin at the axis point + z0
    xd = e1 * math.cos(phi0) + e2 * math.sin(phi0)
    loc = cq.Location(cq.Plane(origin=cq.Vector(*AXIS_POINT) + a * z0, xDir=xd, normal=a))
    h = h.moved(loc)
    prof = cq.Wire.makeCircle(WIRE_R, h.startPoint(), h.tangentAt(0))
    s = cq.Solid.sweep(prof, [], h, True, True)
    lo, hi = CUTS
    if lo is not None or hi is not None:
        big = 10 * (R + pitch + abs(z1 - z0))
        zl = lo if lo is not None else z0 - big
        zh = hi if hi is not None else z1 + big
        box = cq.Solid.makeBox(2 * big, 2 * big, zh - zl, pnt=cq.Vector(-big, -big, zl),
                               dir=cq.Vector(0, 0, 1))
        box = box.moved(cq.Location(cq.Plane(origin=cq.Vector(*AXIS_POINT), xDir=e1, normal=a)))
        s = s.intersect(box)
    return s


def build():
    return cq.Compound.makeCompound([wire_solid(*w) for w in WIRES])


result = build()

if __name__ == "__main__":
    import sys
    out = sys.argv[1] if len(sys.argv) > 1 else 'stent_exact.step'
    cq.exporters.export(result, out)
    print('solids', len(result.Solids()), 'sum of volumes', round(sum(s.Volume() for s in result.Solids()), 6))
