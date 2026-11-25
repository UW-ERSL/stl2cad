"""Test stand-in for the SolidWorks COM API, built on CadQuery.

It implements only the calls that stl2cad's SolidWorks scripts make, with the documented
SolidWorks conventions: lengths in metres, Front/Top/Right planes with normals +z/+y/+x, sketch
coordinates from ModelToSketchTransform, CreateArc direction +1 counter-clockwise. Default
extrusion directions are a parameter (CUT_REVERSED), so a test can check that a script corrects
a wrong guess. It checks the script's logic and geometry, not SolidWorks itself.
"""
import math
import os

import numpy as np
import cadquery as cq

CUT_REVERSED = os.environ.get("SWMOCK_CUT_REVERSED", "1") == "1"   # cut goes against the normal
BOSS_REVERSED = os.environ.get("SWMOCK_BOSS_REVERSED", "0") == "1"
FRAMES = {"Front Plane": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
          "Top Plane": ((1, 0, 0), (0, 0, -1), (0, 1, 0)),
          "Right Plane": ((0, 0, -1), (0, 1, 0), (1, 0, 0))}
LAST = {}
K = 1000.0                          # API metres -> mm, so the body is built and saved in mm


class CDispatch:
    """Stands for a COM object returned by late binding."""


class VARIANT:
    def __init__(self, vt, value):
        self.vt, self.value = vt, value


class Point:
    def __init__(self, xyz):
        self.ArrayData = tuple(float(v) for v in xyz)

    def MultiplyTransform(self, T):
        return Point(T @ np.array(self.ArrayData))


class MathUtility:
    def CreatePoint(self, v):
        return Point(v.value if isinstance(v, VARIANT) else v)


class Transform:
    """Model -> sketch. ArrayData as in SolidWorks: rotation R (row-major, applied as p R),
    translation, scale; a point on the sketch plane maps to (a, b, 0)."""
    def __init__(self, frame, offset=0.0, origin=None):
        u, v, n = (np.array(a, float) for a in frame)
        self.M = np.vstack([u, v, n])
        R = np.column_stack([u, v, n])
        t = (0.0, 0.0, -offset) if origin is None else tuple(-(np.asarray(origin, float) @ R))
        self.ArrayData = tuple(R.ravel()) + tuple(float(x) for x in t) + (1.0, 0.0, 0.0, 0.0)

    def __matmul__(self, p):
        return self.M @ p


def earc_edge(e, P, n):
    """CadQuery edge of an elliptical arc ('earc', centre, major point, minor point, start, end, direction) in
    sketch coordinates; P maps a sketch point to a model Vector (mm), n is the sketch normal."""
    c, pa, pb, p1, p2 = (np.array(q, float) for q in e[1:6])
    a, b = np.linalg.norm(pa - c), np.linalg.norm(pb - c)
    X = (pa - c) / a
    Y = np.array([-X[1], X[0]])                     # 90 deg counter-clockwise in the sketch
    t = lambda p: math.degrees(math.atan2((p - c) @ Y / b, (p - c) @ X / a))
    t1, t2 = t(p1), t(p2)
    if e[6] < 0:
        t1, t2 = t2, t1
    t2 = t1 + (t2 - t1) % 360.0
    C = P(c)
    return cq.Edge.makeEllipse(K * a, K * b, C, cq.Vector(*n), P(pa) - C, t1, t2)


class Segment:
    def __init__(self, sm, edge):
        self.sm, self.edge = sm, edge
        self.ConstructionGeometry = False
        self.deleted = False

    def Select4(self, append, data):
        m = self.sm.model
        sk = getattr(self, "sk", None)
        if sk is not None and self.sm._sk is not sk:        # sketch closed: selected as an external segment
            if not append:
                m.ext_sel = []
            pts = [sk.to_model(q) for q in self.edge[1:]]
            V = [cq.Vector(*(K * q)) for q in pts]
            edge = cq.Edge.makeLine(*V) if self.edge[0] == "line" else cq.Edge.makeThreePointArc(*V)
            m.ext_sel.append(("segment", edge, getattr(data, "Mark", 0)))
            m.ext_segs = getattr(m, "ext_segs", []) + [(self, getattr(data, "Mark", 0))]
            return True
        m.selected = self
        if not append:
            m.sel_list = []
        m.sel_list.append((self.edge[0], self))
        return True


class SkPoint:
    """A sketch point made by CreatePoint; once its sketch is closed, Select4 selects it as an external point."""
    def __init__(self, sm, xy):
        self.sm, self.sk, self.xy = sm, sm._sk, (float(xy[0]), float(xy[1]))
        self.X, self.Y, self.Z = self.xy[0], self.xy[1], 0.0

    def Select4(self, append, data):
        m = self.sm.model
        if not append:
            m.ext_sel = []
        m.ext_sel.append(("point", self.sk.to_model(self.xy), getattr(data, "Mark", 0)))
        return True

    def _end(self, i):
        if not hasattr(self, "_ends"):
            self._ends = [SketchPoint(self.sm.model, self.edge[1]), SketchPoint(self.sm.model, self.edge[-1])]
        return self._ends[i]

    def GetStartPoint2(self):
        return self._end(0)

    def GetEndPoint2(self):
        return self._end(1)

    def _pt(self, i):                                  # ellipse: centre, major, minor point
        if not hasattr(self, "_pts"):
            self._pts = [SketchPoint(self.sm.model, q) for q in self.edge[1:4]]
        return self._pts[i]

    def GetCenterPoint2(self):
        if self.edge[0] == "circle":
            if not hasattr(self, "_c"):
                self._c = SketchPoint(self.sm.model, self.edge[1])
            return self._c
        if self.edge[0] == "arc":                      # (start, mid, end): the circle through them
            if not hasattr(self, "_c"):
                p, q, r = (np.array(v, float) for v in self.edge[1:4])
                a = 2 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1]))
                c = (((p @ p) * (q[1] - r[1]) + (q @ q) * (r[1] - p[1]) + (r @ r) * (p[1] - q[1])) / a,
                     ((p @ p) * (r[0] - q[0]) + (q @ q) * (p[0] - r[0]) + (r @ r) * (q[0] - p[0])) / a)
                self._c = SketchPoint(self.sm.model, c)
            return self._c
        return self._pt(0)

    def GetMajorPoint2(self):
        return self._pt(1)

    def GetMinorPoint2(self):
        return self._pt(2)


class SketchPoint:
    def __init__(self, model, xy):
        self.model, self.X, self.Y, self.Z = model, float(xy[0]), float(xy[1]), 0.0

    def Select4(self, append, data):
        if not append:
            self.model.sel_list = []
        self.model.sel_list.append(("point", (self.X, self.Y)))
        return True


FACE_TWIST = float(os.environ.get("SWMOCK_FACE_TWIST", "0"))   # extra turn (deg) of face sketch axes


def face_frame(n):
    """Sketch axes on a planar face with outward normal n, as seen in SolidWorks 2025 for the
    face x = const (normal +x): sketch y along the projection of +y, sketch x = y x n (-z there);
    SWMOCK_FACE_TWIST turns them, so a test can check that scripts read the transform."""
    n = np.asarray(n, float)
    up = np.array([0.0, 1.0, 0.0]) if abs(n[1]) < 0.9 else np.array([0.0, 0.0, -1.0])
    v = up - (up @ n) * n
    v /= np.linalg.norm(v)
    u = np.cross(v, n)
    t = math.radians(FACE_TWIST)
    u, v = math.cos(t) * u + math.sin(t) * v, -math.sin(t) * u + math.cos(t) * v
    return tuple(u), tuple(v), tuple(n)


class Sketch:
    def __init__(self, plane, frame=None, offset=0.0, origin=None):
        self.plane = plane
        self.frame = FRAMES[plane] if frame is None else frame
        self.offset = offset                       # m, plane position along its normal
        self.origin = None if origin is None else np.asarray(origin, float)   # m, a general plane's origin
        self.ModelToSketchTransform = Transform(self.frame, offset, self.origin)
        self.edges = []

    def GetSketchPoints2(self):
        """The sketch's points (SolidWorks: shared end points once, arc and ellipse centres, an ellipse's axis
        points, points made by CreatePoint), as selectable points."""
        m = getattr(self, "model", None)
        out, seen = [], []

        def add(q):
            q = (float(q[0]), float(q[1]))
            if not any(abs(q[0] - a) < 1e-12 and abs(q[1] - b) < 1e-12 for a, b in seen):
                seen.append(q)
                out.append(SketchPoint(m, q))
        for seg in self.edges:
            if seg.deleted:
                continue
            e = seg.edge
            if e[0] == "line":
                add(e[1]); add(e[2])
            elif e[0] == "arc":
                p_, q_, r_ = (np.array(v, float) for v in e[1:4])
                a = 2 * (p_[0] * (q_[1] - r_[1]) + q_[0] * (r_[1] - p_[1]) + r_[0] * (p_[1] - q_[1]))
                c = (((p_ @ p_) * (q_[1] - r_[1]) + (q_ @ q_) * (r_[1] - p_[1]) + (r_ @ r_) * (p_[1] - q_[1])) / a,
                     ((p_ @ p_) * (r_[0] - q_[0]) + (q_ @ q_) * (p_[0] - r_[0]) + (r_ @ r_) * (q_[0] - p_[0])) / a)
                add(e[1]); add(e[3]); add(c)
            elif e[0] == "ellipse":
                add(e[1]); add(e[2]); add(e[3])
            elif e[0] == "earc":                       # centre, axis points, ends
                add(e[1]); add(e[2]); add(e[3]); add(e[4]); add(e[5])
            elif e[0] == "circle":
                add(e[1])
        for q in getattr(self, "points", []):         # CreatePoint: its own point even where an end point is
            out.append(SketchPoint(m, q))
        return out

    def to_model(self, q):
        """Sketch point (m) -> model point (m)."""
        u, v, n = (np.array(a, float) for a in self.frame)
        o = self.origin if self.origin is not None else self.offset * n
        return q[0] * u + q[1] * v + o

    def wire(self):
        """The sketch's edges as one CadQuery wire in model coordinates (mm)."""
        u, v, n = (np.array(a, float) for a in self.frame)
        P = lambda q: cq.Vector(*(K * self.to_model(q)))
        es = []
        for seg in self.edges:
            if seg.deleted or seg.ConstructionGeometry:
                continue
            e = seg.edge
            if e[0] == "line":
                es.append(cq.Edge.makeLine(P(e[1]), P(e[2])))
            elif e[0] == "arc":
                es.append(cq.Edge.makeThreePointArc(P(e[1]), P(e[2]), P(e[3])))
            elif e[0] == "earc":
                es.append(earc_edge(e, P, n))
            elif e[0] == "ellipse":
                c, pa, pb = (np.array(q, float) for q in e[1:4])
                a, b = np.linalg.norm(pa - c), np.linalg.norm(pb - c)
                es.append(cq.Edge.makeEllipse(K * a, K * b, P(c), cq.Vector(*n), P(pa) - P(c)))
            else:
                es.append(cq.Edge.makeCircle(K * e[2], P(e[1]), cq.Vector(*n)))
        return cq.Wire.assembleEdges(es)


class Sketch3D:
    """A 3D sketch (Insert3DSketch): entities in model coordinates (m); its transform is the identity."""
    def __init__(self):
        self.plane, self.frame, self.offset, self.origin = "3D", ((1, 0, 0), (0, 1, 0), (0, 0, 1)), 0.0, None
        self.ModelToSketchTransform = Transform(self.frame)
        self.edges, self.points3 = [], []
        self.is3d = True

    def to_model(self, q):
        return np.asarray(q, float)

    def GetSketchPoints2(self):
        m = getattr(self, "model", None)
        out = []
        for seg in self.edges:
            if not seg.deleted:
                out += [SketchPoint3(m, seg.edge[1]), SketchPoint3(m, seg.edge[2])]
        out += [SketchPoint3(m, q) for q in self.points3]
        return out

    def wire(self):
        es = [cq.Edge.makeLine(cq.Vector(*(K * np.asarray(seg.edge[1]))), cq.Vector(*(K * np.asarray(seg.edge[2]))))
              for seg in self.edges if not seg.deleted and not seg.ConstructionGeometry]
        return cq.Wire.assembleEdges(es)


class SketchPoint3:
    """A point of an open 3D sketch (model coordinates, m)."""
    def __init__(self, model, xyz):
        self.model = model
        self.X, self.Y, self.Z = (float(v) for v in xyz)

    def Select4(self, append, data):
        if not append:
            self.model.sel_list = []
        self.model.sel_list.append(("point", (self.X, self.Y, self.Z)))
        return True


class SkPoint3:
    """A point made by CreatePoint in a 3D sketch; once the sketch is closed, Select4 selects it as an external point."""
    def __init__(self, sm, xyz):
        self.sm, self.sk = sm, sm._sk
        self.X, self.Y, self.Z = (float(v) for v in xyz)

    def Select4(self, append, data):
        m = self.sm.model
        if self.sm._sk is self.sk:                     # its own sketch is open: an ordinary sketch point
            return SketchPoint3(m, (self.X, self.Y, self.Z)).Select4(append, data)
        if not append:
            m.ext_sel = []
        m.ext_sel.append(("point", np.array([self.X, self.Y, self.Z]), getattr(data, "Mark", 0)))
        return True


class Feature:
    def __init__(self, model, kind, data=None):
        self.model, self.kind, self.data, self.Name = model, kind, data, kind

    def Select2(self, append, mark):
        self.model.selected = self
        if not append:
            self.model.ext_sel = []
            self.model.ext_segs = []
        self.model.ext_sel.append(("feature", self, mark))
        if self.kind in ("Extrusion", "Revolve", "Fillet"):
            sl = getattr(self.model, "sel_list", None)
            if sl is None or not append:
                self.model.sel_list = sl = []
            sl.append(("feature", self))
        return True

    def GetTypeName2(self):
        return self.kind

    def GetSpecificFeature2(self):
        class S:
            def GetConstrainedStatus(self):
                return 1
        return S()

    def GetNextFeature(self):
        fs = self.model.features
        i = fs.index(self)
        return fs[i + 1] if i + 1 < len(fs) else None


ADDTODB_OPEN = os.environ.get("SWMOCK_ADDTODB_OPEN", "0") == "1"   # AddToDB sketches do not extrude
FRONT_ONLY = os.environ.get("SWMOCK_FRONT_ONLY", "1") == "1"   # sketches land on the Front Plane (seen)
LINES_NEED_ENTITY = os.environ.get("SWMOCK_LINES_NEED_ENTITY", "1") == "1"   # as seen in SolidWorks 2025
LAZY = int(os.environ.get("SWMOCK_LAZY_SKETCH", "0"))   # first N creation calls in a sketch fail
HIDE_SKETCH = os.environ.get("SWMOCK_HIDE_SKETCH", "0") == "1"   # ActiveSketch reads as None


class SketchManager:
    def __init__(self, model):
        self.model = model
        self._sk = None
        self.AddToDB = False
        self.DisplayWhenAdded = True

    @property
    def ActiveSketch(self):
        return None if HIDE_SKETCH else self._sk
        self.AddToDB = False
        self.DisplayWhenAdded = True

    def InsertSketch(self, flag):
        m = self.model
        if self._sk is None:
            assert m.selected is not None and m.selected.kind in ("RefPlane", "Face"), "no plane selected"
            if m.selected.kind == "Face":              # a face sketch keeps its face (seen)
                fr, off = m.selected.data
                self._sk = Sketch("face", fr, off)
            elif isinstance(m.selected.data, dict):    # a plane made by InsertRefPlane keeps its plane (seen)
                self._sk = Sketch("refplane", m.selected.data["frame"], origin=m.selected.data["origin"])
            else:
                self._sk = Sketch("Front Plane" if FRONT_ONLY else m.selected.data)
            self._sk.open = False
            self._sk.model = m
            self._lazy = LAZY
            m.selected = None
        else:
            m.features.append(Feature(m, "ProfileFeature", self._sk))
            self._sk = None

    def Insert3DSketch(self, flag):
        m = self.model
        if self._sk is None:
            self._sk = Sketch3D()
            self._sk.open = False
            self._sk.model = m
            self._lazy = 0
            m.selected = None
        else:
            m.features.append(Feature(m, "ProfileFeature", self._sk))
            self._sk = None

    def CreateLine(self, x1, y1, z1, x2, y2, z2):
        if getattr(self._sk, "is3d", False):
            return self._add(("line", (x1, y1, z1), (x2, y2, z2)))
        if self._lazy > 0:
            self._lazy -= 1
            return None
        self._sk.open |= ADDTODB_OPEN and self.AddToDB
        if LINES_NEED_ENTITY and not self._sk.edges:
            return None
        return self._add(("line", (x1, y1), (x2, y2)))

    def CreateCircleByRadius(self, xc, yc, zc, r):
        if self._lazy > 0:
            self._lazy -= 1
            return None
        return self._add(("circle", (xc, yc), r))

    def CreateArc(self, xc, yc, zc, x1, y1, z1, x2, y2, z2, d):
        if self._lazy > 0:
            self._lazy -= 1
            return None
        a1 = math.atan2(y1 - yc, x1 - xc)
        a2 = math.atan2(y2 - yc, x2 - xc)
        sw = (a2 - a1) % (2 * math.pi)
        if d < 0:
            sw -= 2 * math.pi
        r = math.hypot(x1 - xc, y1 - yc)
        am = a1 + sw / 2
        return self._add(("arc", (x1, y1), (xc + r * math.cos(am), yc + r * math.sin(am)), (x2, y2)))

    def CreateCenterLine(self, x1, y1, z1, x2, y2, z2):
        seg = self._add(("line", (x1, y1), (x2, y2)))
        seg.ConstructionGeometry = True
        seg.centerline = True
        return seg

    def CreateCircularSketchStepAndRepeat(self, R, A, n, spacing, rotate, delete, *dims):
        """Copies of the selected segments turned about a centre at the seeds' bounding-box centre
        plus R (cos A, sin A); A outside 0..2 pi is refused (seen in SolidWorks 2025)."""
        seeds = [seg for kind, seg in self.model.sel_list if isinstance(seg, Segment)]
        if not seeds or not (0.0 <= A < 2 * math.pi):
            return False
        xs, ys = [], []
        for seg in seeds:
            e = seg.edge
            pts = [np.array(q, float) for q in e[1:]]
            xs += [q[0] for q in pts]
            ys += [q[1] for q in pts]
            if e[0] == "arc":
                p, q, r = pts
                den = 2 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1]))
                c = np.array([((p @ p) * (q[1] - r[1]) + (q @ q) * (r[1] - p[1]) + (r @ r) * (p[1] - q[1])) / den,
                              ((p @ p) * (r[0] - q[0]) + (q @ q) * (p[0] - r[0]) + (r @ r) * (q[0] - p[0])) / den])
                rad = np.linalg.norm(p - c)
                ang = lambda v: math.atan2(v[1] - c[1], v[0] - c[0])
                a0, a1, am = ang(p), ang(r), ang(q)
                sw = (a1 - a0) % (2 * math.pi)
                ccw = (am - a0) % (2 * math.pi) < sw
                for k in range(4):
                    t = k * math.pi / 2
                    d = (t - a0) % (2 * math.pi) if ccw else (a0 - t) % (2 * math.pi)
                    if d < (sw if ccw else 2 * math.pi - sw):
                        xs.append(c[0] + rad * math.cos(t))
                        ys.append(c[1] + rad * math.sin(t))
        bc = np.array([(max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2])
        cen = bc + R * np.array([math.cos(A), math.sin(A)])
        for m in range(1, n):
            t = m * spacing
            rot = lambda q: tuple(cen + np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]]) @ (np.array(q) - cen))
            for seg in seeds:
                e = seg.edge
                self._add((e[0],) + tuple(rot(q) for q in e[1:]))
        return True

    def CreateLinearSketchStepAndRepeat(self, nx, ny, sx, sy, ax, ay, delete, *dims):
        """Copies of the selected segments (seen in SolidWorks 2025: copies share end points)."""
        seeds = [seg for kind, seg in self.model.sel_list if isinstance(seg, Segment)]
        if not seeds or ny != 1:
            return False
        dx, dy = sx * math.cos(ax), sx * math.sin(ax)
        for m in range(1, nx):
            for seg in seeds:
                e = seg.edge
                self._add((e[0],) + tuple((q[0] + m * dx, q[1] + m * dy) for q in e[1:]))
        return True

    def CreateEllipse(self, xc, yc, zc, xa, ya, za, xb, yb, zb):
        if self._lazy > 0:
            self._lazy -= 1
            return None
        return self._add(("ellipse", (xc, yc), (xa, ya), (xb, yb)))

    def _add(self, edge):
        seg = Segment(self, edge)
        seg.sk = self._sk
        self._sk.edges.append(seg)
        return seg

    def CreatePoint(self, x, y, z):
        if getattr(self._sk, "is3d", False):
            self._sk.points3.append((x, y, z))
            return SkPoint3(self, (x, y, z))
        self._sk.points = getattr(self._sk, "points", []) + [(x, y)]
        return SkPoint(self, (x, y))

    def CreateEllipticalArc(self, xc, yc, zc, xa, ya, za, xb, yb, zb, x1, y1, z1, x2, y2, z2, d):
        return self._add(("earc", (xc, yc), (xa, ya), (xb, yb), (x1, y1), (x2, y2), d))


class SelectData:
    Mark = 0


class Body:
    def __init__(self, model):
        self.model = model

    def Select2(self, append, data):
        self.model.selected = self
        return True


class MassProperty:
    def __init__(self, v, c):
        self.Volume, self.CenterOfMass = v, c


class Extension:
    def __init__(self, model):
        self.model = model
        self._mark = 0

    def SelectByID2(self, name, typ, x, y, z, append, mark, callout, opt):
        m = self.model
        self._mark = mark
        if not append:
            m.sel_list = []
        if typ in ("SKETCHSEGMENT", "SKETCHPOINT") and m.SketchManager._sk is not None:
            P = np.array([x, y])
            for seg in m.SketchManager._sk.edges:
                e = seg.edge
                if seg.deleted or seg.ConstructionGeometry or e[0] == "ellipse":
                    continue
                if e[0] == "circle":
                    c, r = np.array(e[1]), e[2]
                    d = np.linalg.norm(P - c)
                    if typ == "SKETCHSEGMENT" and abs(d - r) < 1e-7:
                        m.sel_list.append(("circle", seg))
                        return True
                    if typ == "SKETCHPOINT" and d < 1e-7:
                        m.sel_list.append(("point", (x, y)))
                        return True
                    continue
                a, b = np.array(e[1]), np.array(e[-1])
                if typ == "SKETCHPOINT" and min(np.linalg.norm(P - a), np.linalg.norm(P - b)) < 1e-7:
                    m.sel_list.append(("point", (x, y)))
                    return True
                if typ == "SKETCHSEGMENT":
                    if e[0] == "line":
                        d = b - a
                        t = (P - a) @ d / (d @ d)
                        if 0 < t < 1 and np.linalg.norm(a + t * d - P) < 1e-7:
                            m.sel_list.append(("line", seg))
                            return True
                    else:                                               # arc picked at any point on it
                        p, q, r = (np.array(v) for v in e[1:])
                        den = 2 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1]))
                        c = np.array([((p @ p) * (q[1] - r[1]) + (q @ q) * (r[1] - p[1]) + (r @ r) * (p[1] - q[1])) / den,
                                      ((p @ p) * (r[0] - q[0]) + (q @ q) * (p[0] - r[0]) + (r @ r) * (q[0] - p[0])) / den])
                        rad = np.linalg.norm(p - c)
                        ang = lambda v: math.atan2(v[1] - c[1], v[0] - c[0])
                        a0, a1, am, ax = ang(p), ang(r), ang(q), ang(P)
                        s1, sm, sx = (a1 - a0) % (2 * math.pi), (am - a0) % (2 * math.pi), (ax - a0) % (2 * math.pi)
                        inside = (sx <= s1) if sm <= s1 else (sx >= s1)
                        if abs(np.linalg.norm(P - c) - rad) < 1e-7 and inside:
                            m.sel_list.append(("arc", seg))
                            return True
            return False
        if typ == "FACE" and m.body is not None:
            # as seen in SolidWorks 2025: the first face along the line of sight through the point
            # (front view, looking along -z), not necessarily the face containing it
            v = np.array([0.0, 0.0, 1.0])
            return self._ray_face(np.array([x, y, z]) * K + 1000.0 * v, -v)
        if typ == "EDGE" and m.body is not None:
            if os.environ.get("SWMOCK_EDGE_SIGHT") == "1":
                # as seen in SolidWorks (TO junction fillets): the first edge along the line of sight
                # through the point (front view, looking along -z), which may lie in front of it
                v = np.array([0.0, 0.0, 1.0])
                return self._ray_edge(np.array([x, y, z]) * K + 1000.0 * v, -v, 1e-6)
            v = cq.Vertex.makeVertex(*(np.array([x, y, z]) * K))     # the body edge through the point
            best = min(m.body.Edges(), key=lambda e: e.distance(v))
            if best.distance(v) > 1e-6:
                return False
            m.sel_list.append(("edge", best))
            return True
        return self._select_other(name, typ, x, y, z, append)

    def SelectByRay(self, x, y, z, dx, dy, dz, radius, typ, append, mark, opt):
        if typ not in (1, 2) or self.model.body is None:         # swSelEDGES, swSelFACES
            return False
        self.model.last_ray = (np.array([x, y, z]) * K, np.array([dx, dy, dz], float))
        d = np.array([dx, dy, dz], float)
        if typ == 1:
            if not append:
                self.model.sel_list = []
            return self._ray_edge(np.array([x, y, z]) * K, d / np.linalg.norm(d), max(radius * K, 1e-6))
        if self._ray_face(np.array([x, y, z]) * K, d / np.linalg.norm(d)):
            return True
        v = cq.Vertex.makeVertex(*(np.array([x, y, z]) * K))     # a curved face just ahead (pattern axis)
        return any(fc.geomType() != "PLANE" and fc.distance(v) < 1.0 for fc in self.model.body.Faces())

    def _ray_edge(self, o, d, tol):
        """Select the body edge nearest to o along the ray o + t d (t > 0) that passes within tol of it."""
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        m = self.model
        ray = cq.Edge.makeLine(cq.Vector(*o), cq.Vector(*(o + 2000.0 * d)))
        best = None
        for e in m.body.Edges():
            ex = BRepExtrema_DistShapeShape(e.wrapped, ray.wrapped)
            if not ex.IsDone() or ex.Value() > tol:
                continue
            q = ex.PointOnShape2(1)
            t = float((np.array([q.X(), q.Y(), q.Z()]) - o) @ d)
            if best is None or t < best[0]:
                best = (t, e)
        if best is None:
            return False
        m.sel_list.append(("edge", best[1]))
        return True

    def _ray_face(self, p0, d):
        """Select the first planar face hit by the ray p0 + t d (mm)."""
        m = self.model
        best = None
        for fc in m.body.Faces():
            if fc.geomType() != "PLANE":
                continue
            c = np.array(fc.Center().toTuple())
            nrm = np.array(fc.normalAt(cq.Vector(*c)).toTuple())
            den = d @ nrm
            if abs(den) < 1e-12:
                continue
            t = ((c - p0) @ nrm) / den
            if t <= 0:
                continue
            q = p0 + t * d
            if cq.Vertex.makeVertex(*q).distance(fc) < 1e-6 and (best is None or t < best[0]):
                best = (t, fc, nrm, c)
        if best is None:
            return False
        _, fc, nrm, c = best
        m.selected = Feature(m, "Face", (face_frame(nrm), float(nrm @ c) / K))   # the face's own plane
        m.selected_normal = tuple(nrm)
        return True

    def _select_other(self, name, typ, x, y, z, append):
        m = self.model
        if not append:
            m.ext_sel = []
        if typ in ("EXTSKETCHSEGMENT", "EXTSKETCHPOINT") and name != "Point1@Origin":
            P = np.array([x, y, z])
            for f in reversed(m.features):
                if f.kind != "ProfileFeature":
                    continue
                for seg in f.data.edges:
                    if seg.deleted or seg.ConstructionGeometry or seg.edge[0] not in ("line", "arc"):
                        continue
                    e = seg.edge
                    ends = [f.data.to_model(e[1]), f.data.to_model(e[-1])]
                    if typ == "EXTSKETCHPOINT":
                        if min(np.linalg.norm(P - q) for q in ends) < 1e-7:
                            m.ext_sel.append(("point", P, self._mark))
                            return True
                        continue
                    pts = [f.data.to_model(q) for q in e[1:]]
                    edge = (cq.Edge.makeLine(cq.Vector(*(K * pts[0])), cq.Vector(*(K * pts[1]))) if e[0] == "line" else
                            cq.Edge.makeThreePointArc(*(cq.Vector(*(K * q)) for q in pts)))
                    if edge.distance(cq.Vertex.makeVertex(*(K * P))) < 1e-4:
                        m.ext_sel.append(("segment", edge, self._mark))
                        return True
            return False
        if typ == "EXTSKETCHPOINT" and name == "Point1@Origin":
            m.sel_list.append(("point", (0.0, 0.0)))
            return True
        for f in self.model.features:
            if f.kind == "RefPlane" and (f.data == name or f.Name == name):
                self.model.selected = f
                m.ext_sel.append(("feature", f, self._mark))
                return True
        return False

    def CreateMassProperty(self):
        b = self.model.body
        if b is None:
            return MassProperty(0.0, (0.0, 0.0, 0.0))
        frame = any(f.kind == "ProfileFeature" and getattr(f.data, "is3d", False) for f in self.model.features)
        if frame or any(f.geomType() == "BSPLINE" for f in b.Faces()):
            # lofted members, or a 3D frame (swframe: a 3D skeleton sketch, then many spheres and cones):
            # CadQuery's integration, as the reference; stl2cad's fine-tessellation check does not finish on
            # lofted faces and runs out of memory on a frame
            cc = b.Center()
            v, c = b.Volume(), (cc.x, cc.y, cc.z)
        else:
            from stl2cad.core import mass_props      # accurate, as SolidWorks
            v, c = mass_props(b)
        return MassProperty(v / K ** 3, tuple(x / K for x in c))   # m^3, m

    def DeleteSelection2(self, opt):
        f = self.model.selected
        assert f is self.model.features[-1] and f.kind in ("Extrusion", "ProfileFeature", "MoveCopyBody", "Revolve",
                                                           "Fillet", "Loft", "RefPlane")
        self.model.features.pop()
        if f.kind in ("Extrusion", "MoveCopyBody", "Revolve", "Fillet", "Loft"):
            self.model.body = f.data["before"]
        return True


def _loops(edges):
    """Closed wires of the edges, grouped by shared end points (Wire.combine joined edges of
    separate contours into one wire with a 53 mm vertex tolerance)."""
    key = lambda p: tuple(round(x, 6) for x in p.toTuple())
    ends = [(key(e.startPoint()), key(e.endPoint())) for e in edges]
    left, out = set(range(len(edges))), []
    while left:
        comp, todo = set(), [left.pop()]
        while todo:
            i = todo.pop()
            comp.add(i)
            for j in list(left):
                if set(ends[i]) & set(ends[j]):
                    left.discard(j)
                    todo.append(j)
        out.append(cq.Wire.assembleEdges([edges[i] for i in sorted(comp)]))
    return out


SPHERE_AXIS = (0.2672612, 0.5345225, 0.8017837)


def _sphere(es, a, b):
    """A semicircle closed by its diameter, revolved a full turn about that diameter: the sphere, built with a
    generic axis (frame3d's nodes). OCC's Booleans of a revolved sphere (poles on the revolve axis) with a member
    cone meeting it near a pole dropped the overlap. None for any other profile."""
    if len(es) != 2:
        return None
    arc = next((e for e in es if e.geomType() == "CIRCLE"), None)
    ln = next((e for e in es if e.geomType() == "LINE"), None)
    if arc is None or ln is None:
        return None
    c, r = arc.arcCenter(), arc.radius()
    p, q = ln.startPoint(), ln.endPoint()
    ax = b - a
    on_axis = lambda v: (v - a).cross(ax).Length < 1e-6 * max(1.0, ax.Length)
    if not (on_axis(c) and on_axis(p) and on_axis(q) and abs((p - q).Length - 2 * r) < 1e-6
            and abs(arc.Length() - math.pi * r) < 1e-6):
        return None
    return cq.Solid.makeSphere(r, c, cq.Vector(*SPHERE_AXIS), angleDegrees1=-90, angleDegrees2=90)


def _cone(profs):
    """Two circle sections on parallel planes with their centres on the common normal: the loft is a cone (or a
    cylinder), built as one (OCC's ruled loft of circles is a B-spline that later Booleans failed on); None
    otherwise."""
    if len(profs) != 2:
        return None
    es = [w.Edges() for w in profs]
    if any(len(e) != 1 or e[0].geomType() != "CIRCLE" for e in es):
        return None
    (a,), (b,) = es
    ca, cb = np.array(a.arcCenter().toTuple()), np.array(b.arcCenter().toTuple())
    na = np.array(a.normal().toTuple()); nb = np.array(b.normal().toTuple())
    d = cb - ca
    h = float(np.linalg.norm(d))
    if h < 1e-9 or abs(abs(na @ nb) - 1) > 1e-9 or np.linalg.norm(np.cross(d / h, na)) > 1e-9:
        return None
    ra, rb = a.radius(), b.radius()
    if abs(ra - rb) < 1e-12:
        return cq.Solid.makeCylinder(ra, h, cq.Vector(*ca), cq.Vector(*(d / h)))
    return cq.Solid.makeCone(ra, rb, h, cq.Vector(*ca), cq.Vector(*(d / h)))


def _turn(shape, axis, ang):
    """Shape turned by ang (rad) about the axis through the origin; a quarter turn is applied as an
    exact matrix (a rotation by angle leaves 1e-16 tilts, and OCC booleans then fail between
    coaxial cylinders)."""
    q = ang / (math.pi / 2)
    if abs(q - round(q)) > 1e-12:
        return shape.rotate(cq.Vector(0, 0, 0), cq.Vector(*axis), math.degrees(ang))
    from OCP.gp import gp_Trsf
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    c, s_ = [1, 0, -1, 0][int(round(q)) % 4], [0, 1, 0, -1][int(round(q)) % 4]
    x, y, z = axis
    K_ = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], float)
    R = np.eye(3) + s_ * K_ + (1 - c) * (K_ @ K_)
    t = gp_Trsf()
    t.SetValues(*R[0], 0.0, *R[1], 0.0, *R[2], 0.0)
    return cq.Shape.cast(BRepBuilderAPI_Transform(shape.wrapped, t, True).Shape())


def _boolean(a, b, cut):
    """a - b or a + b. OCC booleans with near-coincident faces can return an invalid solid that
    later booleans get wrong (seen: a boss whose circle is a flange's end arc, then a cut beside
    the body removing material); a result that is invalid or breaks the volume identity
    vol(a) +- vol(b) -+ vol(a & b) is recomputed with a fuzzy tolerance."""
    def run(tol, clean=True, glue=False):
        r = a.cut(b, tol=tol) if cut else a.fuse(b, tol=tol, glue=glue)
        return r.clean() if clean else r
    want = None
    for t_ in (None, 1e-5, 1e-4):                  # the identity's common part, fuzzy if need be
        try:
            common = a.intersect(b, tol=t_).Volume()
            want = a.Volume() - common if cut else a.Volume() + b.Volume() - common
            break
        except Exception:
            continue

    def good(r):
        try:
            return r.isValid() and (want is None or abs(r.Volume() - want) <= 1e-6 * max(1.0, abs(want)))
        except Exception:
            return False
    try:
        r = run(None)
    except Exception:
        r = None
    if r is not None and good(r):
        return r
    for tol, cl, gl in ((None, False, False), (None, True, True), (1e-5, True, False), (1e-5, False, False),
                        (1e-4, True, False), (1e-3, True, False)):
        try:
            r2 = run(tol, cl, gl and not cut)  # near-coincident faces (a boss on a flange's end arc)
        except Exception:
            continue
        if good(r2):
            return r2
        r2 = r2.fix()
        if good(r2):
            return r2
        r = r if r is not None else r2
    return r


class FeatureManager:
    def __init__(self, model):
        self.model = model

    def InsertRefPlane(self, c1, d1, c2, d2, c3, d3):
        """Two cases: a plane at distance d1 from a selected plane (c1 = 8, swRefPlaneReferenceConstraint_Distance;
        256 added: the other side), and a plane normal to a selected sketch segment through a selected point
        (c1 = 2 Perpendicular, c2 = 4 Coincident; normal along the curve's direction at the point, sketch axes
        as on faces: seen in sw_probe_sweep)."""
        m = self.model
        sel = getattr(m, "ext_sel", [])
        if c1 & 8:
            base = next((o for k, o, _ in sel if k == "feature" and o.kind == "RefPlane"), None)
            if base is None:
                return None
            if isinstance(base.data, dict):
                fr, o = base.data["frame"], base.data["origin"]
            else:
                fr, o = FRAMES[base.data], np.zeros(3)
            n = np.array(fr[2], float)
            data = dict(frame=fr, origin=o + (-1 if c1 & 256 else 1) * d1 * n)
        elif c1 == 1 and c2 == 4:                      # parallel to a selected plane, through a selected point
            base = next((o for k, o, _ in sel if k == "feature" and o.kind == "RefPlane"), None)
            pt = next((o for k, o, mk in sel if k == "point"), None)
            if base is None or pt is None:
                return None
            fr = base.data["frame"] if isinstance(base.data, dict) else FRAMES[base.data]
            data = dict(frame=fr, origin=np.asarray(pt, float))
        elif c1 == 2 and c2 == 4:
            seg = next((o for k, o, mk in sel if k == "segment"), None)
            pt = next((o for k, o, mk in sel if k == "point"), None)
            if seg is None or pt is None:
                return None
            v = cq.Vertex.makeVertex(*(K * pt))
            ps = [seg.positionAt(t) for t in np.linspace(0, 1, 201)]
            t = np.linspace(0, 1, 201)[int(np.argmin([(p - v.Center()).Length for p in ps]))]
            d = np.array(seg.tangentAt(t).toTuple()); d /= np.linalg.norm(d)
            data = dict(frame=face_frame(d), origin=np.asarray(pt, float))
        else:
            return None
        f = Feature(m, "RefPlane", data)
        f.Name = "Plane%d" % (sum(1 for g in m.features if g.kind == "RefPlane") - 2)
        m.features.append(f)
        return f

    def InsertProtrusionBlend2(self, *args):
        """Loft through the selected sketches (Mark 1, in order) with a centerline sketch (Mark 4): the
        sections swept along the centerline (CadQuery's sweep_multi; within 0.04 % of SolidWorks' loft in
        sw_probe_sweep)."""
        assert len(args) == 18, len(args)
        m = self.model
        sel = getattr(m, "ext_sel", [])
        profs = [o.data.wire() for k, o, mk in sel if k == "feature" and o.kind == "ProfileFeature" and mk == 1]
        cl = [o.data.wire() for k, o, mk in sel if k == "feature" and o.kind == "ProfileFeature" and mk == 4]
        if len(profs) < 2 or len(cl) > 1:
            return None
        if cl:
            body = cq.Solid.sweep_multi(profs, cl[0], True, False)
        else:                                # straight loft between the sections
            body = _cone(profs) or cq.Solid.makeLoft(profs, True)
        before = m.body
        if before is None:
            m.body = body
        elif not cl:                         # a frame member (swframe): fused with the volume-identity check
            m.body = _boolean(before, body, cut=False)
        else:                                # as the reference (swsweep.steps): plain fuse, fuzzy if that fails
            def fuse(tol):
                r = before.fuse(body, tol=tol) if tol else before.fuse(body)
                try:
                    return r.clean()
                except Exception:            # OCC's clean failed (a lofted frame member): the fuse as it is
                    return r
            u = fuse(None)
            if not (u.isValid() and len(u.Solids()) <= len(before.Solids()) + 1):
                u = fuse(1e-4)
            m.body = u
        f = Feature(m, "Loft", dict(before=before))
        m.features.append(f)
        return f

    def _solid(self, sk, rev, h, t0, off, flip, cut, draft=None):
        u, v, n = (np.array(a, float) for a in sk.frame)
        o = sk.origin if getattr(sk, "origin", None) is not None else sk.offset * n
        seen = []

        def canon(q):                        # end points within 1e-9 m are one point (pattern
            q = np.asarray(q, float)         # copies are turned in floating point)
            for w in seen:
                if np.linalg.norm(w - q) < 1e-9:
                    return w
            seen.append(q)
            return q
        P = lambda q: cq.Vector(*(K * (canon(q)[0] * u + canon(q)[1] * v + o)))   # geometry kept in mm
        es = []
        for seg in sk.edges:
            if seg.deleted or seg.ConstructionGeometry:
                continue
            e = seg.edge
            if e[0] == "line":
                es.append(cq.Edge.makeLine(P(e[1]), P(e[2])))
            elif e[0] == "arc":
                es.append(cq.Edge.makeThreePointArc(P(e[1]), P(e[2]), P(e[3])))
            elif e[0] == "earc":
                es.append(earc_edge(e, P, n))
            elif e[0] == "ellipse":
                c, pa, pb = (np.array(q, float) for q in e[1:4])
                a, b = np.linalg.norm(pa - c), np.linalg.norm(pb - c)
                es.append(cq.Edge.makeEllipse(K * a, K * b, P(c), cq.Vector(*n), P(pa) - P(c)))
            else:
                es.append(cq.Edge.makeCircle(K * e[2], P(e[1]), cq.Vector(*n)))
        wires = cq.Wire.combine(es)
        start = K * (off if t0 == 3 else 0.0) * (-1 if flip else 1)
        h = K * h
        d = -1.0 if (CUT_REVERSED if cut else BOSS_REVERSED) else 1.0
        if rev:
            d = -d
        solids = []
        if draft is not None:                        # (angle rad, outward): a tapered extrusion
            ang, outward = draft
            out = None
            for w in wires:
                es = w.Edges()
                if len(es) == 1 and es[0].geomType() == "CIRCLE":
                    # a circle tapered to (nearly) a point is a cone (OCC's tapered prism
                    # crashes there)
                    r0 = es[0].radius()
                    r1 = max(r0 + (1 if outward else -1) * h * math.tan(ang), 0.0)
                    c0 = np.array(es[0].arcCenter().toTuple()) + n * start
                    s = cq.Solid.makeCone(r0, r1, h, cq.Vector(*c0), cq.Vector(*(n * d)))
                else:
                    f = cq.Face.makeFromWires(w, []).translate(cq.Vector(*(n * start)))
                    s = cq.Solid.extrudeLinear(f, cq.Vector(*(n * d * h)),
                                               taper=(-1 if outward else 1) * math.degrees(ang))
                out = s if out is None else out.fuse(s)
            return out
        # nested contours as SolidWorks reads them: a contour inside an odd number of others is a
        # hole in the contour that directly contains it (even-odd nesting)
        from shapely.geometry import Polygon as _Poly
        def _poly(w):                        # the contour in the sketch plane, as a polygon
            pts = [e.positionAt(t).toTuple() for e in w.Edges() for t in np.linspace(0, 1, 16, endpoint=False)]
            return _Poly([(np.dot(p, u), np.dot(p, v)) for p in pts]).buffer(0)
        polys = [_poly(w) for w in wires]
        def inside(i, j):                    # contour i lies inside contour j
            return polys[i].area < polys[j].area and polys[j].contains(polys[i].representative_point())
        parents = [[j for j in range(len(wires)) if j != i and inside(i, j)] for i in range(len(wires))]
        for i, w in enumerate(wires):
            if len(parents[i]) % 2:
                continue                     # a hole: cut from the contour that directly holds it
            holes = [wires[k] for k in range(len(wires)) if i in parents[k] and len(parents[k]) == len(parents[i]) + 1]
            f = cq.Face.makeFromWires(w, holes)
            # always extruded along +n from the lower end: the same solid, but OCC fuses
            # coincident faces (a boss on a flange's end arc) reliably only when they are
            # parametrized alike
            s = cq.Solid.extrudeLinear(f, cq.Vector(*(n * h)))
            solids.append(s.translate(cq.Vector(*(n * (start if d > 0 else start - h)))))
        out = solids[0]
        for s in solids[1:]:
            out = out.fuse(s)
        return out

    def _apply(self, rev, h, t0, off, flip, cut, draft=None):
        m = self.model
        sk = m.selected.data
        if getattr(sk, "open", False):
            return None
        body = self._solid(sk, rev, h, t0, off, flip, cut, draft)
        before = m.body
        if cut:
            if before is None:
                return None
            m.body = _boolean(before, body, cut=True)
        else:
            m.body = body if before is None else _boolean(before, body, cut=False)
        f = Feature(m, "Extrusion", dict(before=before, tool=body, cut=cut))
        m.features.append(f)
        LAST["features"] = m.features
        return f

    def _pattern(self, n):
        """Circular pattern (n instances, equal spacing) of the selected Extrusion features about the
        axis of the cylindrical face nearest to the last ray's start."""
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.GeomAbs import GeomAbs_Cylinder
        m = self.model
        ray = getattr(m, "last_ray", None)
        feats = [f for kind, f in (getattr(m, "sel_list", None) or []) if kind == "feature" and f.kind == "Extrusion"]
        if ray is None or not feats or m.body is None:
            return None
        v0 = cq.Vertex.makeVertex(*ray[0])
        best = None
        for fc in m.body.Faces():
            ad = BRepAdaptor_Surface(fc.wrapped)
            if ad.GetType() != GeomAbs_Cylinder:
                continue
            dd = fc.distance(v0)
            if best is None or dd < best[0]:
                ax = ad.Cylinder().Axis()
                best = (dd, ax.Location(), ax.Direction())
        if best is None:
            return None
        loc = cq.Vector(best[1].X(), best[1].Y(), best[1].Z())
        dirn = cq.Vector(best[2].X(), best[2].Y(), best[2].Z())
        before = m.body
        for f in feats:
            for k in range(1, n):
                t = f.data["tool"].rotate(loc, loc + dirn, 360.0 * k / n)
                m.body = _boolean(m.body, t, cut=f.data["cut"])
        f = Feature(m, "Extrusion", dict(before=before, tool=None, cut=False))
        m.features.append(f)
        return f

    def FeatureCircularPattern4(self, n, spacing, flip, dname, geometry, equal, vary):
        return self._pattern(int(n))

    def FeatureRevolve2(self, single, solid, thin, cut, rev, both, t1, t2, a1, a2, *rest):
        """Revolve the selected sketch's profile about its single construction centerline."""
        m = self.model
        sk = m.selected.data
        u, v, n = (np.array(a, float) for a in sk.frame)
        o = sk.origin if getattr(sk, "origin", None) is not None else sk.offset * n
        seen = []

        def canon(q):                                # end points within 1e-9 m are one point
            q = np.asarray(q, float)
            for w in seen:
                if np.linalg.norm(w - q) < 1e-9:
                    return w
            seen.append(q)
            return q
        P = lambda q: cq.Vector(*np.round(K * (canon(q)[0] * u + canon(q)[1] * v + o), 9))
        axes = [seg for seg in sk.edges if not seg.deleted and seg.ConstructionGeometry and seg.edge[0] == "line"
                and getattr(seg, "centerline", False)]
        if not axes:                                 # a line of the profile selected as the axis (Mark 4 or 16)
            axes = [g for g, mk in getattr(m, "ext_segs", []) if mk in (4, 16) and g in sk.edges and g.edge[0] == "line"]
        if len(axes) != 1:
            return None
        es = []
        for seg in sk.edges:
            if seg.deleted or seg.ConstructionGeometry:
                continue
            e = seg.edge
            if e[0] == "line":
                es.append(cq.Edge.makeLine(P(e[1]), P(e[2])))
            elif e[0] == "arc":
                es.append(cq.Edge.makeThreePointArc(P(e[1]), P(e[2]), P(e[3])))
            elif e[0] == "earc":
                es.append(earc_edge(e, lambda q: cq.Vector(*np.round(K * (q[0] * u + q[1] * v + o), 9)), n))
            else:
                es.append(cq.Edge.makeCircle(K * e[2], P(e[1]), cq.Vector(*n)))
        wires = sorted(_loops(es), key=lambda w: -cq.Face.makeFromWires(w, []).Area())
        # disjoint contours are separate regions; a contour inside another is its hole
        outers = []
        for w in wires:
            f_ = cq.Face.makeFromWires(w, [])
            host = next((o for o in outers if o[0].distance(cq.Vertex.makeVertex(*f_.Center().toTuple())) < 1e-9), None)
            if host is None:
                outers.append((f_, w, []))
            else:
                host[2].append(w)
        a, b = axes[0].edge[1], axes[0].edge[2]
        before = m.body
        body = before
        sph = _sphere(es, P(a), P(b)) if abs(a1 - 2 * math.pi) < 1e-9 else None
        for _, w, holes in outers:                  # one region at a time (coincident faces)
            part = sph or cq.Solid.revolve(cq.Face.makeFromWires(w, holes), math.degrees(a1), P(a), P(b))
            body = part if body is None else _boolean(body, part, cut=bool(cut))
        m.body = body
        f = Feature(m, "Revolve", dict(before=before))
        m.features.append(f)
        LAST["features"] = m.features
        return f

    def FeatureExtrusion3(self, Sd, Flip, Dir, T1, T2, D1, D2, *rest):
        t0, off, flip = rest[-3], rest[-2], rest[-1]
        if not Sd:                                   # both directions from the sketch plane
            m = self.model
            sk = m.selected
            f = self._apply(Dir, D1, t0, off, flip, False)
            if f is None or D2 <= 0:
                return f
            m.features.pop()
            m.selected = sk
            g = self._apply(not Dir, D2, t0, off, flip, False)
            g.data["before"] = f.data["before"]
            return g
        return self._apply(Dir, D1, t0, off, flip, False)

    def InsertMoveCopyBody2(self, tx, ty, tz, td, px, py, pz, ax, ay, az, copy, n):
        m = self.model
        if not isinstance(m.selected, Body) or m.body is None:
            return None
        before = m.body
        b = before
        order = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        if os.environ.get("SWMOCK_ANGLE_ORDER", "zyx") == "zyx":   # observed: first angle about z
            order = ((0, 0, 1), (0, 1, 0), (1, 0, 0))
        for ang, axis in zip((ax, ay, az), order):
            if ang:
                sgn = -1 if os.environ.get("SWMOCK_ROT_SIGN") == "-1" else 1
                b = _turn(b, axis, sgn * ang)
        # translation: the vector itself (default), or a direction times TransDist
        if os.environ.get("SWMOCK_TRANS", "vector") == "unit":
            t = (tx * td, ty * td, tz * td)
        else:
            t = (tx, ty, tz)
        if any(t):
            b = b.translate(cq.Vector(*(v * K for v in t)))
        m.body = b
        f = Feature(m, "MoveCopyBody", dict(before=before))
        m.features.append(f)
        return f

    def _edge_feature(self, make):
        m = self.model
        es = [e for t, e in getattr(m, "sel_list", []) if t == "edge"]
        if not es or m.body is None:
            return None
        try:
            body = make(m.body, es)
        except Exception:
            return None
        if body is None or not body.isValid():
            return None
        before = m.body
        m.body = body
        f = Feature(m, "Fillet", dict(before=before))
        m.features.append(f)
        m.sel_list = []
        return f

    def FeatureFillet3(self, opts, r1, r2, rho, ftyp, overflow, conic, *arrays):
        return self._edge_feature(lambda b, es: b.fillet(r1 * K, es))

    def InsertFeatureChamfer(self, opts, typ, width, angle, other, *rest):
        if typ == 2:                                  # distance-distance
            return self._edge_feature(lambda b, es: b.chamfer(width * K, other * K if other else None, es))
        if typ == 1:                                  # angle-distance
            return self._edge_feature(lambda b, es: b.chamfer(width * K, width * K * math.tan(angle), es))
        return None

    def FeatureCut4(self, Sd, Flip, Dir, T1, T2, D1, D2, Dchk1, Dchk2, Ddir1, Ddir2, Dang1, *rest):
        t0, off, flip = rest[-4], rest[-3], rest[-2]
        return self._apply(Dir, D1, t0, off, flip, True, (Dang1, Ddir1) if Dchk1 else None)


class SelectionManager:
    def __init__(self, model):
        self.model = model

    def CreateSelectData(self):
        return SelectData()

    def GetSelectedObject6(self, i, mark):
        class F:
            pass
        f = F()
        f.Normal = getattr(self.model, "selected_normal", (0.0, 0.0, 0.0))
        return f

    def GetSelectedObjectCount2(self, mark):
        n = len(getattr(self.model, "sel_list", []) or [])
        return n if n else (0 if self.model.selected is None else 1)


class Model:
    def __init__(self):
        self.features = [Feature(self, "RefPlane", p) for p in ("Front Plane", "Top Plane", "Right Plane")]
        self.selected = None
        self.sel_list = []
        self.ext_sel = []                          # (kind, object, mark): sketches, planes, closed-sketch entities
        self.dims = []
        self.rels = []
        self.body = None
        self.Extension = Extension(self)
        self.SketchManager = SketchManager(self)
        self.FeatureManager = FeatureManager(self)
        self.SelectionManager = SelectionManager(self)
        LAST["model"] = self

    def _ddim(self):
        """Display dimension of the last dimension added: name D<n>@Sketch<k>, value in m."""
        k = sum(f.kind == "ProfileFeature" for f in self.features) + 1
        self._dn = getattr(self, "_dn", 0) + 1
        kind, v = self.dims[-1]
        val = v * 1e-3 if kind != "angle" else v

        class Dim:
            FullName = f"D{self._dn}@Sketch{k}@Part1.Part"
            SystemValue = val

        class DDim:
            def GetDimension2(self, i):
                return Dim()
        return DDim()

    def Parameter(self, name):
        class D:
            SystemValue = 0.0
        return D()

    def EditUndo2(self, steps):
        return True

    def GetEquationMgr(self):
        if not hasattr(self, "_eq"):
            class EqMgr:
                def __init__(self):
                    self.eqs = []

                def Add2(self, i, text, solve):
                    self.eqs.append(text)
                    return len(self.eqs) - 1

                def GetCount(self):
                    return len(self.eqs)

                def Delete(self, i):
                    del self.eqs[i]
                    return 1
            self._eq = EqMgr()
        return self._eq

    def EditRebuild3(self):
        return True

    def GetActiveSketch2(self):
        return None if HIDE_SKETCH else self.SketchManager._sk

    def GetBodies2(self, typ, visible):
        return [Body(self)] if self.body is not None else None

    def FeatureByName(self, name):
        for f in self.features:
            if f.kind == "RefPlane" and f.data == name:
                return f
        for f in self.features:
            if f.kind != "RefPlane" and f.Name == name:
                return f
        return None

    def BlankRefGeom(self):
        pass

    def BlankSketch(self):
        pass

    def ClearSelection2(self, flag):
        self.selected = None
        self.sel_list = []
        self.ext_sel = []
        self.ext_segs = []

    def AddDimension2(self, x, y, z):
        s = self.sel_list
        planes = [o for k, o, _ in self.ext_sel if k == "feature" and o.kind == "RefPlane"]
        if len(s) == 1 and s[0][0] == "point" and len(s[0][1]) == 3 and len(planes) == 1:
            pl = planes[0]
            fr, o = (pl.data["frame"], pl.data["origin"]) if isinstance(pl.data, dict) else (FRAMES[pl.data], np.zeros(3))
            d = abs(float((np.array(s[0][1]) - np.asarray(o, float)) @ np.array(fr[2], float)))
            self.dims.append(("distance", d / 1e-3))
            return self._ddim()
        if self._two_points() is not None:
            p, q = self._two_points()
            self.dims.append(("distance", float(np.linalg.norm(p - q)) / 1e-3))
            return self._ddim()
        if len(s) == 1 and s[0][0] == "circle":
            self.dims.append(("diameter", 2 * s[0][1].edge[2] / 1e-3))
            return self._ddim()
        if len(s) == 1 and s[0][0] == "line":
            e = s[0][1].edge
            self.dims.append(("length", float(np.linalg.norm(np.array(e[2]) - np.array(e[1]))) / 1e-3))
            return self._ddim()
        if len(s) == 1 and s[0][0] == "arc":
            e = s[0][1].edge
            p, q, r = (np.array(v) for v in e[1:])
            a = 2 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1]))
            ux = ((p @ p) * (q[1] - r[1]) + (q @ q) * (r[1] - p[1]) + (r @ r) * (p[1] - q[1])) / a
            uy = ((p @ p) * (r[0] - q[0]) + (q @ q) * (p[0] - r[0]) + (r @ r) * (q[0] - p[0])) / a
            self.dims.append(("radius", float(np.hypot(p[0] - ux, p[1] - uy)) / 1e-3))
            return self._ddim()
        return None

    def SketchAddConstraints(self, rel):
        s = self.sel_list
        if rel == "sgATPIERCE":                       # a point and a segment of another (closed) sketch
            assert [k for k, _ in s] == ["point"] and any(k == "segment" for k, _, _ in self.ext_sel), (rel, s)
            self.rels.append(rel)
            return
        ext_pts = [o for k, o, _ in self.ext_sel if k == "point"]
        if rel == "sgCOINCIDENT" and [k for k, _ in s] == ["point"] and len(ext_pts) == 1:
            sk = self.SketchManager._sk
            q = s[0][1]
            p = sk.to_model(q) if len(q) == 2 else np.asarray(q, float)
            assert np.linalg.norm(p - ext_pts[0]) < 1e-9, ("external point not at the sketch point", p, ext_pts[0])
            self.rels.append(rel)
            return
        if rel == "sgCOINCIDENT" and sorted(k for k, _ in s) == ["arc", "point"]:     # a point on an arc's circle
            self.rels.append(rel)
            return
        if rel in ("sgSYMMETRIC", "sgATMIDDLE") or (rel == "sgCOINCIDENT" and sorted(k for k, _ in s) == ["line", "point"]):
            need = {"sgSYMMETRIC": ["line", "point", "point"], "sgATMIDDLE": ["line", "point"],
                    "sgCOINCIDENT": ["line", "point"]}[rel]
            assert sorted(k for k, _ in s) == need, (rel, s)
            self.rels.append(rel)
            return
        if rel == "sgPARALLEL":
            assert sorted(k for k, _ in s) == ["line", "line"], (rel, s)
            self.rels.append(rel)
            return
        if rel == "sgCONCENTRIC":                    # a circle and an arc sharing their centre
            ks = sorted(k for k, _ in s)
            assert ks in (["arc", "circle"], ["circle", "circle"], ["arc", "arc"]), (rel, s)
            self.rels.append(rel)
            return
        if rel == "sgSAMELENGTH":                    # Equal: two lines or two arcs/circles
            ks = sorted(k for k, _ in s)
            assert ks in (["line", "line"], ["arc", "arc"], ["arc", "circle"], ["circle", "circle"]), (rel, s)
            self.rels.append(rel)
            return
        if rel in ("sgHORIZONTAL2D", "sgVERTICAL2D"):      # as SolidWorks: the line would turn
            (k0, seg), = s
            (x1, y1), (x2, y2) = seg.edge[1], seg.edge[2]
            along, across = (abs(x2 - x1), abs(y2 - y1)) if rel == "sgHORIZONTAL2D" else (abs(y2 - y1), abs(x2 - x1))
            assert across <= 1e-3 * along, f"{rel} on a line that is not {'horizontal' if rel == 'sgHORIZONTAL2D' else 'vertical'}"
        need = {"sgHORIZONTAL2D": ["line"], "sgVERTICAL2D": ["line"], "sgTANGENT": ["arc", "line"],
                "sgCOINCIDENT": ["point", "point"], "sgHORIZONTALPOINTS2D": ["point", "point"],
                "sgVERTICALPOINTS2D": ["point", "point"]}[rel]
        assert sorted(k for k, _ in s) == sorted(need), (rel, s)
        self.rels.append(rel)

    def _two_points(self):
        s = self.sel_list
        if len(s) == 2 and all(k == "point" for k, _ in s):
            return np.array(s[0][1]), np.array(s[1][1])
        return None

    def AddHorizontalDimension2(self, x, y, z):
        p = self._two_points()
        if p is None:
            return None
        self.dims.append(("horizontal", abs(p[0][0] - p[1][0]) / 1e-3))
        return self._ddim()

    def AddVerticalDimension2(self, x, y, z):
        p = self._two_points()
        if p is None:
            return None
        self.dims.append(("vertical", abs(p[0][1] - p[1][1]) / 1e-3))
        return self._ddim()

    def EditDelete(self):
        if os.environ.get("SWMOCK_NO_DELETE") == "1":
            return
        if isinstance(self.selected, Segment):
            self.selected.deleted = True
            self.selected = None

    def FirstFeature(self):
        return self.features[0]

    def FeatureByPositionReverse(self, i):
        return self.features[-1 - i]

    def ForceRebuild3(self, flag):
        return True

    def SaveAs3(self, path, a, b):
        import json
        with open(path + ".dims.json", "w") as fh:
            json.dump(dict(dims=self.dims, rels=self.rels), fh)
        self.body.exportStep(path)                                        # mm
        self.body.exportBrep(path + ".brep")     # exact copy: STEP lost a node of revolution (swsweep)
        return 0


class App:
    Visible = False
    prefs = {}

    def GetUserPreferenceToggle(self, k):
        return self.prefs.get(k, True)

    def SetUserPreferenceToggle(self, k, v):
        self.prefs[k] = v

    def CloseAllDocuments(self, include_unsaved):
        return True

    def NewPart(self):
        return Model()

    def GetMathUtility(self):
        return MathUtility()


def Dispatch(name):
    assert name == "SldWorks.Application"
    return App()
