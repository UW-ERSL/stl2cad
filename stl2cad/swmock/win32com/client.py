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
    def __init__(self, frame, offset=0.0):
        u, v, n = (np.array(a, float) for a in frame)
        self.M = np.vstack([u, v, n])
        R = np.column_stack([u, v, n])
        self.ArrayData = tuple(R.ravel()) + (0.0, 0.0, -offset, 1.0, 0.0, 0.0, 0.0)

    def __matmul__(self, p):
        return self.M @ p


class Segment:
    def __init__(self, sm, edge):
        self.sm, self.edge = sm, edge
        self.ConstructionGeometry = False
        self.deleted = False

    def Select4(self, append, data):
        m = self.sm.model
        m.selected = self
        if not append:
            m.sel_list = []
        m.sel_list.append((self.edge[0], self))
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
    def __init__(self, plane, frame=None, offset=0.0):
        self.plane = plane
        self.frame = FRAMES[plane] if frame is None else frame
        self.offset = offset                       # m, plane position along its normal
        self.ModelToSketchTransform = Transform(self.frame, offset)
        self.edges = []


class Feature:
    def __init__(self, model, kind, data=None):
        self.model, self.kind, self.data, self.Name = model, kind, data, kind

    def Select2(self, append, mark):
        self.model.selected = self
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
            else:
                self._sk = Sketch("Front Plane" if FRONT_ONLY else m.selected.data)
            self._sk.open = False
            self._lazy = LAZY
            m.selected = None
        else:
            m.features.append(Feature(m, "ProfileFeature", self._sk))
            self._sk = None

    def CreateLine(self, x1, y1, z1, x2, y2, z2):
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
        self._sk.edges.append(seg)
        return seg


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

    def SelectByID2(self, name, typ, x, y, z, append, mark, callout, opt):
        m = self.model
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
            P = cq.Vector(x * K, y * K, z * K)
            for fc in m.body.Faces():
                if fc.geomType() != "PLANE":
                    continue
                if (cq.Vertex.makeVertex(*P.toTuple())).distance(fc) < 1e-6:
                    nrm = np.array(fc.normalAt(P).toTuple())
                    fr = face_frame(nrm)
                    off = float(nrm @ np.array(fc.Center().toTuple())) / K     # the face's own plane
                    m.selected = Feature(m, "Face", (fr, off))
                    return True
            return False
        if typ == "EXTSKETCHPOINT" and name == "Point1@Origin":
            m.sel_list.append(("point", (0.0, 0.0)))
            return True
        for f in self.model.features:
            if f.kind == "RefPlane" and f.data == name:
                self.model.selected = f
                return True
        return False

    def CreateMassProperty(self):
        b = self.model.body
        if b is None:
            return MassProperty(0.0, (0.0, 0.0, 0.0))
        c = b.Center()
        return MassProperty(b.Volume() / K ** 3, (c.x / K, c.y / K, c.z / K))   # m^3, m

    def DeleteSelection2(self, opt):
        f = self.model.selected
        assert f is self.model.features[-1] and f.kind in ("Extrusion", "ProfileFeature", "MoveCopyBody")
        self.model.features.pop()
        if f.kind in ("Extrusion", "MoveCopyBody"):
            self.model.body = f.data["before"]
        return True


def _boolean(a, b, cut):
    """a - b or a + b. OCC booleans with near-coincident faces can return an invalid solid that
    later booleans get wrong (seen: a boss whose circle is a flange's end arc, then a cut beside
    the body removing material); a result that is invalid or breaks the volume identity
    vol(a) +- vol(b) -+ vol(a & b) is recomputed with a fuzzy tolerance."""
    def run(tol, clean=True, glue=False):
        r = a.cut(b, tol=tol) if cut else a.fuse(b, tol=tol, glue=glue)
        return r.clean() if clean else r
    common = a.intersect(b).Volume()
    want = a.Volume() - common if cut else a.Volume() + b.Volume() - common
    good = lambda r: r.isValid() and abs(r.Volume() - want) <= 1e-6 * max(1.0, abs(want))
    r = run(None)
    if good(r):
        return r
    for tol, cl, gl in ((None, False, False), (None, True, True), (1e-5, True, False), (1e-5, False, False),
                        (1e-4, True, False), (1e-3, True, False)):
        r2 = run(tol, cl, gl and not cut)  # near-coincident faces (a boss on a flange's end arc)
        if good(r2):
            return r2
        r2 = r2.fix()
        if good(r2):
            return r2
    return r


class FeatureManager:
    def __init__(self, model):
        self.model = model

    def _solid(self, sk, rev, h, t0, off, flip, cut):
        u, v, n = (np.array(a, float) for a in sk.frame)
        o = sk.offset * n
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
        for w in wires:
            f = cq.Face.makeFromWires(w, [])
            # always extruded along +n from the lower end: the same solid, but OCC fuses
            # coincident faces (a boss on a flange's end arc) reliably only when they are
            # parametrized alike
            s = cq.Solid.extrudeLinear(f, cq.Vector(*(n * h)))
            solids.append(s.translate(cq.Vector(*(n * (start if d > 0 else start - h)))))
        out = solids[0]
        for s in solids[1:]:
            out = out.fuse(s)
        return out

    def _apply(self, rev, h, t0, off, flip, cut):
        m = self.model
        sk = m.selected.data
        if getattr(sk, "open", False):
            return None
        body = self._solid(sk, rev, h, t0, off, flip, cut)
        before = m.body
        if cut:
            if before is None:
                return None
            m.body = _boolean(before, body, cut=True)
        else:
            m.body = body if before is None else _boolean(before, body, cut=False)
        f = Feature(m, "Extrusion", dict(before=before))
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
                b = b.rotate(cq.Vector(0, 0, 0), cq.Vector(*axis), sgn * math.degrees(ang))
        m.body = b
        f = Feature(m, "MoveCopyBody", dict(before=before))
        m.features.append(f)
        return f

    def FeatureCut4(self, Sd, Flip, Dir, T1, T2, D1, D2, *rest):
        t0, off, flip = rest[-4], rest[-3], rest[-2]
        return self._apply(Dir, D1, t0, off, flip, True)


class SelectionManager:
    def __init__(self, model):
        self.model = model

    def CreateSelectData(self):
        return SelectData()

    def GetSelectedObjectCount2(self, mark):
        return 0 if self.model.selected is None else 1


class Model:
    def __init__(self):
        self.features = [Feature(self, "RefPlane", p) for p in ("Front Plane", "Top Plane", "Right Plane")]
        self.selected = None
        self.sel_list = []
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
        return None

    def ClearSelection2(self, flag):
        self.selected = None
        self.sel_list = []

    def AddDimension2(self, x, y, z):
        s = self.sel_list
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
        if rel in ("sgSYMMETRIC", "sgATMIDDLE") or (rel == "sgCOINCIDENT" and sorted(k for k, _ in s) == ["line", "point"]):
            need = {"sgSYMMETRIC": ["line", "point", "point"], "sgATMIDDLE": ["line", "point"],
                    "sgCOINCIDENT": ["line", "point"]}[rel]
            assert sorted(k for k, _ in s) == need, (rel, s)
            self.rels.append(rel)
            return
        if rel == "sgSAMELENGTH":                    # Equal: two lines or two arcs/circles
            ks = sorted(k for k, _ in s)
            assert ks in (["line", "line"], ["arc", "arc"], ["arc", "circle"], ["circle", "circle"]), (rel, s)
            self.rels.append(rel)
            return
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
        return 0


class App:
    Visible = False
    prefs = {}

    def GetUserPreferenceToggle(self, k):
        return self.prefs.get(k, True)

    def SetUserPreferenceToggle(self, k, v):
        self.prefs[k] = v

    def NewPart(self):
        return Model()

    def GetMathUtility(self):
        return MathUtility()


def Dispatch(name):
    assert name == "SldWorks.Application"
    return App()
