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
    def __init__(self, frame):
        u, v, n = (np.array(a, float) for a in frame)
        self.M = np.vstack([u, v, n])

    def __matmul__(self, p):
        return self.M @ p


class Segment:
    def __init__(self, sm, edge):
        self.sm, self.edge = sm, edge
        self.ConstructionGeometry = False
        self.deleted = False

    def Select4(self, append, data):
        self.sm.model.selected = self
        return True


class Sketch:
    def __init__(self, plane):
        self.plane = plane
        self.ModelToSketchTransform = Transform(FRAMES[plane])
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
            assert m.selected is not None and m.selected.kind == "RefPlane", "no plane selected"
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
                if seg.deleted or seg.ConstructionGeometry:
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


class FeatureManager:
    def __init__(self, model):
        self.model = model

    def _solid(self, sk, rev, h, t0, off, flip, cut):
        u, v, n = (np.array(a, float) for a in FRAMES[sk.plane])
        P = lambda q: cq.Vector(*(K * (q[0] * u + q[1] * v)))      # geometry kept in mm
        es = []
        for seg in sk.edges:
            if seg.deleted or seg.ConstructionGeometry:
                continue
            e = seg.edge
            if e[0] == "line":
                es.append(cq.Edge.makeLine(P(e[1]), P(e[2])))
            elif e[0] == "arc":
                es.append(cq.Edge.makeThreePointArc(P(e[1]), P(e[2]), P(e[3])))
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
            s = cq.Solid.extrudeLinear(f, cq.Vector(*(n * h * d)))
            solids.append(s.translate(cq.Vector(*(n * start))))
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
            m.body = before.cut(body)
        else:
            m.body = body if before is None else before.fuse(body).clean()
        f = Feature(m, "Extrusion", dict(before=before))
        m.features.append(f)
        LAST["features"] = m.features
        return f

    def FeatureExtrusion3(self, Sd, Flip, Dir, T1, T2, D1, D2, *rest):
        t0, off, flip = rest[-3], rest[-2], rest[-1]
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
        if len(s) == 1 and s[0][0] == "circle":
            self.dims.append(("diameter", 2 * s[0][1].edge[2] / 1e-3))
            return object()
        if len(s) == 1 and s[0][0] == "line":
            e = s[0][1].edge
            self.dims.append(("length", float(np.linalg.norm(np.array(e[2]) - np.array(e[1]))) / 1e-3))
            return object()
        if len(s) == 1 and s[0][0] == "arc":
            e = s[0][1].edge
            p, q, r = (np.array(v) for v in e[1:])
            a = 2 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1]))
            ux = ((p @ p) * (q[1] - r[1]) + (q @ q) * (r[1] - p[1]) + (r @ r) * (p[1] - q[1])) / a
            uy = ((p @ p) * (r[0] - q[0]) + (q @ q) * (p[0] - r[0]) + (r @ r) * (q[0] - p[0])) / a
            self.dims.append(("radius", float(np.hypot(p[0] - ux, p[1] - uy)) / 1e-3))
            return object()
        return None

    def SketchAddConstraints(self, rel):
        s = self.sel_list
        need = {"sgHORIZONTAL2D": ["line"], "sgVERTICAL2D": ["line"], "sgTANGENT": ["arc", "line"],
                "sgCOINCIDENT": ["point", "point"]}[rel]
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
        return object()

    def AddVerticalDimension2(self, x, y, z):
        p = self._two_points()
        if p is None:
            return None
        self.dims.append(("vertical", abs(p[0][1] - p[1][1]) / 1e-3))
        return object()

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
