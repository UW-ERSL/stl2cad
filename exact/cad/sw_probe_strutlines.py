"""Probe: for each strut of truss-connector_solidworks.py in the open SolidWorks part, the 3D sketch line that
SolidWorks holds (its segments and end points) and the cylinder faces of the sweep, against the data.
Run with the truss-connector part open (after truss-connector_solidworks.py):
    python sw_probe_strutlines.py            (all struts)
    python sw_probe_strutlines.py 20 27      (these struts)"""
import json
import sys

import numpy as np
import win32com.client

src = open("truss-connector_solidworks.py").read()
DATA = json.loads(src.split("DATA_JSON = r'''")[1].split("'''")[0])


def call0(obj, name):
    """A property or a method without arguments (late binding gives either); a COM object a property
    returns is callable too, and is not called."""
    a = getattr(obj, name)
    if isinstance(a, win32com.client.CDispatch) or not callable(a):
        return a
    return a()


def pt(p):
    return np.array([call0(p, "X"), call0(p, "Y"), call0(p, "Z")]) * 1000.0


app = win32com.client.Dispatch("SldWorks.Application")
model = app.ActiveDoc
ks = [int(a) for a in sys.argv[1:]] or list(range(1, len(DATA["struts"]) + 1))
bad = 0
for k in ks:
    s = DATA["struts"][k - 1]
    p0, p1 = np.array(s["p0"]), np.array(s["p1"])
    u = (p1 - p0) / np.linalg.norm(p1 - p0)
    L = np.linalg.norm(p1 - p0)
    lines = []
    fs = model.FeatureByName(f"Strut {k} line")
    if fs is None:
        print(f"strut {k}: no sketch 'Strut {k} line'")
        continue
    sk = call0(fs, "GetSpecificFeature2")
    segs = call0(sk, "GetSketchSegments") or ()
    for g in segs:
        try:
            a, b = pt(call0(g, "GetStartPoint2")), pt(call0(g, "GetEndPoint2"))
        except Exception:
            lines.append(("non-line", None, None))
            continue
        lines.append(("line", a, b))
    msg = []
    if len(lines) != 1:
        msg.append(f"{len(lines)} segments in the sketch")
    for kind, a, b in lines:
        if a is None:
            continue
        d = min(np.linalg.norm(a - p0) + np.linalg.norm(b - p1), np.linalg.norm(a - p1) + np.linalg.norm(b - p0))
        if d > 1e-3:
            ta, tb = (a - p0) @ u / L, (b - p0) @ u / L
            msg.append(f"line {np.round(a, 3)} -> {np.round(b, 3)} (t {ta:.3f} .. {tb:.3f}, off by {d:.3f} mm)")
    ff = model.FeatureByName(f"Strut {k}")
    if ff is None:
        msg.append("no feature")
    else:
        for f in call0(ff, "GetFaces") or ():
            srf = call0(f, "GetSurface")
            if not call0(srf, "IsCylinder"):
                continue
            prm = np.array(call0(srf, "CylinderParams"))
            box = np.array(call0(f, "GetBox")) * 1000.0
            C = np.array([[x, y, z] for x in box[[0, 3]] for y in box[[1, 4]] for z in box[[2, 5]]])
            t = (C - p0) @ u / L
            r = prm[6] * 1000.0
            if abs(r - s["r"]) > 1e-3:
                continue
            lo, hi = t.min(), t.max()
            pad = 1.5 * s["r"] / L
            if lo < -pad or hi > 1 + pad:
                msg.append(f"cylinder face box spans t {lo:.3f} .. {hi:.3f} (strut 0 .. 1)")
    if msg:
        bad += 1
        print(f"strut {k} p0 {np.round(p0, 3)} p1 {np.round(p1, 3)} r {s['r']}: " + "; ".join(msg))
    elif len(ks) < 10:
        print(f"strut {k}: sketch line and sweep match the data")
print(f"{bad} of {len(ks)} strut(s) differ (face boxes: a merged face can include a neighbour's)")

# exact check: tessellated faces of the part; points on a strut's cylinder (distance r from its axis line)
# beyond its ends, with the feature that owns the face
if len(sys.argv) > 1:
    bodies = model.GetBodies2(0, True)
    faces = []
    for b in bodies or ():
        faces += list(call0(b, "GetFaces") or ())
    tess = []
    for f in faces:
        tri = f.GetTessTriangles(True)
        X = np.array(tri, float).reshape(-1, 3) * 1000.0 if tri else np.zeros((0, 3))
        try:
            owner = call0(call0(f, "GetFeature"), "Name")
        except Exception:
            owner = "?"
        tess.append((owner, X))
    for k in ks:
        s = DATA["struts"][k - 1]
        p0, p1 = np.array(s["p0"]), np.array(s["p1"])
        L = np.linalg.norm(p1 - p0)
        u = (p1 - p0) / L
        for owner, X in tess:
            if not len(X):
                continue
            t = (X - p0) @ u / L
            rad = np.linalg.norm((X - p0) - np.outer(t * L, u), axis=1)
            on = (np.abs(rad - s["r"]) < 0.02 * s["r"]) & ((t < -0.02) | (t > 1.02))
            if on.sum() >= 3:
                print(f"strut {k}: face of '{owner}' on its cylinder beyond the ends: t {t[on].min():.3f} .. "
                      f"{t[on].max():.3f}, z {X[on, 2].min():.2f} .. {X[on, 2].max():.2f} ({on.sum()} points)")
