"""Mesh loading, unit handling, deviation metrics, and shared helpers."""

from dataclasses import dataclass, field

import numpy as np
import trimesh

UNITS_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0, "in": 25.4}


def load_mesh(path, units="mm"):
    """Load an STL and return a cleaned trimesh in millimetres."""
    mesh = trimesh.load(path, force="mesh")
    mesh.merge_vertices()
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.update_faces(mesh.unique_faces())
    mesh.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(mesh)
    mesh.apply_scale(UNITS_TO_MM[units])
    return mesh


def merge_overlapping(mesh, log=print):
    """An STL exported as several overlapping closed shells (a hub, a web and a toothed ring
    written as separate bodies) describes their union. Return the union when shells overlap,
    otherwise the mesh unchanged. Needs the manifold3d package; without it the mesh is returned."""
    if mesh.body_count < 2:
        return mesh
    bodies = mesh.split(only_watertight=False)
    if not all(b.is_watertight for b in bodies):
        return mesh
    try:
        u = trimesh.boolean.union(bodies, engine="manifold")
    except Exception as e:                       # manifold3d missing or failed
        log(f"  {len(bodies)} closed shells; union not computed ({type(e).__name__})")
        return mesh
    vsum = sum(b.volume for b in bodies)
    if not u.is_watertight or vsum - u.volume < 1e-6 * vsum:
        return mesh
    log(f"  {len(bodies)} overlapping closed shells merged into {u.body_count} body/bodies "
        f"(overlap {vsum - u.volume:.1f} mm^3)")
    return u


def default_tolerance(mesh):
    """Fit tolerance: 0.05 % of the bounding-box diagonal (suits CAD-exported STL)."""
    return 5e-4 * float(np.linalg.norm(mesh.extents))


def closest_distance(mesh, pts, chunk=500):
    """Distance from each point to the mesh surface, queried in chunks: one trimesh query over
    many points against long thin facets built candidate arrays of several GB (WedgeHoles)."""
    pts = np.asarray(pts, float)
    out = np.empty(len(pts))
    for i in range(0, len(pts), chunk):
        _, out[i:i + chunk], _ = trimesh.proximity.closest_point(mesh, pts[i:i + chunk])
    return out


def deviation(reference, model, n=30000, seed=0):
    """Two-sided surface deviation statistics between two meshes (mm)."""
    out = {}
    for name, a, b in (("ref_to_model", reference, model), ("model_to_ref", model, reference)):
        pts, _ = trimesh.sample.sample_surface_even(a, n, seed=seed)
        d = closest_distance(b, pts)
        out[name] = dict(mean=float(d.mean()), median=float(np.median(d)),
                         p95=float(np.percentile(d, 95)), max=float(d.max()))
    out["volume_ref"] = float(reference.volume)
    out["volume_model"] = float(model.volume)
    return out


def vertex_deviation(reference, model):
    """Distance from each vertex of `model` to `reference`, for color maps."""
    return closest_distance(reference, model.vertices)


def heal_via_step(shape):
    """Write the shape to STEP and read it back: the STEP reader's shape healing gives clean
    topology and tolerances (the in-memory boolean result can mesh with gaps)."""
    import os
    import tempfile
    import cadquery as cq
    fd, path = tempfile.mkstemp(suffix=".step")
    os.close(fd)
    try:
        shape.exportStep(path)
        healed = cq.importers.importStep(path).val()
    except Exception:
        return shape
    finally:
        os.remove(path)
    if healed.Solids() and abs(healed.Volume() - shape.Volume()) < 1e-4 * abs(shape.Volume()) and (
            healed.isValid() or not shape.isValid()):
        return healed.Solids()[0] if len(healed.Solids()) == 1 else healed
    return shape


def shape_to_mesh(shape, tolerance=0.05, angular=0.2):
    """CadQuery shape -> watertight trimesh (OCC incremental mesher, as used for STL export)."""
    import os
    import tempfile
    from OCP.BRepTools import BRepTools
    BRepTools.Clean_s(shape.wrapped, True)    # drop stale triangulations left on shared faces
    fd, path = tempfile.mkstemp(suffix=".stl")
    os.close(fd)
    try:
        shape.exportStl(path, tolerance=tolerance, angularTolerance=angular)
        mesh = trimesh.load(path, force="mesh")
    finally:
        os.remove(path)
    mesh.merge_vertices()
    return mesh


def format_deviation(dev):
    lines = []
    for key in ("ref_to_model", "model_to_ref"):
        d = dev[key]
        lines.append(f"  {key:13s} mean {d['mean']:.4f}  median {d['median']:.4f}  "
                     f"p95 {d['p95']:.4f}  max {d['max']:.4f} mm")
    vr, vm = dev["volume_ref"], dev["volume_model"]
    lines.append(f"  volume ref {vr:.2f}  model {vm:.2f} mm^3  ({100 * (vm / vr - 1):+.2f} %)")
    return "\n".join(lines)


class Log:
    """Minimal logger that forwards to a callback (GUI) or stdout (CLI)."""

    def __init__(self, sink=None):
        self.sink = sink or (lambda s: print(s, flush=True))

    def __call__(self, msg):
        self.sink(str(msg))


@dataclass
class Result:
    shape: object = None            # CadQuery shape
    mesh: object = None             # tessellated result (trimesh)
    mode: str = ""
    report: str = ""
    deviation: dict = field(default_factory=dict)
    extras: dict = field(default_factory=dict)
