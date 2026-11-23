"""High-level pipeline shared by the GUI and the command line."""

import json
import os
import time
from dataclasses import asdict, dataclass, field

import numpy as np

from .core import Log, default_tolerance, deviation, format_deviation, load_mesh, shape_to_mesh
from .symmetry import detect, group_matrices, spec_from_options


@dataclass
class Options:
    units: str = "mm"                  # mm, cm, m, in
    mode: str = "auto"                 # auto, exact, smooth
    tolerance: float = 0.0             # mm; 0 = automatic
    # symmetry
    symmetry: str = "auto"             # off, auto, manual
    mirror_x: bool = False
    mirror_y: bool = False
    mirror_z: bool = False
    cyclic: int = 1
    cyclic_axis: str = "z"             # frame axis of the rotation
    dihedral: bool = False
    # exact mode
    surfaces: list = field(default_factory=lambda: ["plane", "cylinder", "cone", "sphere", "torus"])
    algebraic_only: bool = False       # no spline extrusions (exact) / cones and spheres only (smooth)
    allow_extrusion: bool = True
    snap_directions_deg: float = 1.0
    snap_angle_step_deg: float = 15.0
    tangency: bool = True
    separators: bool = False           # experimental: strips at near-tangent junctions
    fillets: bool = True               # rebuild tangent blends as fillets on a sharp solid
    prism25: bool = True               # try stacked prisms (2.5D) before the surface arrangement
    export_script: bool = True         # write a verified CadQuery script (exact mode)
    export_solidworks: bool = True     # write a SolidWorks feature-tree script (2.5D results)
    grid: float = 0.0                  # mm; 0 = off
    # smooth mode
    pitch: float = 0.0                 # voxel size, mm; 0 = automatic
    blur: float = 0.0                  # mm; 0 = 3 x pitch
    prune: float = 1.0                 # spur pruning, multiples of local radius
    min_copy_frac: float = 0.5
    holes: bool = True
    keep_free_ends: bool = True
    clip: str = "thin"                 # none, thin, all
    area_radius: bool = True
    # output
    out_dir: str = ""
    export_step: bool = True
    export_stl: bool = True

    def to_json(self):
        return json.dumps(asdict(self), indent=2)


def coplanar_fraction(mesh):
    """Share of surface area on triangles with an exactly coplanar neighbour (CAD-like if high)."""
    ang = mesh.face_adjacency_angles
    adj = mesh.face_adjacency
    flat = np.zeros(len(mesh.faces), bool)
    f = adj[ang < np.radians(0.05)]
    flat[f[:, 0]] = True
    flat[f[:, 1]] = True
    return float(mesh.area_faces[flat].sum() / mesh.area)


def resolve_symmetry(mesh, o, tol, log):
    if o.symmetry == "off":
        return None, "symmetry: off"
    if o.symmetry == "auto":
        spec, rep = detect(mesh, tol, log=lambda s: None)
        log("symmetry (detected): " + rep.splitlines()[-2])
        return (spec if len(group_matrices(spec)) > 1 else None), rep
    mirrors = tuple(i for i, f in enumerate((o.mirror_x, o.mirror_y, o.mirror_z)) if f)
    axis = "xyz".index(o.cyclic_axis)
    spec = spec_from_options(mesh, mirrors=mirrors, cyclic=o.cyclic, axis=axis, dihedral=o.dihedral, log=log)
    return (spec if len(group_matrices(spec)) > 1 else None), spec.describe()


def run(path, o: Options, log=None, mesh=None):
    log = Log(log)
    t0 = time.time()
    if mesh is None:
        mesh = load_mesh(path, o.units)
    log(f"loaded {os.path.basename(path)}: {len(mesh.faces)} triangles, extents "
        f"{np.round(mesh.extents, 3).tolist()} mm, watertight={mesh.is_watertight}")
    from .core import merge_overlapping
    mesh = merge_overlapping(mesh, log)
    tol = o.tolerance or default_tolerance(mesh)
    spec, sym_report = resolve_symmetry(mesh, o, tol, log)
    mode = o.mode
    cf_ = coplanar_fraction(mesh)
    if mode == "auto":
        mode = "exact" if cf_ > 0.3 else "smooth"
        log(f"auto mode: coplanar-area fraction {cf_:.2f} -> {mode}")

    res = None
    if mode == "exact":
        from . import exact
        kinds = [k for k in o.surfaces]
        opts = dict(tol=tol, kinds=kinds, allow_extrusion=o.allow_extrusion and not o.algebraic_only,
                    allow_facets=True, snap_deg=o.snap_directions_deg, angle_step=o.snap_angle_step_deg,
                    tangency=o.tangency, grid=o.grid, impose_symmetry=True, separators=o.separators,
                    fillets=o.fillets, prism25=o.prism25)
        # exact mode reproduces the part; a detected (not imposed) symmetry is only reported
        exact_spec = spec if o.symmetry == "manual" else None
        if spec is not None and exact_spec is None:
            log("exact mode: detected symmetry is reported, not imposed (choose 'Impose' to enforce it)")
        try:
            opts["symmetry_hint"] = spec          # detected group: used only to propose fits
            res = exact.reconstruct(mesh, opts, spec=exact_spec, log=log)
            if not res.extras.get("passed", False):
                log("WARNING: exact result failed the volume/validity check; inspect the deviation map")
                if o.mode == "auto" and cf_ < 0.6:
                    log("falling back to smooth mode")
                    res = None
                    mode = "smooth"
        except Exception as e:
            log(f"exact mode failed: {e}")
            if o.mode != "auto" or cf_ >= 0.6:
                raise
            mode = "smooth"
    if res is None and mode == "smooth":
        from . import smooth
        if spec is not None and o.symmetry == "auto" and spec.cyclic > 1 and (spec.dihedral or spec.mirrors):
            from .symmetry import SymmetrySpec
            full = spec.describe()
            spec = SymmetrySpec(center=spec.center, axes=spec.axes, cyclic=spec.cyclic, cyclic_axis=spec.cyclic_axis)
            log(f"smooth mode: detected {full}; imposing the rotational part {spec.describe()} "
                f"(mirror planes are reported, not imposed)")
        opts = dict(pitch=o.pitch, blur=o.blur, prune=o.prune, min_copy_frac=o.min_copy_frac, holes=o.holes,
                    keep_free_ends=o.keep_free_ends, clip=o.clip, area_radius=o.area_radius,
                    algebraic_only=o.algebraic_only)
        res = smooth.reconstruct(mesh, opts, spec=spec, log=log)

    diag = float(np.linalg.norm(mesh.extents))
    from .core import heal_via_step
    res.shape = heal_via_step(res.shape)
    from .brep_helpers import clean_mesh
    clean_mesh(res.shape)
    res.valid = res.shape.isValid()          # checked before meshing (the mesh leaves edge polygons)
    res.mesh = shape_to_mesh(res.shape, tolerance=max(2e-4 * diag, 0.25 * tol), angular=0.25)
    res.deviation = deviation(mesh, res.mesh)
    res.extras["input_mesh"] = mesh
    res.extras["symmetry"] = sym_report
    res.extras["tolerance"] = tol
    used = spec if (res.mode == "smooth" or o.symmetry == "manual") else None
    lines = [f"mode: {res.mode}", f"tolerance: {tol:.4f} mm",
             f"symmetry imposed: {used.describe() if used else 'none'}"
             + (f" (detected: {spec.describe()})" if spec is not None and used is None else ""),
             f"solid valid: {res.valid}, faces: {len(res.shape.Faces())}",
             "deviation:", format_deviation(res.deviation), f"time: {time.time() - t0:.1f} s"]
    res.report = "\n".join(lines)
    log(res.report)
    if o.out_dir:
        save(res, path, o, log)
    return res


def save(res, path, o, log=print, out_dir=None):
    out_dir = out_dir or o.out_dir
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(path))[0] + f"_{res.mode}"
    files = []
    if o.export_step:
        f = os.path.join(out_dir, stem + ".step")
        res.shape.exportStep(f)
        files.append(f)
    if o.export_stl:
        f = os.path.join(out_dir, stem + ".stl")
        res.mesh.export(f)
        files.append(f)
    if o.export_script and res.mode == "exact":
        from .codegen import generate
        t0 = time.time()
        code, kind, msg = generate(res, stem, res.extras.get("tolerance", 0.1), log=log)
        if code is not None:
            f = os.path.join(out_dir, stem + ".py")
            with open(f, "w") as fh:
                fh.write(code)
            files.append(f)
            log(f"script: {kind} script, {msg} ({time.time() - t0:.1f}s)")
        else:
            log(f"script: not written ({msg})")
        res.report += f"\nscript: {kind + ' script, ' + msg if code else 'not written (' + msg + ')'}"
    if o.export_solidworks and res.mode == "exact" and (res.extras.get("prism25") is not None
                                                         or res.extras.get("revolve") is not None):
        from .swgen import solidworks_script, solidworks_revolve_script, verify_with_standin
        t0 = time.time()
        swname = os.path.splitext(os.path.basename(path))[0] + "_solidworks"
        tol_ = res.extras.get("tolerance", 0.1)
        try:
            if res.extras.get("revolve") is not None:
                code = solidworks_revolve_script(res.extras["revolve"], res.shape, swname)
            else:
                code = solidworks_script(res.extras["prism25"], swname)
        except ValueError as e:                  # a feature the generator cannot express
            code, ok1, ok2, m1 = None, False, False, str(e)
        if code is not None:
            # both default-direction conventions for cuts and bosses must give the same part
            ok1, m1 = verify_with_standin(code, res.shape, tol_, cut_reversed=True, boss_reversed=False)
            ok2, m2 = verify_with_standin(code, res.shape, tol_, cut_reversed=False, boss_reversed=True)
        if ok1 and ok2:
            f = os.path.join(out_dir, swname + ".py")
            with open(f, "w") as fh:
                fh.write(code)
            files.append(f)
            msg = f"checked on the CadQuery API stand-in ({m1}) ({time.time() - t0:.1f}s)"
        else:
            msg = "not written: " + (m1 if not ok1 else m2)
        log("SolidWorks script: " + msg)
        res.report += "\nSolidWorks script: " + msg
    f = os.path.join(out_dir, stem + "_report.txt")
    with open(f, "w") as fh:
        fh.write(res.report + "\n\noptions:\n" + o.to_json() + "\n\nsymmetry analysis:\n" +
                 str(res.extras.get("symmetry", "")) + "\n")
    files.append(f)
    log("saved: " + ", ".join(files))
    return files
