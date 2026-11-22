"""Command-line interface: stl2cad-cli part.stl [options]."""

import argparse

from .pipeline import Options, run


def main(argv=None):
    o = Options()
    ap = argparse.ArgumentParser(description="Reconstruct a parametric CAD solid (STEP) from an STL mesh.")
    ap.add_argument("stl")
    ap.add_argument("-o", "--out-dir", default=".")
    ap.add_argument("--units", default=o.units, choices=["mm", "cm", "m", "in"])
    ap.add_argument("--mode", default=o.mode, choices=["auto", "exact", "smooth"])
    ap.add_argument("--tol", type=float, default=0.0, help="fit tolerance in mm (0 = auto)")
    ap.add_argument("--symmetry", default=o.symmetry, choices=["off", "auto", "manual"])
    ap.add_argument("--mirror", default="", help="mirror planes to impose, e.g. 'xy' (manual symmetry)")
    ap.add_argument("--cyclic", type=int, default=1, help="cyclic order to impose (manual symmetry)")
    ap.add_argument("--axis", default="z", choices=["x", "y", "z"], help="frame axis of the rotation")
    ap.add_argument("--dihedral", action="store_true")
    ap.add_argument("--algebraic-only", action="store_true",
                    help="exact: no spline extrusions; smooth: members from cones and spheres")
    ap.add_argument("--surfaces", default=",".join(o.surfaces), help="exact mode surface types")
    ap.add_argument("--no-extrusion", action="store_true")
    ap.add_argument("--no-tangency", action="store_true")
    ap.add_argument("--no-fillets", action="store_true", help="do not rebuild tangent blends as fillets")
    ap.add_argument("--no-prism", action="store_true", help="do not try the stacked-prism (2.5D) method")
    ap.add_argument("--no-solidworks", action="store_true", help="do not write the SolidWorks feature-tree script")
    ap.add_argument("--no-script", action="store_true", help="do not write a CadQuery script")
    ap.add_argument("--snap-deg", type=float, default=o.snap_directions_deg)
    ap.add_argument("--angle-step", type=float, default=o.snap_angle_step_deg)
    ap.add_argument("--grid", type=float, default=0.0, help="round dimensions to this grid (mm)")
    ap.add_argument("--pitch", type=float, default=0.0, help="smooth: voxel size (mm)")
    ap.add_argument("--blur", type=float, default=0.0, help="smooth: smoothing radius (mm)")
    ap.add_argument("--prune", type=float, default=o.prune)
    ap.add_argument("--min-copy-frac", type=float, default=o.min_copy_frac)
    ap.add_argument("--no-holes", action="store_true")
    ap.add_argument("--drop-free-ends", action="store_true")
    ap.add_argument("--clip", default=o.clip, choices=["none", "thin", "all"])
    ap.add_argument("--no-stl", action="store_true")
    a = ap.parse_args(argv)
    o = Options(units=a.units, mode=a.mode, tolerance=a.tol, symmetry=a.symmetry,
                mirror_x="x" in a.mirror, mirror_y="y" in a.mirror, mirror_z="z" in a.mirror,
                cyclic=a.cyclic, cyclic_axis=a.axis, dihedral=a.dihedral,
                surfaces=[s for s in a.surfaces.split(",") if s], algebraic_only=a.algebraic_only,
                allow_extrusion=not a.no_extrusion, snap_directions_deg=a.snap_deg,
                snap_angle_step_deg=a.angle_step, tangency=not a.no_tangency, fillets=not a.no_fillets, prism25=not a.no_prism, export_script=not a.no_script, export_solidworks=not a.no_solidworks, grid=a.grid,
                pitch=a.pitch, blur=a.blur, prune=a.prune, min_copy_frac=a.min_copy_frac, holes=not a.no_holes,
                keep_free_ends=not a.drop_free_ends, clip=a.clip, out_dir=a.out_dir, export_stl=not a.no_stl)
    if (a.mirror or a.cyclic > 1) and a.symmetry == "auto":
        o.symmetry = "manual"
    run(a.stl, o)


if __name__ == "__main__":
    main()
