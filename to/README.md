# TO reconstruction (topology-optimized parts)

SolidWorks scripts for topology-optimized (TO) parts, given the original part and the TO result as STL
(case 2). The part is built as plates of the original: each plate is the original between two levels along
an axis, minus cuts (voids and outline notches) fitted to the TO section with straight sides and tangent
fillets. Corner radii come from a few global values, linked in SolidWorks through the Equation Manager.

Modules (run from this folder; they import each other by name and need the stl2cad package installed):

| file | role |
|---|---|
| run_plates.py | driver |
| plates.py | plate detection (prism test, voxel run lengths), sections |
| tofit.py | cuts of one plate section; detail levels 1-5 |
| lf2.py | side and fillet fitting primitives |
| radii.py | global corner radii |
| codegen3d.py | CadQuery reference model (checks) |
| swplates.py | SolidWorks script and its check on the API stand-in |

Usage, from this folder:

    python run_plates.py ../models/to/EdgeCantilever.STL ../models/to/EdgeCantilever_TO_3.stl cantilever3_to --units m --sw --plot

Options: `--detail 1-5` (1: fewest independent dimensions, larger model error; 5: smallest model error;
default 4), `--radii N` (global corner radii; 0 = one per corner; default from the detail level),
`--fillet R` (radius on re-entrant edges where plates meet, default 1 mm), `--fit-only`.

Outputs in the current folder: NAME_solidworks.py, NAME.py and NAME.step (reference), NAME_standin.log,
NAME_sections.png. `run_all.sh` runs the pairs in ../models/to and copies the SolidWorks scripts to
../results/sw_scripts/to/.

Limits: plates are flat and axis-aligned; material at plate junctions (thickened nodes) is not modelled;
a TO section that does not reach the original outline is not handled; TO-only input (case 1) is not handled.

## 3D frames (hubs on several axes)

When the TO joins hubs on different axes with members running in 3D (the Knuckle: a clevis whose two eyes are on
the x axis and whose boss is on the y axis), neither plates nor the single-axis swept members apply. The part is
then built as a 3D frame:

| file | role |
|---|---|
| hubs3d.py | hubs on any axis: full 360 deg holes (fillets excluded), coaxial holes split, outermost TO ring |
| frame3d.py | 3D skeleton -> joints, straight tapered round members, sphere nodes; mirror symmetry; CadQuery model |
| swframe.py | SolidWorks script (3D skeleton sketch, lofts, revolved nodes, hubs and holes) and its stand-in check |
| sw_probe_frame.py | probe of the SolidWorks calls swframe adds (3D sketch, point-to-plane dimensions, ...) |
| sw_probe_bars.py | probe of the bar calls (3D spline, planes normal to it, lofts along it, fillet by CreateDefinition) |

Usage, from this folder (the original's exact model first, once):

    python -m stl2cad.cli ../models/to/Knuckle.STL --units m --no-solidworks -o out_Knuckle
    python swframe.py out_Knuckle/Knuckle_exact.py ../models/to/Knuckle_TO_1.stl knuckle_frame_d3 --units m --detail 3 --plot

Outputs: NAME_solidworks.py, NAME.step and NAME.stl (CadQuery model), NAME_standin.log, NAME_frame.png.
`--detail 1-5` sets the Douglas-Peucker tolerance of the members, the spur length pruned, and the number of
member and node radius classes (1: fewest dimensions ... 5: smallest model error).

Method (frame3d.fit): voxels of the TO; mirror planes through the original's box centre where the TO is
symmetric (>= 95 % overlap), imposed exactly; the medial curves of the TO outside the hubs as a graph (spurs
pruned, close joints merged); a chain whose inscribed tube comes within reach of a hub ends inside that hub, at a
point deep enough that the member's end section lies inside the hub; chains cut into straight members; joints on
the original's faces pushed inward until their node fits; member radii from the TO's distance transform,
scaled to the volume of the TO voxels nearest each member (area-equivalent), in classes.

SolidWorks model: mirror planes; one 3D sketch "Frame skeleton" with every joint dimensioned to three planes
(a mirror plane along a symmetric axis), so mirror images share one global variable per coordinate; hubs and
holes as extrusions on offset planes (coaxial hubs share their centre variables); nodes revolved on planes through
their skeleton points (they follow the joints); members as lofts between circles on planes normal to their
skeleton lines. Run `sw_probe_frame.py` once on a new SolidWorks version and check its log first.

Overlap (frame3d.JUNCTION = 'overlap', the default): no spheres. Each straight tapered member runs past its free
joints by r cot(a/2), a the smallest angle to another member there (at least 0.1 r, at most r and 0.4 x its
length), so neighbouring members overlap and the outside of a bend is filled; one fillet of radius R_FILLET (half
the smallest member radius) on every edge where members and hubs meet. SolidWorks: the section points stand on
the extensions of the skeleton lines; each member is a straight Loft (the calls of the sphere model, which built
in SolidWorks 2025). Joints are 3D sketch points edited by their X, Y, Z parameters (JOINT_DIMS = False in the
script; True restores the x, y dimensions and the ZLevel planes for z). Reference volumes: each feature adds its
exact volume (OCC's integration) minus its overlap with the earlier features whose boxes meet it, by a fine mesh
intersection (frame3d.steps_mesh); OCC's union of these frames was wrong by several mm^3 on some members (a Monte
Carlo check agreed with the mesh). The stand-in computes its volumes the same way for frame models.

Bars (frame3d.JUNCTION = 'bars'; 'spheres' gives the member and node model above): each chain of
members between junctions or hubs is one bar, a Loft along a 3D spline through its joints (CreateSpline2 with
natural ends, in the skeleton sketch). Its sections are circles on planes normal to the spline at the two ends
and at the interior joints where the radius departs from the linear change along the bar by more than a
tolerance set by `--detail`; between sections the radius changes linearly with the arc length, as SolidWorks
lofts it (`sw_probe_bars.py`: 0.03 % for two sections, 0.12 % for three). No spheres: the bars end at the
junction points, and one fillet of radius R_FILLET (half the smallest section radius) covers every edge where
bars and hubs meet. The reference volumes after each feature come from mesh Booleans (manifold3d) of the CadQuery
features with a first-order correction for the tessellation error, since OCC's Booleans fail where several
swept bars end at one point; the STEP file is then a compound of the features. In SolidWorks 2025 the loft of a
sharply bent bar (Knuckle level 3, Bar 4: end planes 96 deg apart) failed for every variant probed
(sw_probe_bar4*.py: with and without centerline, guide curve, Merge off, new circles, aligned connectors), so
bars are not the default.

Limits: members are round (a flat TO wall becomes a round member of equal area); a hollow hub in the TO (the
Knuckle's boss) is modelled solid; the hubs must come from full holes of the original; TO-only input is not
handled. The stand-in check of the bar model runs with a volume check after every feature, since OCC's fuse
fails on some bars there; the fillet is not checked by the stand-in.

