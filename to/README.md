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

Limits: members are round (a flat TO wall becomes a round member of equal area); junctions are spheres, so the
re-entrant edges where members meet nodes and hubs are not filleted; a hollow hub in the TO (the Knuckle's boss)
is modelled solid; the hubs must come from full holes of the original; TO-only input is not handled.

