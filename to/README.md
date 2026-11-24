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
