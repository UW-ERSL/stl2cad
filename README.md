# stl2cad

Reconstruct a parametric CAD solid (STEP) from an STL mesh, with a Tkinter GUI
and a command-line interface. The program uses deterministic geometry
processing only: least-squares fitting, RANSAC with a fixed seed, region
growing, voxel skeletons, and OpenCascade (through CadQuery) for the B-rep.
It contains no machine-learning models.

Two reconstruction modes are provided:

1. **Exact** (CAD-like parts): fits planes, cylinders, cones, spheres, tori,
   and optionally spline extrusions, then builds an exact B-rep bounded only by
   the fitted surfaces.
2. **Smooth** (topology-optimized and other organic parts): builds a smooth,
   symmetric approximation from the medial skeleton, with exact bolt holes and
   pads.

**Auto** mode chooses between them from the share of surface area on exactly
coplanar neighbouring triangles. This share is above 0.8 for CAD exports and
about 0.02 for the topology-optimized rotor tested here.

## Installation

```
pip install -e .            # Python 3.10+; installs cadquery, trimesh, scipy, scikit-image, ...
```

The GUI needs Tkinter, which ships with the python.org and conda builds of
Python. On Debian or Ubuntu system Python, install it with
`sudo apt install python3-tk`.

## Usage

```
stl2cad                                   # GUI (or: python -m stl2cad)
stl2cad-cli part.stl --units in -o out/   # command line, automatic everything
stl2cad-cli rotor.stl --units m --mode smooth --cyclic 6 -o out/
stl2cad-cli rotor.stl --units m --mode smooth --cyclic 6 --algebraic-only -o out/
stl2cad-cli bracket.stl --mode exact --mirror y --grid 0.5 -o out/
```

Each run writes `<name>_<mode>.step`, a tessellated `.stl` of the result, and
`_report.txt` with the deviation statistics, the options used, and the
symmetry analysis. In exact mode it also writes `<name>_exact.py`, a CadQuery
script that rebuilds the solid (see below); `--no-script` turns this off.

## CadQuery script output (`codegen.py`)

A script is written only after it has been run in a fresh interpreter and its
solid compared with the reconstruction: equal volume within 1e-4, and surfaces
within half the tolerance (area-weighted samples projected onto the exact
surfaces, then measured with `BRepExtrema`). Two forms are tried in order:

1. **2.5D feature script** (results of the stacked-prism method below): the
   part as a sequence of `prism(outline, lo, hi)` extrusions and cuts in the
   order they are applied, with named levels `z0, z1, ...`, round holes as a
   table with named diameters, and an outline that repeats about a centre (gear
   teeth) written as one period plus a `count`. About 50 to 100 lines.
2. **Feature script** (prisms): when two parallel, world-aligned planar caps
   bound the part and every other face is a plane or cylinder parallel to the
   extrusion direction spanning the full thickness. The script uses
   `cadquery.func`: the cap outline as `cf.segment` and three-point arcs,
   `cf.extrude` by a named `thickness`, and round through-holes listed with
   named diameters. About 40 to 90 lines; edit any number and re-run.
3. **Construction script** (everything else): the fitted surfaces as an
   editable table (`SURFACES`), the bounding box, one interior point per kept
   cell (`INSIDE`), the fillet operations (`FILLETS`: radius and edge midpoint)
   and the exact tori, followed by verbatim copies of the construction helpers
   from `brep_helpers.py`. It needs CadQuery only. Values are kept at full
   precision: tangent and coincident surfaces depend on them, and rounding even
   to 1e-6 mm changes the arrangement. Edit a radius or an offset and re-run;
   keep each edit small enough that the `INSIDE` points stay inside the part.

Run any of them as `python <name>_exact.py [out.step]`.

## SolidWorks feature tree (`swgen.py`, 2.5D results)

For a part built as stacked prisms, `<name>_solidworks.py` builds the same part
in SolidWorks through its COM API (Windows, SolidWorks, `pip install pywin32`):
one sketch and one Extruded Boss or Extruded Cut per feature, in order, with the
round holes of one height range in one cut. The result is a native feature tree.
Run `python <name>_solidworks.py [C:\path\part.SLDPRT]` with SolidWorks installed.

Sketch points use SolidWorks' base-plane orientations (Front (x, y), Top
(x, −z), Right (−z, y)); the active sketch object is not reliably available over
late-bound COM. The script does not rely on the default extrusion direction:
after each feature the part's volume and centre of mass are compared with the
values stl2cad computed; if they differ, the feature is deleted and rebuilt with
the other direction or start-offset side, and the script stops with an error if
no combination matches.

Observed in SolidWorks 2025 over pywin32: `SketchManager.CreateLine` creates nothing until
another entity exists in the sketch, and the active sketch object is not returned. Each sketch
therefore starts with a 0.5 mm helper circle, deleted at the end (made construction geometry if
the delete fails); `swmock` reproduces both behaviours.

Check before writing: the script is run against `swmock`, a CadQuery stand-in for
the subset of the SolidWorks API it uses, under both conventions for the default
cut and boss direction, and the body is compared with the reconstruction. This
checks the script's geometry and direction handling, not SolidWorks itself; the
scripts have not yet been run in SolidWorks. Parts built by the surface
arrangement have no feature sequence and get no SolidWorks script.
Hole sketches get a diameter dimension per hole and the horizontal and vertical distance of its
centre from the origin (selected by coordinates; the Input Dimension Value dialog is switched off
during the run and restored). The stand-in checks the dimension values.
CLI `--no-solidworks`, GUI checkbox.

## Options

| Option | Meaning |
|---|---|
| Units | Units of the STL (mm, cm, m, in). Everything internal and the STEP file are in mm. |
| Mode | auto / exact / smooth |
| Tolerance | Fit tolerance in mm; 0 = 0.05 % of the bounding-box diagonal |
| Symmetry | off / auto-detect / impose. Impose takes mirror planes (frame x, y, z), cyclic order *n* about a frame axis, and a dihedral flag. |
| Algebraic surfaces only | Exact: no spline extrusions. Smooth: members are built from cones and spheres instead of lofted B-spline tubes. |
| Exact: surface types | Any subset of plane, cylinder, cone, sphere, torus, plus spline extrusions |
| Exact: snapping | Surface directions are snapped to the part frame (and to a 15° lattice) within the given angle. Tangency is enforced between cylinders and planes. Dimensions can optionally be rounded to a grid. |
| Smooth: pitch, blur | Voxel size and Gaussian smoothing radius (0 = automatic) |
| Smooth: prune | Removes skeleton spurs shorter than this multiple of the local radius |
| Smooth: min copy fraction | A member is kept if it is found in at least this fraction of the symmetric copies |
| Smooth: holes and pads | Detects concave cylinders (bolt holes) and rebuilds each as an exact hole in a round pad |
| Smooth: trim | Trims the result to the design-domain slab (`thin`: the thinnest frame axis only) or to the full bounding box (`all`) |

How symmetry is used depends on the mode. In smooth mode, symmetry is
enforced exactly. If the detected symmetry is dihedral (D_n), automatic mode
imposes only its rotational part C_n and reports the mirror planes. Imposing
D_n manually works, but it is experimental: on the brake rotor, D_n gave
−17 % volume, against +0.8 % for C_n. In exact mode, a detected symmetry is
only reported, because exact mode reproduces the part as it is; choose
**Impose** to enforce it. Imposing completes every surface orbit, refits each
surface on all of its symmetric copies, and classifies each cell by majority
vote over the group.

## Method

### Stacked prisms, 2.5D (`prism25.py`)

Tried first in exact mode (`--no-prism` to skip). A part qualifies when at
least 99.9 % of its surface area faces along one world axis or across it.
1. **Levels:** the heights of the cap faces, clustered within the tolerance.
2. **Sections:** one section per slab between levels, as closed loops of mesh
   vertices (exact surface points), not of plane crossings (chord points).
   Loops thinner than the tolerance are gaps between faceted shells and are
   dropped.
3. **Loop fitting:** greedy lines and circular arcs within the tolerance; an
   arc needs facet steps of similar size, so a far point is not absorbed. A loop
   that goes once around a centre and stays within the tolerance of a circle
   (after dropping chord points) is a full circle. A loop unchanged by a
   rotation of 2π/N (N from the strongest harmonic of radius against angle,
   then confirmed point by point) is fitted once per period on the points of
   all periods pooled.
4. **Features:** each loop is followed through consecutive slabs while it stays
   the same, giving one extrusion over its whole height range. Its nesting
   depth decides whether it adds (even) or removes (odd) material; features are
   applied in depth order, which reproduces every section.
5. **Consolidation:** circle and pattern centres within the tolerance become
   one centre, radii within the tolerance one radius.
6. **Snapping (`snap2d.py`):** lines within tolerance of horizontal or vertical become
   exactly so; line offsets, centres and radii are clustered and rounded to the coarsest clean
   value within tol/4, on mm or inch grids chosen per part; arcs tangent to axis-parallel lines
   are made exactly tangent; corners are recomputed. A loop keeps its snapped form only if every
   section point lies within the tolerance of it; a part that then fails the check is rebuilt
   unsnapped.
7. **Check:** valid single solid, every mesh vertex within 1.5 tol of the model,
   the model within 4 tol of the mesh (facet sag plus fitting error), volume
   within 3 % (a faceted thin ring is 1 % smaller than its true circles).
   Otherwise the surface arrangement below is used.

An STL written as several overlapping closed shells is first merged into their
union (`manifold3d`).

### Exact mode (`segment.py`, `exact.py`)

1. **Segmentation:** efficient-RANSAC-style candidate generation from
   oriented minimal samples: one point for a plane, two for a cylinder or
   sphere, three for a cone, and a local Plücker normal-line fit for a torus.
   Each candidate's largest connected inlier region is scored by area, then
   refit by least squares. Inliers are tested on vertices and on face
   centroids, with an allowance for the chord sag of correctly tessellated
   curved surfaces.
2. **Clean-up:** small fragments are absorbed by larger surfaces that already
   explain them. Region competition reassigns boundary faces to the
   best-fitting neighbour. Poorly fitting curved fragments are merged into one
   cylinder, or into one spline extrusion when extrusions are allowed.
3. **Snapping and cleanup:** directions are snapped to the frame, coplanar and
   coaxial surfaces are unified, and cylinder–plane tangency is enforced by
   constrained least squares. Spline revolutions and extrusions whose profile
   consists of lines and circular arcs are decomposed into planes, cylinders,
   cones, tori and spheres; coaxial cones lying on one torus become that torus;
   regions on one curved surface are merged.
4. **Blends:** a cylinder or torus spanning less than a half turn whose spine
   lies at distance r from two neighbouring surfaces is a fillet. It is removed,
   its neighbours are extended until they meet, and after the sharp solid is
   built the fillet is re-created on that edge with the OCC fillet operation and
   the fitted radius (B-spline fillet faces on a torus are rebuilt as exact
   tori). Plane–cylinder tangencies parallel to the cylinder axis are handled
   by tangent chains instead (step 5). Option `--no-fillets` turns this off.
5. **B-rep by cell arrangement:** every surface becomes a tool face extended
   a margin past its region. Chains of mutually tangent planes and cylinders
   become one extruded line-and-arc profile shell, because separate faces
   leave gaps at tangent contacts. A box is split by all tool faces
   (`BOPAlgo_Splitter`), and each cell is classified inside or outside the
   STL. The true boundary lies on the union of the tool faces, so every cell
   is entirely inside or entirely outside. The inside cells are glued and
   their faces unified.

### Smooth mode (`smooth.py`)

1. Voxelize the part, blur it, and average it over the symmetry group.
2. Skeletonize the result into a graph of junctions and members.
3. Detect holes with concave-only cylinder RANSAC. Each hole gets a pad whose
   radius is the 20th percentile of the material ring around it.
4. Map each member into a fundamental domain of the group (a Dirichlet-style
   choice of group element), cluster the copies of each member, and average
   them.
5. Size each member to the median cross-sectional area of its copies in the
   original mesh.
6. Build members as lofted tubes along smoothing splines, or as cone chains
   in algebraic-only mode. Junctions become spheres.
7. Pattern the domain by the group, trim to the design domain, and cut the
   holes.

## Test results

Exact mode, default tolerance (5e-4 of the part diagonal); run times with 3 to 4
runs in parallel. Script "2.5D" means the stacked-prism method built the part.
*Gear: model to mesh; the mesh to model maximum (10 mm) is at thin slits
between the faceted web and the faceted ring of the input, which the model
fills.

| Part | Triangles | Faces | Max deviation | Volume | Script | Time |
|---|---|---|---|---|---|---|
| C-bracket | 56 | 8 | 0.038 mm | −0.01 % | 2.5D | 7 s |
| HeatBlock | 186 | 15 | 0.047 mm | +0.00 % | construction | 22 s |
| BicycleCrank | 280 | 16 | 0.057 mm | −0.01 % | construction | 21 s |
| SimpleBracket | 308 | 8 | 0.011 mm | −0.01 % | 2.5D | 11 s |
| CircularPlateHole | 392 | 4 | 0.120 mm | −0.10 % (faceted rim) | 2.5D | 8 s |
| SwingArm | 396 | 16 | 0.046 mm | −0.01 % | 2.5D | 14 s |
| ThreeHoleBracket | 452 | 8 | 0.019 mm | −0.01 % | 2.5D | 19 s |
| CrankArm | 588 | 11 | 0.058 mm | +0.04 % | 2.5D | 15 s |
| angle_block (trimesh) | 704 | 26 | 0.026 mm | −0.02 % | construction | 29 s |
| idler_riser (trimesh) | 1,572 | 42 | 0.115 mm | +0.13 % | construction | 37 s |
| plate_holes (trimesh) | 1,252 | 29 (7 fillets) | 0.237 mm (coarse seams) | −0.03 % | construction | 35 s |
| 1002_tray_bottom (trimesh) | 4,520 | 48 | 0.187 mm | −0.01 % | 2.5D | 24 s |
| Gear (3 overlapping shells) | 7,300 | 302 (48 teeth) | 0.214 mm* | +0.10 % | 2.5D | 25 s |
| GE jet-engine bracket | 5,448 | 58 (4 fillets) | 0.079 mm | +0.00 % | construction | 79 s |
| AlcoaGrabCAD | 6,754 | 52 | 0.062 mm | +0.01 % | construction | 69 s |
| featuretype (trimesh) | 14,954 | 79 | 0.076 mm | +0.03 % | construction | 96 s |

Smooth mode: the topology-optimized brake rotor (28,696 triangles), C6 imposed,
18 holes and pads found, mean deviation 0.86 mm, volume +1.0 %, 138 s.

Failures: BatteryHolder (blend network), the arrowHead lattice (needs a lattice
mode); the exact brakeRotor is not yet resolved.

## Limitations

1. **Blends:** constant-radius cylinder and torus fillets are handled; branched
   networks of plane–cylinder fillets (several fillets meeting at a corner, as
   in a filleted shell) and variable-radius blends are not. A missed surface
   leaves a gap in the arrangement, and the check then reports a failure.
2. **Freeform surfaces:** only linear extrusions of splines are supported.
   General B-spline patches are not fitted, so in exact mode such regions fall
   back to planar facets.
3. **Validity:** occasionally the result is geometrically correct but OCC's
   checker flags it (as for `idler_riser`). The STEP usually still imports,
   but inspect it before downstream use.
4. **Speed:** symmetry detection and hole detection on large meshes (about
   30k triangles) each take one to two minutes. Smooth mode on the rotor takes
   about three minutes in total.
5. **Smooth mode is an approximation by design:** member cross-sections are
   circular (area-equivalent), and junction geometry is spherical.
