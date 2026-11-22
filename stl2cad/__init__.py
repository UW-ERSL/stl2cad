"""stl2cad: reconstruct parametric CAD solids (STEP) from STL meshes.

Deterministic geometry processing only: RANSAC with a fixed seed, least squares,
region growing, voxel skeletons, and OpenCascade (via CadQuery) for the B-rep.
"""

__version__ = "0.1.0"
