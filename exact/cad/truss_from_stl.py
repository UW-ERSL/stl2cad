"""STL of a truss (truss-connector) to the CadQuery script, STEP and SolidWorks script, with the time of each
step: tile repair (if the STL is tiled), strut recognition (stl2cad.struts.detect), script generation.
    python truss_from_stl.py ../truss-connector.stl truss-connector"""
import os
import sys
import time

import trimesh

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from stl2cad import struts, tiles                      # noqa: E402
from stl2cad.core import default_tolerance             # noqa: E402


def check_packages():
    """The packages the chain needs, checked before the run (pip names)."""
    import importlib.util
    need = {"networkx": "networkx", "manifold3d": "manifold3d", "fast_simplification": "fast-simplification",
            "rtree": "rtree", "scipy": "scipy", "cadquery": "cadquery"}
    missing = [pip for mod, pip in need.items() if importlib.util.find_spec(mod) is None]
    if missing:
        sys.exit("missing packages: pip install " + " ".join(missing))
    if importlib.util.find_spec("embreex") is None:
        print("embreex not installed: rays use trimesh's own intersector (about 2 minutes more); "
              "pip install embreex")


def main():
    check_packages()
    src, name = sys.argv[1], sys.argv[2]
    T = t = time.time()

    def lap(label):
        nonlocal t
        print(f"[{label}: {time.time() - t:.1f} s]")
        t = time.time()

    mesh = trimesh.load(src)
    if tiles.is_tiled(mesh):
        mesh = tiles.repair_tiles(mesh)
        mesh.export(f"{name}-clean.stl")
        lap("tile repair")
    tol = default_tolerance(mesh)
    print(f"tolerance {tol:.4f} mm")
    res = struts.detect(mesh, tol)
    if res is None:
        sys.exit("no truss of straight struts found")
    lap("recognition")
    open(f"{name}_exact.py", "w").write(struts.cadquery_code(res, name))
    shape = struts.build(res, log=lambda *a: None)
    shape.exportStep(f"{name}_exact.step")
    print(f"CadQuery model: valid {shape.isValid()}, volume {shape.Volume():.3f} mm^3 (mesh {mesh.volume:.3f})")
    lap("CadQuery script and STEP")
    open(f"{name}_solidworks.py", "w").write(struts.solidworks_code(res, name))
    lap("SolidWorks script")
    print(f"total {time.time() - T:.1f} s")


if __name__ == "__main__":     # required: the segmentation runs in worker processes (spawned on Windows)
    main()
