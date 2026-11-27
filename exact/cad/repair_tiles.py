"""Repair of a tiled STL (stl2cad.tiles): python repair_tiles.py in.stl out.stl"""
import os
import sys

import trimesh

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from stl2cad.tiles import repair_tiles                 # noqa: E402

if __name__ == "__main__":
    repair_tiles(trimesh.load(sys.argv[1])).export(sys.argv[2])
