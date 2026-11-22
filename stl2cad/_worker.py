"""Child-process entry for OCC operations that can crash or exhaust memory.

    python -m stl2cad._worker split <in.brep> <out.brep> <fuzzy> <nondestructive>
    python -m stl2cad._worker clean <in.brep> <out.brep> <mem_bytes>
"""

import sys


def main():
    op, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    import cadquery as cq
    if op == "split":
        from .brep_helpers import _split
        fuzzy, nd = float(sys.argv[4]), sys.argv[5] == "1"
        parts = list(cq.Shape.importBrep(src))
        cells = _split(parts[0], [t.wrapped for t in parts[1:]], fuzzy, nd)
        if cells:
            cq.Compound.makeCompound(cells).exportBrep(dst)
    elif op == "clean":
        try:
            import resource
            mem = int(sys.argv[4])
            resource.setrlimit(resource.RLIMIT_AS, (mem, mem))
        except Exception:
            pass
        from cadquery import func as cf
        cf.clean(cq.Shape.importBrep(src)).exportBrep(dst)


if __name__ == "__main__":
    main()
