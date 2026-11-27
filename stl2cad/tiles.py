"""Repair of an STL exported in tiles (truss-connector): each tile is a separate shell, open where it was cut,
with seam vertices not shared between tiles. Each tile's planar boundary loops are capped (manifold3d.triangulate) and pushed
1 um past the cut so neighbouring tiles overlap; each tile is decimated (quadric, 12 % of its faces); the
tiles are united (manifold3d) and the largest body is kept."""
import numpy as np
import trimesh


def _cap(b, overlap=1e-3, log=print):
    """The tile with its planar boundary loops capped (cap triangles facing out of the tile). Needs manifold3d
    (also used for the union) and networkx."""
    import manifold3d
    import networkx as nx
    b.merge_vertices()
    u, c = np.unique(b.edges_sorted, axis=0, return_counts=True)
    bd = u[c == 1]
    if len(bd) == 0:
        return b
    G = nx.Graph()
    G.add_edges_from(map(tuple, bd))
    caps = []
    for comp in nx.connected_components(G):
        try:
            cyc = nx.cycle_basis(G.subgraph(comp))
        except Exception:
            cyc = []
        for cy in cyc:
            P = b.vertices[cy]
            ax = int(np.argmin(np.ptp(P, axis=0)))           # the cut plane's normal axis
            if np.ptp(P[:, ax]) > 1e-3:
                log(f"  tiles: non-planar boundary loop {np.ptp(P, axis=0)}; left open")
                continue
            idx = [k for k in range(3) if k != ax]
            xc = P[:, ax].mean()
            out_sign = 1.0 if b.vertices[:, ax].mean() < xc else -1.0
            xc = xc + out_sign * overlap                      # overlap the neighbouring tile
            b.vertices[cy, ax] = xc
            V2 = P[:, idx]
            x, y = V2[:, 0], V2[:, 1]
            if (x * np.roll(y, -1) - np.roll(x, -1) * y).sum() < 0:    # counter-clockwise for triangulate
                V2 = V2[::-1]
            F = np.asarray(manifold3d.triangulate([V2]), int).reshape(-1, 3)
            V3 = np.zeros((len(V2), 3))
            V3[:, idx] = V2
            V3[:, ax] = xc
            # the tile's own triangles keep the STL's winding, so the caps are oriented directly (a global
            # normal repair took 90 % of the run)
            nrm = np.cross(V3[F[:, 1]] - V3[F[:, 0]], V3[F[:, 2]] - V3[F[:, 0]])[:, ax]
            F = np.where((nrm * out_sign < 0)[:, None], F[:, ::-1], F)
            caps.append(trimesh.Trimesh(V3, F, process=False))
    t = trimesh.util.concatenate([b] + caps)
    t.merge_vertices(digits_vertex=6)
    if not t.is_winding_consistent or t.volume < 0:
        trimesh.repair.fix_normals(t)
    return t


def is_tiled(mesh):
    """Several shells, at least one open: an STL exported in tiles."""
    parts = mesh.split(only_watertight=False)
    return len(parts) > 1 and not all(p.is_watertight for p in parts)


def repair_tiles(mesh, keep=0.12, log=print):
    parts = [b for b in mesh.split(only_watertight=False) if len(b.faces) > 10]
    dec = []
    for b in parts:
        t = _cap(b, log=log)
        t2 = t.simplify_quadric_decimation(face_count=max(2000, int(len(t.faces) * keep)))
        dec.append(t2 if t2.is_watertight else t)
    U = trimesh.boolean.union(dec, engine="manifold")
    bodies = U.split(only_watertight=False)
    if len(bodies) > 1:
        U = max(bodies, key=lambda x: abs(x.volume))
    log(f"  tiles: {len(parts)} tile(s) capped and united: {len(U.faces)} faces, watertight {U.is_watertight}, "
        f"volume {U.volume:.3f} mm^3")
    return U
