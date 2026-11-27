import sys, numpy as np, trimesh
from shapely.geometry import Polygon
m = trimesh.load(sys.argv[1])
parts = [b for b in m.split(only_watertight=False) if len(b.faces) > 10]
out = []
for b in parts:
    b.merge_vertices()
    e = b.edges_sorted; u, c = np.unique(e, axis=0, return_counts=True)
    bd = u[c == 1]
    if len(bd) == 0:
        out.append(b); continue
    import networkx as nx
    G = nx.Graph(); G.add_edges_from(map(tuple, bd))
    caps = []
    for comp in nx.connected_components(G):
        H = G.subgraph(comp)
        try:
            cyc = nx.cycle_basis(H)
        except Exception:
            cyc = []
        for cy in cyc:
            P = b.vertices[cy]
            ax = int(np.argmin(np.ptp(P, axis=0)))           # the cut plane's normal axis
            if np.ptp(P[:, ax]) > 1e-3:
                print("non-planar loop", np.ptp(P, axis=0)); continue
            idx = [k for k in range(3) if k != ax]
            xc = P[:, ax].mean()
            out_sign = 1.0 if b.vertices[:, ax].mean() < xc else -1.0
            xc = xc + out_sign * 1e-3                      # overlap the neighbouring tile slightly
            b.vertices[cy, ax] = xc
            poly = Polygon(P[:, idx])
            if not poly.is_valid: poly = poly.buffer(0)
            V2, F = trimesh.creation.triangulate_polygon(poly, engine="earcut")
            V3 = np.zeros((len(V2), 3)); V3[:, idx] = V2; V3[:, ax] = xc
            cap = trimesh.Trimesh(V3, F)
            caps.append(cap)
    t = trimesh.util.concatenate([b] + caps)
    t.merge_vertices(digits_vertex=6)
    trimesh.repair.fix_normals(t)
    out.append(t)
    print(len(b.faces), "caps", len(caps), "watertight", t.is_watertight, "vol", round(t.volume, 2))
dec = []
for t in out:
    t2 = t.simplify_quadric_decimation(face_count=max(2000, int(len(t.faces) * 0.12)))
    dec.append(t2 if t2.is_watertight else t)
U = trimesh.boolean.union(dec, engine="manifold")
print("union", len(U.faces), U.is_watertight, "bodies", U.body_count, "volume", U.volume)
U.export(sys.argv[2])
