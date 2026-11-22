"""Symmetry detection and symmetry groups.

A symmetry is described in a frame (center c, orthonormal axes as columns of R).
Supported groups: any combination of mirror planes normal to the frame axes,
and cyclic symmetry C_n about one frame axis, optionally dihedral (D_n: mirror
planes containing the axis). Groups are lists of 4x4 world-space matrices.
"""

from dataclasses import dataclass, field

import numpy as np
import trimesh
from scipy.optimize import minimize, minimize_scalar


@dataclass
class SymmetrySpec:
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))
    axes: np.ndarray = field(default_factory=lambda: np.eye(3))   # columns: frame x, y, z
    mirrors: tuple = ()          # subset of (0, 1, 2): mirror planes normal to frame axes
    cyclic: int = 1              # order n of rotation symmetry
    cyclic_axis: int = 2         # frame axis of the rotation
    dihedral: bool = False       # add mirror planes containing the cyclic axis
    dihedral_angle: float = 0.0  # angle of the first such plane, radians, from the next frame axis

    def describe(self):
        parts = []
        if self.cyclic > 1:
            parts.append(f"{'D' if self.dihedral else 'C'}{self.cyclic} about frame axis "
                         f"{'xyz'[self.cyclic_axis]}")
        if self.mirrors:
            parts.append("mirror " + ",".join("xyz"[i] for i in self.mirrors))
        return " + ".join(parts) if parts else "none"


def _rot(axis, angle):
    return trimesh.transformations.rotation_matrix(angle, axis)[:3, :3]


def _mirror(normal):
    n = np.asarray(normal, float) / np.linalg.norm(normal)
    return np.eye(3) - 2 * np.outer(n, n)


def _world(M, spec):
    """Frame-space 3x3 linear map -> world 4x4 about spec.center."""
    R = spec.axes
    L = R @ M @ R.T
    T = np.eye(4)
    T[:3, :3] = L
    T[:3, 3] = spec.center - L @ spec.center
    return T


def group_matrices(spec):
    """All group elements as 4x4 world matrices (closure of the generators)."""
    gens = []
    e = np.eye(3)
    for i in spec.mirrors:
        gens.append(_mirror(e[i]))
    if spec.cyclic > 1:
        a = e[spec.cyclic_axis]
        gens.append(_rot(a, 2 * np.pi / spec.cyclic))
        if spec.dihedral:
            u = e[(spec.cyclic_axis + 1) % 3]
            v = np.cross(a, u)
            ang = spec.dihedral_angle
            gens.append(_mirror(np.cos(ang) * v - np.sin(ang) * u))  # plane contains a and (cos u + sin v)
    elems = [np.eye(3)]
    frontier = [np.eye(3)]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                k = h @ g
                if not any(np.allclose(k, x, atol=1e-9) for x in elems):
                    elems.append(k)
                    new.append(k)
        frontier = new
        if len(elems) > 96:
            break
    return [_world(M, spec) for M in elems]


def transform_points(T, pts):
    return pts @ T[:3, :3].T + T[:3, 3]


class SymmetryProbe:
    """Measures how well a mesh maps onto itself under a transform."""

    def __init__(self, mesh, n=4000, seed=1, n_ref=150000):
        from scipy.spatial import cKDTree
        self.mesh = mesh
        self.pts, _ = trimesh.sample.sample_surface_even(mesh, n, seed=seed)
        self.scale = float(np.linalg.norm(mesh.extents))
        self.max_frac = 0.02          # worst 0.5 % of points must stay within 2 % of the diagonal
        # dense reference sample with face normals: nearest sample + point-to-plane distance
        # replaces the exact closest-point query (about 100 times faster)
        ref, fid = trimesh.sample.sample_surface(mesh, n_ref, seed=seed + 1)
        self.ref, self.ref_n = ref, mesh.face_normals[fid]
        self.spacing = np.sqrt(mesh.area / n_ref)
        self.tree = cKDTree(ref)

    def distance(self, P):
        d, j = self.tree.query(P)
        dp = np.abs(np.einsum("ij,ij->i", P - self.ref[j], self.ref_n[j]))
        return np.where(d < 3 * self.spacing, np.minimum(d, dp), d)

    def error(self, T):
        """Mean distance; returns inf when the worst points are far off (outline mismatch)."""
        d = self.distance(transform_points(T, self.pts))
        self.last_max = float(np.percentile(d, 99.5))
        if self.last_max > self.max_frac * self.scale:
            return float(d.mean()) + self.last_max
        return float(d.mean())


def principal_frame(mesh, snap_deg=3.0):
    """Center (area centroid of the surface) and principal axes, snapped to world axes if close."""
    tri_c = mesh.triangles_center
    w = mesh.area_faces
    c = (tri_c * w[:, None]).sum(0) / w.sum()
    X = tri_c - c
    C = (X * w[:, None]).T @ X / w.sum()
    evals, evecs = np.linalg.eigh(C)
    R = evecs[:, ::-1].copy()                 # largest spread first
    # snap to world axes when nearly aligned or when the frame is degenerate
    snapped = np.zeros((3, 3))
    used = set()
    ok = True
    for j in range(3):
        k = int(np.argmax(np.abs(R[:, j])))
        if k in used or np.degrees(np.arccos(min(1.0, abs(R[k, j])))) > snap_deg:
            ok = False
            break
        used.add(k)
        snapped[k, j] = 1.0
    if ok:
        R = snapped
    if np.linalg.det(R) < 0:
        R[:, 2] *= -1
    return c, R, evals[::-1]


def detect(mesh, tol, log=print, max_order=12):
    """Detect the symmetry group of a mesh. Returns (SymmetrySpec, report string)."""
    c0, R, evals = principal_frame(mesh)
    probe = SymmetryProbe(mesh)
    lines = [f"frame center {np.round(c0, 4).tolist()}", "frame axes (columns):",
             str(np.round(R, 4))]
    thr = max(tol, 2e-3 * probe.scale)       # accept a symmetry if mean error is below this

    def refine_center_mirror(i, c):
        n = R[:, i]

        def f(s):
            spec = SymmetrySpec(center=c + s * n, axes=R, mirrors=(i,))
            return probe.error(group_matrices(spec)[1])
        r = minimize_scalar(f, bounds=(-0.05 * probe.scale, 0.05 * probe.scale), method="bounded",
                            options=dict(xatol=1e-4 * probe.scale))
        return c + r.x * n, r.fun

    def rot_err(axis, n, c):
        T = np.eye(4)
        T[:3, :3] = R @ _rot(np.eye(3)[axis], 2 * np.pi / n) @ R.T
        T[:3, 3] = c - T[:3, :3] @ c
        return probe.error(T)

    center = c0.copy()
    mirrors, merr = [], {}
    for i in range(3):
        cm, e = refine_center_mirror(i, center)
        merr[i] = e
        lines.append(f"mirror normal to frame {'xyz'[i]}: mean error {e:.4f} mm")
        if e < thr:
            mirrors.append(i)
            center = cm

    best = (1, 2, np.inf)
    for axis in range(3):
        u = [k for k in range(3) if k != axis]

        def f(s, axis=axis, n=2):
            c = center + s[0] * R[:, u[0]] + s[1] * R[:, u[1]]
            return rot_err(axis, n, c)
        errs = {}
        for n in range(2, max_order + 1):
            errs[n] = rot_err(axis, n, center)
        good = [n for n, e in errs.items() if e < thr]
        lines.append(f"rotation about frame {'xyz'[axis]}: " +
                     " ".join(f"C{n}:{e:.3f}" for n, e in errs.items()))
        if good:
            n = max(good)
            if n > best[0]:
                best = (n, axis, errs[n])
    spec = SymmetrySpec(center=center, axes=R, mirrors=tuple(mirrors))
    n, axis, e = best
    if n > 1:
        # refine the center in the plane normal to the rotation axis
        u = [k for k in range(3) if k != axis]

        def g(s):
            c = center + s[0] * R[:, u[0]] + s[1] * R[:, u[1]]
            return rot_err(axis, n, c)
        r = minimize(g, np.zeros(2), method="Nelder-Mead",
                     options=dict(xatol=1e-4 * probe.scale, fatol=1e-6, maxiter=80))
        center = center + r.x[0] * R[:, u[0]] + r.x[1] * R[:, u[1]]
        spec.center, spec.cyclic, spec.cyclic_axis = center, n, axis
        # mirrors normal to in-plane axes are absorbed by dihedral planes
        inplane = [i for i in spec.mirrors if i != axis]
        if inplane:
            spec.dihedral = True
            i = inplane[0]
            a = np.eye(3)[axis]
            uax = np.eye(3)[(axis + 1) % 3]
            v = np.cross(a, uax)
            nrm = np.eye(3)[i]
            # plane normal = cos*v - sin*u  =>  angle from components
            spec.dihedral_angle = float(np.arctan2(-nrm @ uax, nrm @ v))
            spec.mirrors = tuple(k for k in spec.mirrors if k == axis)
        else:
            # look for a dihedral mirror plane at another angle
            a = np.eye(3)[axis]
            uax = np.eye(3)[(axis + 1) % 3]
            v = np.cross(a, uax)

            def h(ang):
                M = _mirror(np.cos(ang) * v - np.sin(ang) * uax)
                T = _world(M, spec)
                return probe.error(T)
            angs = np.linspace(0, np.pi / n, 24, endpoint=False)
            vals = [h(x) for x in angs]
            k = int(np.argmin(vals))
            r = minimize_scalar(h, bounds=(angs[k] - np.pi / n / 24, angs[k] + np.pi / n / 24),
                                method="bounded")
            if r.fun < thr:
                spec.dihedral, spec.dihedral_angle = True, float(r.x)
                lines.append(f"dihedral mirror plane at {np.degrees(r.x):.2f} deg: error {r.fun:.4f}")
    G = group_matrices(spec)
    errs = [probe.error(T) for T in G[1:]]
    lines.append(f"detected: {spec.describe()}  (group order {len(G)}, "
                 f"max element error {max(errs) if errs else 0:.4f} mm, threshold {thr:.4f})")
    lines.append(f"center {np.round(spec.center, 4).tolist()}")
    report = "\n".join(lines)
    log(report)
    return spec, report


def spec_from_options(mesh, mirrors=(), cyclic=1, axis=2, dihedral=False, center=None, tol=None,
                      log=print):
    """Build a user-specified symmetry; center and in-plane orientation are refined on the mesh."""
    c0, R, _ = principal_frame(mesh)
    spec = SymmetrySpec(center=c0 if center is None else np.asarray(center, float), axes=R,
                        mirrors=tuple(mirrors), cyclic=int(cyclic), cyclic_axis=int(axis),
                        dihedral=bool(dihedral))
    if center is None:
        probe = SymmetryProbe(mesh)
        scale = probe.scale

        def total(x):
            s = SymmetrySpec(center=c0 + R @ x, axes=R, mirrors=spec.mirrors, cyclic=spec.cyclic,
                             cyclic_axis=spec.cyclic_axis, dihedral=False)
            G = group_matrices(s)
            return np.mean([probe.error(T) for T in G[1:]]) if len(G) > 1 else 0.0
        if len(group_matrices(spec)) > 1:
            r = minimize(total, np.zeros(3), method="Nelder-Mead",
                         options=dict(xatol=1e-4 * scale, fatol=1e-6, maxiter=150))
            spec.center = c0 + R @ r.x
    if spec.dihedral and spec.cyclic > 1:
        probe = SymmetryProbe(mesh)
        a = np.eye(3)[spec.cyclic_axis]
        uax = np.eye(3)[(spec.cyclic_axis + 1) % 3]
        v = np.cross(a, uax)

        def h(ang):
            return probe.error(_world(_mirror(np.cos(ang) * v - np.sin(ang) * uax), spec))
        angs = np.linspace(0, np.pi / spec.cyclic, 24, endpoint=False)
        k = int(np.argmin([h(x) for x in angs]))
        r = minimize_scalar(h, bounds=(angs[k] - np.pi / spec.cyclic / 24, angs[k] + np.pi / spec.cyclic / 24),
                            method="bounded")
        spec.dihedral_angle = float(r.x)
    G = group_matrices(spec)
    probe = SymmetryProbe(mesh)
    errs = [probe.error(T) for T in G[1:]]
    msg = (f"imposed symmetry: {spec.describe()} (order {len(G)}), center {np.round(spec.center, 4).tolist()}, "
           f"mesh departure mean {np.mean(errs) if errs else 0:.4f} / max {max(errs) if errs else 0:.4f} mm")
    log(msg)
    return spec
