"""Tkinter GUI for stl2cad."""

import os
import queue
import threading
import traceback
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np

from .pipeline import Options, coplanar_fraction, run, save

PREVIEW_FACES = 12000


class App(ttk.Frame):
    def __init__(self, master):
        super().__init__(master, padding=6)
        self.master = master
        master.title("stl2cad - STL to parametric CAD")
        master.geometry("1320x860")
        self.q = queue.Queue()
        self.path = None
        self.mesh = None
        self.result = None
        self.busy = False
        self.v = {}
        self._build()
        self.grid(sticky="nsew")
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        self.after(100, self._poll)

    # ---- layout -------------------------------------------------------------
    def _var(self, name, kind, value):
        var = {"s": tk.StringVar, "b": tk.BooleanVar, "d": tk.DoubleVar, "i": tk.IntVar}[kind](value=value)
        self.v[name] = var
        return var

    def _build(self):
        top = ttk.Frame(self)
        top.grid(row=0, column=0, columnspan=2, sticky="ew")
        ttk.Button(top, text="Open STL...", command=self.open_file).pack(side="left")
        self.file_label = ttk.Label(top, text="no file", width=48)
        self.file_label.pack(side="left", padx=6)
        ttk.Label(top, text="Units").pack(side="left")
        ttk.Combobox(top, textvariable=self._var("units", "s", "mm"), values=["mm", "cm", "m", "in"],
                     width=5, state="readonly").pack(side="left", padx=4)
        ttk.Label(top, text="Output folder").pack(side="left", padx=(12, 2))
        ttk.Entry(top, textvariable=self._var("out_dir", "s", os.path.expanduser("~")), width=34).pack(side="left")
        ttk.Button(top, text="...", width=3, command=self.pick_dir).pack(side="left")

        left = ttk.Frame(self)
        left.grid(row=1, column=0, sticky="nsw", pady=6)
        nb = ttk.Notebook(left)
        nb.pack(fill="both", expand=True)
        nb.add(self._tab_general(nb), text="General")
        nb.add(self._tab_symmetry(nb), text="Symmetry")
        nb.add(self._tab_exact(nb), text="Exact")
        nb.add(self._tab_smooth(nb), text="Smooth (TO)")
        btns = ttk.Frame(left)
        btns.pack(fill="x", pady=6)
        self.btn_an = ttk.Button(btns, text="Analyze", command=self.analyze)
        self.btn_an.pack(side="left", expand=True, fill="x")
        self.btn_run = ttk.Button(btns, text="Run", command=self.run)
        self.btn_run.pack(side="left", expand=True, fill="x", padx=4)
        self.btn_save = ttk.Button(btns, text="Save", command=self.save, state="disabled")
        self.btn_save.pack(side="left", expand=True, fill="x")
        self.progress = ttk.Progressbar(left, mode="indeterminate")
        self.progress.pack(fill="x")

        right = ttk.Frame(self)
        right.grid(row=1, column=1, sticky="nsew", padx=(8, 0), pady=6)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure
        self.fig = Figure(figsize=(8, 4.6), dpi=90)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.log_box = tk.Text(right, height=14, font=("Courier", 9), wrap="none")
        self.log_box.pack(fill="both", expand=False)

    def _row(self, parent, r, label, widget):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=2)
        widget.grid(row=r, column=1, sticky="ew", pady=2)

    def _tab_general(self, nb):
        f = ttk.Frame(nb, padding=8)
        ttk.Label(f, text="Mode").grid(row=0, column=0, sticky="w")
        mf = ttk.Frame(f)
        mf.grid(row=0, column=1, sticky="w")
        mode = self._var("mode", "s", "auto")
        for t, val in (("Auto", "auto"), ("Exact (CAD-like)", "exact"), ("Smooth (TO / organic)", "smooth")):
            ttk.Radiobutton(mf, text=t, variable=mode, value=val).pack(anchor="w")
        self._row(f, 1, "Tolerance [mm] (0 = auto)", ttk.Entry(f, textvariable=self._var("tolerance", "d", 0.0), width=10))
        ttk.Checkbutton(f, text="Algebraic surfaces only", variable=self._var("algebraic_only", "b", False)).grid(
            row=2, column=0, columnspan=2, sticky="w", pady=4)
        msg = ("Exact: fits planes, cylinders, cones, spheres, tori (and optionally spline\n"
               "extrusions) and builds an exact B-rep.\n"
               "Smooth: skeleton-based members with symmetric averaging, for\n"
               "topology-optimized parts. Algebraic-only uses cones and spheres.")
        ttk.Label(f, text=msg, foreground="#555").grid(row=3, column=0, columnspan=2, sticky="w", pady=8)
        return f

    def _tab_symmetry(self, nb):
        f = ttk.Frame(nb, padding=8)
        sym = self._var("symmetry", "s", "auto")
        for i, (t, val) in enumerate((("Off", "off"), ("Auto-detect", "auto"), ("Impose (below)", "manual"))):
            ttk.Radiobutton(f, text=t, variable=sym, value=val).grid(row=0, column=i, sticky="w")
        ttk.Separator(f).grid(row=1, column=0, columnspan=3, sticky="ew", pady=6)
        for i, ax in enumerate("xyz"):
            ttk.Checkbutton(f, text=f"Mirror {ax}", variable=self._var(f"mirror_{ax}", "b", False)).grid(
                row=2, column=i, sticky="w")
        self._row(f, 3, "Cyclic order n", ttk.Spinbox(f, from_=1, to=24, textvariable=self._var("cyclic", "i", 1), width=6))
        self._row(f, 4, "Rotation axis (part frame)", ttk.Combobox(f, textvariable=self._var("cyclic_axis", "s", "z"),
                                                                   values=["x", "y", "z"], width=4, state="readonly"))
        ttk.Checkbutton(f, text="Dihedral (add mirror planes through the axis)",
                        variable=self._var("dihedral", "b", False)).grid(row=5, column=0, columnspan=3, sticky="w")
        ttk.Button(f, text="Detect symmetry now", command=self.detect).grid(row=6, column=0, columnspan=3, pady=8, sticky="ew")
        ttk.Label(f, text="The part frame is the principal frame of the surface,\n"
                          "snapped to the world axes when nearly aligned.", foreground="#555").grid(
            row=7, column=0, columnspan=3, sticky="w")
        return f

    def _tab_exact(self, nb):
        f = ttk.Frame(nb, padding=8)
        ttk.Label(f, text="Surfaces").grid(row=0, column=0, sticky="w")
        sf = ttk.Frame(f)
        sf.grid(row=0, column=1, sticky="w")
        for k in ("plane", "cylinder", "cone", "sphere", "torus"):
            ttk.Checkbutton(sf, text=k, variable=self._var(f"s_{k}", "b", True)).pack(anchor="w")
        ttk.Checkbutton(f, text="Spline extrusions (non-algebraic)", variable=self._var("allow_extrusion", "b", True)).grid(
            row=1, column=0, columnspan=2, sticky="w")
        self._row(f, 2, "Snap directions within [deg]", ttk.Entry(f, textvariable=self._var("snap_directions_deg", "d", 1.0), width=8))
        self._row(f, 3, "Angle lattice step [deg]", ttk.Entry(f, textvariable=self._var("snap_angle_step_deg", "d", 15.0), width=8))
        ttk.Checkbutton(f, text="Enforce tangency (fillets, slots)", variable=self._var("tangency", "b", True)).grid(
            row=4, column=0, columnspan=2, sticky="w")
        self._row(f, 5, "Round dimensions to grid [mm] (0 = off)", ttk.Entry(f, textvariable=self._var("grid", "d", 0.0), width=8))
        ttk.Checkbutton(f, text="Rebuild blends as fillets on a sharp solid", variable=self._var("fillets", "b", True)).grid(
            row=6, column=0, columnspan=2, sticky="w")
        ttk.Checkbutton(f, text="Write a CadQuery script (verified)", variable=self._var("export_script", "b", True)).grid(
            row=7, column=0, columnspan=2, sticky="w")
        ttk.Checkbutton(f, text="Try stacked prisms (2.5D) first", variable=self._var("prism25", "b", True)).grid(
            row=8, column=0, columnspan=2, sticky="w")
        ttk.Checkbutton(f, text="Write a SolidWorks feature-tree script (2.5D parts)",
                        variable=self._var("export_solidworks", "b", True)).grid(row=9, column=0, columnspan=2, sticky="w")
        return f

    def _tab_smooth(self, nb):
        f = ttk.Frame(nb, padding=8)
        self._row(f, 0, "Voxel pitch [mm] (0 = auto)", ttk.Entry(f, textvariable=self._var("pitch", "d", 0.0), width=8))
        self._row(f, 1, "Smoothing blur [mm] (0 = auto)", ttk.Entry(f, textvariable=self._var("blur", "d", 0.0), width=8))
        self._row(f, 2, "Prune spurs shorter than [x radius]", ttk.Entry(f, textvariable=self._var("prune", "d", 1.0), width=8))
        self._row(f, 3, "Keep members found in >= [fraction of sectors]",
                  ttk.Entry(f, textvariable=self._var("min_copy_frac", "d", 0.5), width=8))
        ttk.Checkbutton(f, text="Detect bolt holes and build pads", variable=self._var("holes", "b", True)).grid(
            row=4, column=0, columnspan=2, sticky="w")
        ttk.Checkbutton(f, text="Keep free-ended members (horns, tabs)", variable=self._var("keep_free_ends", "b", True)).grid(
            row=5, column=0, columnspan=2, sticky="w")
        ttk.Checkbutton(f, text="Size members by cross-section area", variable=self._var("area_radius", "b", True)).grid(
            row=6, column=0, columnspan=2, sticky="w")
        self._row(f, 7, "Trim to design domain", ttk.Combobox(f, textvariable=self._var("clip", "s", "thin"),
                                                              values=["none", "thin", "all"], width=6, state="readonly"))
        return f

    # ---- actions -----------------------------------------------------------
    def options(self):
        g = lambda k: self.v[k].get()
        o = Options()
        for k in ("units", "mode", "symmetry", "mirror_x", "mirror_y", "mirror_z", "cyclic_axis", "dihedral",
                  "algebraic_only", "allow_extrusion", "tangency", "fillets", "export_script", "prism25", "export_solidworks", "holes", "keep_free_ends", "clip", "area_radius",
                  "out_dir"):
            setattr(o, k, g(k))
        for k in ("tolerance", "snap_directions_deg", "snap_angle_step_deg", "grid", "pitch", "blur", "prune",
                  "min_copy_frac"):
            setattr(o, k, float(g(k)))
        o.cyclic = int(g("cyclic"))
        o.surfaces = [k for k in ("plane", "cylinder", "cone", "sphere", "torus") if g(f"s_{k}")]
        o.out_dir = ""          # saving is explicit in the GUI
        return o

    def log(self, msg):
        self.q.put(("log", msg))

    def open_file(self):
        p = filedialog.askopenfilename(filetypes=[("STL", "*.stl *.STL"), ("All", "*.*")])
        if p:
            self.path = p
            self.mesh = None
            self.result = None
            self.file_label.config(text=os.path.basename(p))
            self.btn_save.config(state="disabled")
            self.analyze()

    def pick_dir(self):
        d = filedialog.askdirectory()
        if d:
            self.v["out_dir"].set(d)

    def _work(self, fn):
        if self.busy:
            return
        self.busy = True
        self.progress.start(12)
        for b in (self.btn_an, self.btn_run):
            b.config(state="disabled")

        def target():
            try:
                fn()
            except Exception as e:
                self.q.put(("log", "ERROR: " + str(e) + "\n" + traceback.format_exc()))
                self.q.put(("error", str(e)))
            self.q.put(("done", None))
        threading.Thread(target=target, daemon=True).start()

    def _load(self):
        from .core import load_mesh
        if self.mesh is None or getattr(self, "_mesh_units", None) != self.v["units"].get():
            self.mesh = load_mesh(self.path, self.v["units"].get())
            self._mesh_units = self.v["units"].get()
        return self.mesh

    def analyze(self):
        if not self.path:
            return messagebox.showinfo("stl2cad", "Open an STL file first.")

        def job():
            m = self._load()
            cf_ = coplanar_fraction(m)
            self.log(f"{os.path.basename(self.path)}: {len(m.faces)} triangles, extents "
                     f"{np.round(m.extents, 3).tolist()} mm, watertight={m.is_watertight}, "
                     f"volume {m.volume:.2f} mm^3")
            self.log(f"coplanar-area fraction {cf_:.2f}: looks "
                     f"{'CAD-like (exact mode)' if cf_ > 0.3 else 'organic (smooth mode)'}")
            self.q.put(("preview", (m, None, None)))
        self._work(job)

    def detect(self):
        if not self.path:
            return
        from .core import default_tolerance
        from .symmetry import detect

        def job():
            m = self._load()
            tol = self.v["tolerance"].get() or default_tolerance(m)
            spec, rep = detect(m, tol, log=lambda s: None)
            self.log(rep)
            self.q.put(("symmetry", spec))
        self._work(job)

    def run(self):
        if not self.path:
            return messagebox.showinfo("stl2cad", "Open an STL file first.")
        o = self.options()

        def job():
            m = self._load()
            res = run(self.path, o, log=self.log, mesh=m)
            self.result = res
            self.q.put(("preview", (m, res.mesh, res)))
        self._work(job)

    def save(self):
        if self.result is None:
            return
        out = self.v["out_dir"].get() or filedialog.askdirectory()
        if not out:
            return
        try:
            files = save(self.result, self.path, self.options(), log=self.log, out_dir=out)
            messagebox.showinfo("stl2cad", "Saved:\n" + "\n".join(files))
        except Exception as e:
            messagebox.showerror("stl2cad", str(e))

    # ---- UI updates --------------------------------------------------------
    def _poll(self):
        try:
            while True:
                kind, data = self.q.get_nowait()
                if kind == "log":
                    self.log_box.insert("end", str(data) + "\n")
                    self.log_box.see("end")
                elif kind == "preview":
                    self.draw(*data)
                elif kind == "symmetry":
                    spec = data
                    self.v["mirror_x"].set(0 in spec.mirrors)
                    self.v["mirror_y"].set(1 in spec.mirrors)
                    self.v["mirror_z"].set(2 in spec.mirrors)
                    self.v["cyclic"].set(spec.cyclic)
                    self.v["cyclic_axis"].set("xyz"[spec.cyclic_axis])
                    self.v["dihedral"].set(spec.dihedral)
                elif kind == "error":
                    messagebox.showerror("stl2cad", data)
                elif kind == "done":
                    self.busy = False
                    self.progress.stop()
                    for b in (self.btn_an, self.btn_run):
                        b.config(state="normal")
                    if self.result is not None:
                        self.btn_save.config(state="normal")
        except queue.Empty:
            pass
        self.after(100, self._poll)

    def draw(self, mesh, model=None, res=None):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        import matplotlib.cm as cm
        self.fig.clf()
        panels = [("original", mesh, None)]
        if model is not None:
            from .core import vertex_deviation
            d = vertex_deviation(mesh, model)
            panels.append(("result (deviation, mm)", model, d))
        rng = np.random.default_rng(0)
        lo, hi = mesh.bounds
        ctr, h = (lo + hi) / 2, (hi - lo).max() / 2
        light = np.array([0.4, -0.6, 0.8]) / np.linalg.norm([0.4, -0.6, 0.8])
        for k, (title, m, d) in enumerate(panels):
            ax = self.fig.add_subplot(1, len(panels), k + 1, projection="3d")
            if len(m.faces) > PREVIEW_FACES:
                try:
                    dm = m.simplify_quadric_decimation(face_count=PREVIEW_FACES)
                    if d is not None:
                        _, vi = __import__("scipy.spatial", fromlist=["cKDTree"]).cKDTree(m.vertices).query(dm.vertices)
                        d = d[vi]
                    m = dm
                except Exception:
                    pass
            idx = np.arange(len(m.faces))
            if len(idx) > PREVIEW_FACES:
                idx = np.sort(rng.choice(idx, PREVIEW_FACES, replace=False))
            sh = 0.35 + 0.65 * np.clip(np.abs(m.face_normals[idx] @ light), 0, 1)
            if d is None:
                cols = np.c_[0.62 * sh, 0.76 * sh, 0.93 * sh, np.ones_like(sh)]
            else:
                from matplotlib.colors import LinearSegmentedColormap
                cmap = LinearSegmentedColormap.from_list("dev", [(0.62, 0.76, 0.93), (0.98, 0.85, 0.3), (0.85, 0.1, 0.1)])
                fd = d[m.faces[idx]].mean(1)
                vmax = max(np.percentile(d, 99), 1e-6)
                cols = cmap(np.clip(fd / vmax, 0, 1))
                cols[:, :3] *= sh[:, None]
                sm = cm.ScalarMappable(cmap=cmap)
                sm.set_clim(0, vmax)
                self.fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
            ax.add_collection3d(Poly3DCollection(m.triangles[idx], facecolors=cols, edgecolor="none"))
            ax.set_xlim(ctr[0] - h, ctr[0] + h)
            ax.set_ylim(ctr[1] - h, ctr[1] + h)
            ax.set_zlim(ctr[2] - h, ctr[2] + h)
            ax.set_box_aspect((1, 1, 1))
            ax.view_init(28, -58)
            ax.set_axis_off()
            ax.set_title(title, fontsize=10)
        if res is not None:
            dv = res.deviation["ref_to_model"]
            self.fig.suptitle(f"{res.mode} mode: mean {dv['mean']:.3f} mm, p95 {dv['p95']:.3f} mm, "
                              f"volume {100 * (res.deviation['volume_model'] / res.deviation['volume_ref'] - 1):+.2f} %, "
                              f"valid={res.shape.isValid()}", fontsize=10)
        self.fig.tight_layout()
        self.canvas.draw_idle()


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except Exception:
        pass
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
