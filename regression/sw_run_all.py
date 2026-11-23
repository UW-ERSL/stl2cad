"""Run every SolidWorks script in a folder against real SolidWorks and summarize the results.

Usage (Windows, from the repository root, in the stl2cad environment with pywin32):
    python regression/sw_run_all.py                     # run sw_scripts/*_solidworks.py
    python regression/sw_run_all.py --generate          # first regenerate the scripts from exact-models/
    python regression/sw_run_all.py --only rack lid     # some parts only
    python regression/sw_run_all.py --save parts        # also save each part as parts/<name>.SLDPRT

--generate runs stl2cad on every part listed in regression/run_all.sh (with the units given there)
and copies each <part>_solidworks.py it writes into the scripts folder. Each script runs in its own
Python process with a time limit; after it, all SolidWorks documents are closed (unless
--keep-open). The log of every run goes to <logs>/<name>.log; the summary is printed and written to
<logs>/summary.txt.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ----------------------------------------------------------------------------- parts list
def parts_list():
    """(file, units) of the regression parts, read from run_all.sh."""
    text = open(os.path.join(ROOT, "regression", "run_all.sh")).read()
    block = re.search(r"<<'EOF'\n(.*?)\nEOF", text, re.S).group(1)
    return [tuple(line.split()) for line in block.splitlines() if line.strip()]


def generate(scripts, work, only, timeout):
    """Regenerate the SolidWorks scripts with the current stl2cad."""
    os.makedirs(scripts, exist_ok=True)
    for f, units in parts_list():
        stem = os.path.splitext(f)[0]
        if only and stem.lower() not in only:
            continue
        out = os.path.join(work, stem)
        t0 = time.time()
        try:
            r = subprocess.run([sys.executable, "-m", "stl2cad.cli", os.path.join("exact-models", f),
                                "--units", units, "-o", out], cwd=ROOT, capture_output=True, text=True,
                               timeout=timeout)
            rc = r.returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
        src = os.path.join(out, stem + "_solidworks.py")
        if os.path.exists(src):
            shutil.copy(src, os.path.join(scripts, stem + "_solidworks.py"))
            what = "script written"
        else:
            what = "no SolidWorks script"
        print(f"  generate {stem:20s} rc={rc}  {what}  ({time.time() - t0:.0f}s)", flush=True)


# ----------------------------------------------------------------------------- running
def close_all():
    """Close every SolidWorks document without saving (keeps the session light)."""
    try:
        import win32com.client
        app = win32com.client.Dispatch("SldWorks.Application")
        app.CloseAllDocuments(True)
    except Exception as ex:
        print(f"  (could not close the SolidWorks documents: {type(ex).__name__})")


def summarize(name, rc, text, seconds):
    """One result row from a script's output."""
    row = dict(name=name, seconds=seconds, ok=rc == 0 and "Traceback" not in text, error="")
    if not row["ok"]:
        tail = [ln for ln in text.strip().splitlines() if ln.strip()]
        row["error"] = "timeout" if rc == "timeout" else (tail[-1][:150] if tail else f"exit code {rc}")
    m = re.search(r"dimensions: (\d+) added(?:, (\d+) failed)?", text)
    row["dims"] = (int(m.group(1)), int(m.group(2) or 0)) if m else None
    m = re.search(r"relations: (\d+) added(?:, (\d+) failed)?", text)
    row["rels"] = (int(m.group(1)), int(m.group(2) or 0)) if m else None
    m = re.search(r"sketch status: (.*)", text)
    st = [s.strip() for s in m.group(1).split(",")] if m else []
    row["sketches"] = (sum(s == "fully defined" for s in st), len(st))
    row["features"] = len(re.findall(r"^\S.*: volume [0-9.]+ mm\^3 \(expected", text, re.M))
    row["redrawn"] = len(re.findall(r"wrong shape; drawing the sketch again", text))
    row["notes"] = len(re.findall(r"^NOTE", text, re.M))
    return row


def fmt(row):
    if not row["ok"]:
        return f"{row['name']:22s} FAILED  {row['seconds']:5.0f}s  {row['error']}"
    d, r, (full, n) = row["dims"] or (0, 0), row["rels"] or (0, 0), row["sketches"]
    extra = []
    if d[1]:
        extra.append(f"{d[1]} dimension(s) failed")
    if r[1]:
        extra.append(f"{r[1]} relation(s) failed")
    if row["redrawn"]:
        extra.append(f"{row['redrawn']} sketch(es) redrawn without plan")
    if row["notes"]:
        extra.append(f"{row['notes']} note(s)")
    return (f"{row['name']:22s} ok      {row['seconds']:5.0f}s  {row['features']:2d} feature(s), "
            f"{d[0]} dims, {r[0]} relations, {full}/{n} sketches fully defined"
            + (f"  [{'; '.join(extra)}]" if extra else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scripts", default=os.path.join(ROOT, "sw_scripts"), help="folder of *_solidworks.py")
    ap.add_argument("--logs", default=os.path.join(ROOT, "sw_runs"), help="folder for logs and the summary")
    ap.add_argument("--generate", action="store_true", help="regenerate the scripts from exact-models/ first")
    ap.add_argument("--only", nargs="*", default=[], help="part names (case-insensitive)")
    ap.add_argument("--save", default="", help="folder to save each part as .SLDPRT")
    ap.add_argument("--timeout", type=int, default=900, help="seconds per part")
    ap.add_argument("--keep-open", action="store_true", help="do not close the SolidWorks documents")
    a = ap.parse_args()
    only = {s.lower() for s in a.only}
    os.makedirs(a.logs, exist_ok=True)
    if a.generate:
        print("regenerating the scripts with the current stl2cad:", flush=True)
        generate(a.scripts, os.path.join(a.logs, "generated"), only, a.timeout)
    scripts = sorted(p for p in glob.glob(os.path.join(a.scripts, "*_solidworks.py")))
    scripts = [p for p in scripts if not only or os.path.basename(p)[:-len("_solidworks.py")].lower() in only]
    if not scripts:
        sys.exit(f"no *_solidworks.py in {a.scripts}")
    if a.save:
        os.makedirs(a.save, exist_ok=True)
    print(f"running {len(scripts)} script(s) from {a.scripts}", flush=True)
    rows = []
    for p in scripts:
        name = os.path.basename(p)[:-len("_solidworks.py")]
        cmd = [sys.executable, "-u", p]          # unbuffered: a log survives a timeout
        if a.save:
            cmd.append(os.path.abspath(os.path.join(a.save, name + ".SLDPRT")))
        t0 = time.time()
        log = os.path.join(a.logs, name + ".log")
        with open(log, "w") as fh:               # written as the script runs
            proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=os.path.dirname(p))
            try:
                rc = proc.wait(timeout=a.timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
                rc = "timeout"
        text = open(log, errors="replace").read()
        row = summarize(name, rc, text, time.time() - t0)
        rows.append(row)
        print(fmt(row), flush=True)
        if not a.keep_open:
            close_all()
    n_ok = sum(r["ok"] for r in rows)
    full = sum(r["ok"] and r["sketches"][1] and r["sketches"][0] == r["sketches"][1] for r in rows)
    lines = [fmt(r) for r in rows] + ["", f"{n_ok} of {len(rows)} part(s) built; {full} with every sketch fully defined"]
    with open(os.path.join(a.logs, "summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n" + lines[-1])
    print(f"logs and summary: {a.logs}")


if __name__ == "__main__":
    main()
