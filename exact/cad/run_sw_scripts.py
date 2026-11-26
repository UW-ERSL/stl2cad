"""Run every SolidWorks Python script in a folder, one after the other.

Usage:  python run_sw_scripts.py [FOLDER] [--pattern "*_solidworks.py"] [--close]

FOLDER defaults to this script's folder. Scripts run in name order; each one's output
goes to FOLDER/logs/NAME.log, and a pass/fail summary is printed at the end.
--close closes all open SolidWorks documents (without saving) after each script.
"""
import argparse, glob, os, subprocess, sys, time

ap = argparse.ArgumentParser()
ap.add_argument('folder', nargs='?', default=os.path.dirname(os.path.abspath(__file__)))
ap.add_argument('--pattern', default='*.py')
ap.add_argument('--close', action='store_true')
a = ap.parse_args()

# scripts to run, excluding this runner
me = os.path.abspath(__file__)
scripts = sorted(p for p in glob.glob(os.path.join(a.folder, a.pattern)) if os.path.abspath(p) != me)
logs = os.path.join(a.folder, 'logs')
os.makedirs(logs, exist_ok=True)

if a.close:
    import win32com.client
    sw = win32com.client.Dispatch('SldWorks.Application')

# run each script in its own process
results = []
for i, p in enumerate(scripts, 1):
    name = os.path.basename(p)
    print('[%d/%d] %s ...' % (i, len(scripts), name), end=' ', flush=True)
    t0 = time.time()
    with open(os.path.join(logs, os.path.splitext(name)[0] + '.log'), 'w') as f:
        rc = subprocess.run([sys.executable, p], stdout=f, stderr=subprocess.STDOUT, cwd=a.folder).returncode
    dt = time.time() - t0
    print('ok' if rc == 0 else 'FAILED (exit %d)' % rc, '%.1f s' % dt)
    results.append((name, rc, dt))
    if a.close:
        sw.CloseAllDocuments(True)

# summary
bad = [r for r in results if r[1] != 0]
print('\n%d scripts, %d ok, %d failed' % (len(results), len(results) - len(bad), len(bad)))
for name, rc, dt in bad:
    print('  FAILED  %s  (see logs/%s.log)' % (name, os.path.splitext(name)[0]))
