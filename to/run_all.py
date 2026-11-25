"""TO regression, same as run_all.sh, for any platform (Windows PowerShell included).

Every pair in ../models/to at every detail level; the SolidWorks scripts are copied to
../results/sw_scripts/to/ as NAME_dL_solidworks.py (L = detail level, 1: fewest dimensions ...
5: smallest model error); the summary goes to summary.txt.

Pairs in SWEEP are also built as swept members (NAME_sweep_dL_solidworks.py, swsweep.py).

Usage (from this folder):  python run_all.py [--levels 1 2 3 4 5] [--jobs 4]
"""
import argparse, os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.join(HERE, '..', 'models', 'to')
OUT = os.path.join(HERE, '..', 'results', 'sw_scripts', 'to')
PAIRS = [  # original, TO, name, units
    ('LBracket.STL', 'LBracket_TO_1.stl', 'lbracket_to', 'm'),
    ('EdgeCantilever.STL', 'EdgeCantilever_TO_1.stl', 'cantilever1_to', 'm'),
    ('EdgeCantilever.STL', 'EdgeCantilever_TO_2.stl', 'cantilever2_to', 'm'),
    ('EdgeCantilever.STL', 'EdgeCantilever_TO_3.stl', 'cantilever3_to', 'm'),
    ('CircularPlateHole.STL', 'CircularPlateHole_TO_1.stl', 'circplate_to', 'in'),
    ('CrankArm.STL', 'CrankArm_TO_1.stl', 'crankarm_to', 'in'),
]
PAIRS = [pr for pr in PAIRS if os.path.exists(os.path.join(M, pr[0])) and os.path.exists(os.path.join(M, pr[1]))]
# pairs also built as swept members (swsweep.py: hubs, curved members with tapered ellipse sections, joint nodes);
# the method needs hubs (holes of the original that the TO rings)
SWEEP = {'crankarm_to'}


def run_sweep(o, t, name, units, level):
    base = os.path.join('out_' + os.path.splitext(o)[0], os.path.splitext(o)[0] + '_exact.py')
    with open(name + '.log', 'w') as f:
        subprocess.run([sys.executable, 'swsweep.py', base, os.path.join(M, t), name, '--units', units,
                        '--detail', str(level)], stdout=f, stderr=subprocess.STDOUT, cwd=HERE)


def run(o, t, name, units, extra):
    with open(name + '.log', 'w') as f:
        subprocess.run([sys.executable, 'run_plates.py', os.path.join(M, o), os.path.join(M, t), name,
                        '--units', units] + extra, stdout=f, stderr=subprocess.STDOUT, cwd=HERE)


def grab(path, pattern):
    try:
        m = re.search(pattern, open(path).read())
        return m.group(1) if m else ''
    except OSError:
        return ''


ap = argparse.ArgumentParser()
ap.add_argument('--levels', type=int, nargs='+', default=[1, 2, 3, 4, 5])
ap.add_argument('--jobs', type=int, default=4)
a = ap.parse_args()
os.chdir(HERE)
os.makedirs(OUT, exist_ok=True)

# exact models of the originals first, one run each (parallel runs would build them at the same time)
done = set()
for o, t, n, u in PAIRS:
    if o not in done:
        run(o, t, '_base', u, ['--fit-only'])
        done.add(o)

jobs = [(o, t, '%s_d%d' % (n, L), u, ['--sw', '--detail', str(L), '--plot']) for L in a.levels for o, t, n, u in PAIRS]
sweeps = [(o, t, '%s_sweep_d%d' % (n.replace('_to', ''), L), u, L) for L in a.levels for o, t, n, u in PAIRS
          if n in SWEEP]
with ThreadPoolExecutor(a.jobs) as ex:
    list(ex.map(lambda j: run(*j), jobs))
    list(ex.map(lambda j: run_sweep(*j), sweeps))

lines = ['%-18s %-6s %-10s %-24s %s' % ('pair', 'level', 'stand-in', 'volume vs TO', 'independent dimensions')]
for o, t, n, u in PAIRS:
    for L in a.levels:
        nm = '%s_d%d' % (n, L)
        if os.path.exists(nm + '_solidworks.py'):
            shutil.copy(nm + '_solidworks.py', OUT)
        lines.append('%-18s %-6d %-10s %-24s %s' % (
            n, L, grab(nm + '.log', r'API stand-in: (\w+)') or 'FAILED', grab(nm + '.log', r'(\([-+][0-9.]* %\))'),
            grab(nm + '_standin.log', r'independent dimensions: (\d+)')))
for o, t, n, u, L in sweeps:
    if os.path.exists(n + '_solidworks.py'):
        shutil.copy(n + '_solidworks.py', OUT)
    lines.append('%-18s %-6d %-10s %-24s %s' % (
        n.rsplit('_d', 1)[0], L, grab(n + '.log', r'API stand-in: (\w+)') or 'FAILED',
        grab(n + '.log', r'volume: TO.*(\([-+][0-9.]* %\))'), grab(n + '_standin.log', r'independent dimensions: (\d+)')))
print('\n'.join(lines))
open('summary.txt', 'w').write('\n'.join(lines) + '\n')
