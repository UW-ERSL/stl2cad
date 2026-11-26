"""TO regression, same as run_all.sh, for any platform (Windows PowerShell included).

Every pair in ../models/to at every detail level; the SolidWorks scripts are copied to
../models/to/sw_scripts/ as NAME_dL_solidworks.py (L = detail level, 1: fewest dimensions ...
5: smallest model error); the summary goes to summary.txt.

Pairs in SWEEP are also built as swept members (NAME_sweep_dL_solidworks.py, swsweep.py).
Pairs in FRAME are built as 3D frames only (NAME_dL_solidworks.py, swframe.py): hubs on several axes, straight
round members and sphere nodes in 3D (the Knuckle: a clevis whose TO joins a boss and two eyes on other axes).

Usage (from this folder):  python run_all.py [--levels 1 2 3 4 5] [--jobs 4]
"""
import argparse, os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.join(HERE, '..', 'models', 'to')
OUT = os.path.join(HERE, '..', 'models', 'to', 'sw_scripts')
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
# pairs built as 3D frames (swframe.py); their originals' exact models come from stl2cad.cli directly
FRAME = [  # original, TO, name, units
    ('Knuckle.STL', 'Knuckle_TO_1.stl', 'knuckle1_frame', 'm'),
]
# Knuckle_TO_2.stl is a thick hollowed wall, not a truss: as a frame it becomes a chain of nodes (not listed)
FRAME = [pr for pr in FRAME if os.path.exists(os.path.join(M, pr[0])) and os.path.exists(os.path.join(M, pr[1]))]


def base_of(o):
    stem = os.path.splitext(o)[0]
    return os.path.join('out_' + stem, stem + '_exact.py')


def run_frame(o, t, name, units, level):
    with open(name + '.log', 'w') as f:
        subprocess.run([sys.executable, 'swframe.py', base_of(o), os.path.join(M, t), name, '--units', units,
                        '--detail', str(level), '--plot'], stdout=f, stderr=subprocess.STDOUT, cwd=HERE)


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

for o, t, n, u in FRAME:
    if not os.path.exists(base_of(o)):
        subprocess.run([sys.executable, '-m', 'stl2cad.cli', os.path.join(M, o), '--units', u, '--no-solidworks',
                        '-o', 'out_' + os.path.splitext(o)[0]], stdout=subprocess.DEVNULL, cwd=HERE)

jobs = [(o, t, '%s_d%d' % (n, L), u, ['--sw', '--detail', str(L), '--plot']) for L in a.levels for o, t, n, u in PAIRS]
frames = [(o, t, '%s_d%d' % (n, L), u, L) for L in a.levels for o, t, n, u in FRAME]
sweeps = [(o, t, '%s_sweep_d%d' % (n.replace('_to', ''), L), u, L) for L in a.levels for o, t, n, u in PAIRS
          if n in SWEEP]
with ThreadPoolExecutor(a.jobs) as ex:
    list(ex.map(lambda j: run(*j), jobs))
    list(ex.map(lambda j: run_sweep(*j), sweeps))
    list(ex.map(lambda j: run_frame(*j), frames))

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
for o, t, n, u, L in frames:
    if os.path.exists(n + '_solidworks.py'):
        shutil.copy(n + '_solidworks.py', OUT)
    lines.append('%-18s %-6d %-10s %-24s %s' % (
        n.rsplit('_d', 1)[0], L, grab(n + '.log', r'API stand-in: (\w+)') or 'FAILED',
        grab(n + '.log', r'volume: TO.*(\([-+][0-9.]* %\))'), grab(n + '_standin.log', r'independent dimensions: (\d+)')))
print('\n'.join(lines))
open('summary.txt', 'w').write('\n'.join(lines) + '\n')
