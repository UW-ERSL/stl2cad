#!/usr/bin/env bash
# TO regression: every pair in ../models/to at every detail level. The SolidWorks scripts are copied to
# ../results/sw_scripts/to/ as NAME_dL_solidworks.py (L = detail level, 1: fewest dimensions ... 5: smallest
# model error); the summary goes to summary.txt.
# Usage (from this folder): ./run_all.sh [levels] [jobs]     e.g. ./run_all.sh "1 4" 2   (default "1 2 3 4 5" 4)
cd "$(dirname "$0")"
LEVELS=${1:-1 2 3 4 5}
JOBS=${2:-4}
M=../models/to
OUT=../results/sw_scripts/to
mkdir -p "$OUT"

# original  TO  name  units
PAIRS='LBracket.STL LBracket_TO_1.stl lbracket_to m
EdgeCantilever.STL EdgeCantilever_TO_1.stl cantilever1_to m
EdgeCantilever.STL EdgeCantilever_TO_2.stl cantilever2_to m
EdgeCantilever.STL EdgeCantilever_TO_3.stl cantilever3_to m'

# exact models of the originals first, one run each (parallel runs would build them at the same time)
echo "$PAIRS" | awk '!seen[$1]++' | while read -r O T N U; do
  python run_plates.py "$M/$O" "$M/$T" "_base" --units "$U" --fit-only > /dev/null 2>&1
done

run() {   # original TO name units level
  local n="$3_d$5"
  python run_plates.py "$M/$1" "$M/$2" "$n" --units "$4" --sw --detail "$5" --plot > "$n.log" 2>&1
  cp "${n}_solidworks.py" "$OUT/" 2>/dev/null
}
export -f run; export M OUT

for L in $LEVELS; do echo "$PAIRS" | sed "s/\$/ $L/"; done |
  xargs -P "$JOBS" -L 1 bash -c 'run "$@"' _

{
  printf '%-18s %-6s %-10s %-24s %s\n' pair level stand-in "volume vs TO" "independent dimensions"
  echo "$PAIRS" | while read -r O T N U; do
    for L in $LEVELS; do
      n="${N}_d$L"
      ok=$(grep -o 'API stand-in: [a-zA-Z]*' "$n.log" | awk '{print $3}')
      vol=$(grep -o '([-+][0-9.]* %)' "$n.log")
      dims=$(grep -o 'independent dimensions: [0-9]*' "${n}_standin.log" 2>/dev/null | awk '{print $3}')
      printf '%-18s %-6s %-10s %-24s %s\n' "$N" "$L" "${ok:-FAILED}" "$vol" "$dims"
    done
  done
} | tee summary.txt
