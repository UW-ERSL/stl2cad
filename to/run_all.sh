#!/usr/bin/env bash
# TO regression: every pair in ../models/to, SolidWorks scripts copied to ../results/sw_scripts/to/.
# Usage (from this folder): ./run_all.sh [detail 1-5]
cd "$(dirname "$0")"
D=${1:-4}
M=../models/to
OUT=../results/sw_scripts/to
mkdir -p "$OUT"
# original  TO  name  units
while read -r O T N U; do
  python run_plates.py "$M/$O" "$M/$T" "$N" --units "$U" --sw --detail "$D" --plot > "$N.log" 2>&1
  echo "$N: $(grep -o 'API stand-in: [a-zA-Z]*' "$N.log") $(grep -o 'volume: TO .*' "$N.log") $(grep -o 'independent dimensions: [0-9]*' "${N}_standin.log")"
  cp "${N}_solidworks.py" "$OUT/"
done <<'PAIRS'
LBracket.STL LBracket_TO_1.stl lbracket_to m
EdgeCantilever.STL EdgeCantilever_TO_1.stl cantilever1_to m
EdgeCantilever.STL EdgeCantilever_TO_2.stl cantilever2_to m
EdgeCantilever.STL EdgeCantilever_TO_3.stl cantilever3_to m
PAIRS
