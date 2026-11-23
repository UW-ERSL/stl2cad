#!/usr/bin/env bash
# Regression run over the meshes in exact-models/. Usage: regression/run_all.sh [outdir] [jobs]
cd "$(dirname "$0")/.."
OUT=${1:-/tmp/stl2cad_reg}
JOBS=${2:-4}
mkdir -p "$OUT"

# file  units
cat > "$OUT/parts.txt" <<'EOF'
multipleFeatures.stl in
angle_block.stl in
idler_riser.stl in
tray.stl mm
plate_holes.stl mm
Gear.STL mm
AlcoaGrabCAD.stl in
BicycleCrank.STL m
CrankArm.STL in
threeholebracket.STL m
SwingArm.STL m
SimpleBracket.STL m
LBracket.STL m
Knuckle.STL m
GEGrabCAD.STL in
CircularPlateHole.STL in
HeatBlock.STL m
C-bracket.STL m
StepBlock.STL m
slottedBlock.STL m
bottleopener.STL m
lid.STL m
hexspanner.STL m
hexcap.STL m
lever.STL m
rack.STL m
heatfin.STL m
crossbar.STL m
flange.STL m
standoff.STL m
chesspawn.STL m
topcover.STL m
snapring.STL m
EOF

run_one() {
  f=$1; u=$2; out=$3
  stem=${f%.*}
  start=$(date +%s)
  timeout 900 python -m stl2cad.cli "exact-models/$f" --units "$u" -o "$out/$stem" > "$out/$stem.log" 2>&1
  rc=$?
  echo "$stem rc=$rc t=$(( $(date +%s) - start ))s"
}
export -f run_one
xargs -P "$JOBS" -L 1 bash -c 'run_one "$0" "$1" '"$OUT" < "$OUT/parts.txt"
