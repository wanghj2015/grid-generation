#!/bin/bash
#
# Run the analytical dipole coordinate generator.
#
#   bash scripts/run/run_dipole.sh

# Repository root, resolved relative to this script (no hardcoded home dir).
GRIDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$GRIDIR/build/bin/dipole_coordinates"
RUNDIR="$GRIDIR/run/dipole"

if [ ! -x "$EXE" ]; then
  echo "Executable not found: $EXE"
  echo "Build first:  cmake -S \"$GRIDIR\" -B \"$GRIDIR/build\" && cmake --build \"$GRIDIR/build\""
  exit 1
fi

mkdir -p "$RUNDIR"
cd "$RUNDIR" || exit 1

cp -f "$GRIDIR/data/igrf13coeffs.txt" ./

cat > grid_params.namelist << EOF
&grid_params
epoch   = 2000.0,
nptx    = 101,
nlp     = 45,
nmp     = 1,
theta1  = 45.0,
theta2  = 82.0,
hmin    = 90.0,
ds_step = 5.0,
use_psi = 0,
aa      = 3.0,
use_dipole  = 1,
/
EOF

"$EXE" > dipole.log 2>&1
echo "Done. Output written to $RUNDIR"
