#!/usr/bin/env bash
# build_all.sh -- one-shot build of Giuseppe's with-CD4s thesis SPPARKS on the cluster.
# Run it from anywhere on the cluster:   bash build_all.sh
#
# It: (1) clones the public SPPARKS for the stock files Giuseppe didn't commit,
#     (2) copies in only the missing core files (no foreign apps),
#     (3) regenerates the style_*.h headers to match the files present,
#     (4) compiles with mpicxx -> ./spk_sexbias
set -uo pipefail

REPO=~/CD16_NK92_project/SPPARKS_Sex_Bias
cd "$REPO" || { echo "!! $REPO not found"; exit 1; }
echo "== working in $REPO"

# 1. public SPPARKS (for the unchanged stock files)
if [ ! -d ~/spparks_official ]; then
  echo "== cloning public SPPARKS (one time)"
  git clone https://github.com/spparks/spparks.git ~/spparks_official || {
    echo "!! could not clone SPPARKS"; exit 1; }
fi

# 2. copy ONLY the stock core files Giuseppe left out (never overwrites his files)
echo "== filling missing core files"
for b in comm_lattice comm_off_lattice random_park random_mars potential irregular; do
  cp -n ~/spparks_official/src/$b.cpp . 2>/dev/null
  cp -n ~/spparks_official/src/$b.h  . 2>/dev/null
done

# 3. regenerate style headers so they list only styles whose files are present
echo "== regenerating style_*.h"
mk(){ grep -l "$2(" $3 2>/dev/null | sort | sed 's|^|#include "|; s|$|"|' > "$1"; }
mk style_app.h     AppStyle     'app_*.h'
mk style_command.h CommandStyle '*.h'
mk style_diag.h    DiagStyle    'diag_*.h'
mk style_dump.h    DumpStyle    'dump_*.h'
mk style_pair.h    PairStyle    'pair_*.h'
mk style_region.h  RegionStyle  'region_*.h'
mk style_solve.h   SolveStyle   'solve_*.h'

# 4. compile
echo "== loading MPI + compiling"
module load foss/2023a 2>/dev/null || echo "!! could not module load foss/2023a -- load your MPI module, then rerun"
command -v mpicxx >/dev/null || { echo "!! mpicxx not on PATH"; exit 1; }
if mpicxx -O3 -std=c++14 -w -fpermissive *.cpp -o spk_sexbias; then
  echo "== BUILD_OK"; ls -l spk_sexbias
else
  echo "== BUILD FAILED -- if it's a 'No such file: X.h', add X to the loop in step 2 and rerun"
fi
