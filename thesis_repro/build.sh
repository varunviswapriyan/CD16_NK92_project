#!/usr/bin/env bash
# build.sh -- compile Giuseppe's SPPARKS (the with-CD4s thesis model) with MPI.
#
# Run from the repo root, on the branch Sex_bias_post-Zihai_with_CD4s.
# The committed CMakeLists.txt only lists 4 files and won't link; SPPARKS is
# built by compiling every source at the repo root with an MPI C++ compiler.
set -euo pipefail
cd "$(dirname "$0")/.."     # -> repo root (assumes this kit lives in repo/thesis_repro)

echo "== branch: $(git rev-parse --abbrev-ref HEAD)  commit: $(git rev-parse --short HEAD)"

# --- pick an MPI compiler -------------------------------------------------
# On OSC you typically need a module first, e.g.:
#   module load gcc openmpi          (or: module spider openmpi  to see names)
# Giuseppe historically used GCC + OpenMPI. If mpicxx isn't found, load it.
if ! command -v mpicxx >/dev/null 2>&1; then
  echo "!! mpicxx not found. Load your MPI module first, e.g.:"
  echo "     module spider openmpi      # find the exact name"
  echo "     module load gcc openmpi"
  exit 1
fi
echo "== using $(mpicxx --version | head -1)"

# --- compile -------------------------------------------------------------
# Old SPPARKS source + modern compiler => allow the legacy constructs.
mpicxx -O3 -std=c++14 -w -fpermissive *.cpp -o spk_sexbias
echo "== built ./spk_sexbias"
ls -l spk_sexbias
