#!/usr/bin/env bash
# Read-only status of every fit. Runs no jobs and writes nothing to the output dirs.
set -u

CA_DIR="${CA_DIR:-/home/gddaslab/vxv016/Ca_fit_c02}"
PZ_DIR="${PZ_DIR:-/home/gddaslab/vxv016/boot_pzap}"
SCRIPTS="${SCRIPTS:-$HOME/CD16_NK92_project/filesCC}"
PY="${PY:-python}"

hr(){ printf '\n%s\n%s\n%s\n' "======================================================================" "$1" "======================================================================"; }

hr "QUEUE"
squeue -u "$USER"

hr "Ca bootstrap -- sample counts per variant"
if [ -d "$CA_DIR/out_boot_ca2" ]; then
  n=0
  for d in "$CA_DIR"/out_boot_ca2/*/; do
    [ -d "$d" ] || continue
    c=$(ls "$d"/boot_*.json 2>/dev/null | wc -l)
    printf '  %-12s %s\n' "$(basename "$d")" "$c"
    n=$((n+c))
  done
  [ "$n" -eq 0 ] && echo "  -> no completed samples anywhere; the Ca bootstrap has not run"
else
  echo "  !! $CA_DIR/out_boot_ca2 does not exist"
fi

hr "Ca bootstrap -- report  (must run from $CA_DIR)"
( cd "$CA_DIR" && $PY "$SCRIPTS/bootstrap_ca2.py" report )

hr "pZAP bootstrap -- report"
$PY "$SCRIPTS/bootstrap_pzap.py" report

hr "pZAP bootstrap -- sample counts"
for k in emp par; do
  c=$(ls "$PZ_DIR"/$k*/boot_*.json 2>/dev/null | wc -l)
  printf '  %-5s %s\n' "$k" "$c"
done

hr "POINT 2: ZAP0 x kzp grid"
bash "$SCRIPTS/zap_kzp_map.sh" report

hr "DONE"
echo "Plots:  python plot_ca_fit.py orig   |   python plot_pzap_fit.py"
