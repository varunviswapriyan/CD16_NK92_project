#!/bin/bash
# ============================================================================
# setup_ca_batch.sh  -  build + (optionally) submit the Ca fitting batch c01-c10
#   bash setup_ca_batch.sh          -> build folders ~/Ca_fit_cNN (no submit)
#   bash setup_ca_batch.sh submit   -> build + sbatch all
# Needs CD16_v2 env active (exports pZAP curves for v57 if not already done).
# Each variant = a copy of Indrani's Ca folder + a pZAP input + ca_config.json + fit_ca.py.
# ============================================================================
set -e
cd ~
HOME_DIR=/home/gddaslab/vxv016
SRC=$(cd "$(dirname "$0")" && pwd)
CABASE=$HOME_DIR/Ca_fitting_NK92
SUB=estimate_params_pzap_cleaned_up
MAIL=Varun.Viswapriyan@nationwidechildrens.org

# pZAP inputs: exported curves (10 reps) from the pZAP fits, and Indrani's original
P58=$HOME_DIR/NK92_fit_v58/$SUB/model_output_pzap_v58.csv
P57=$HOME_DIR/NK92_fit_v57/$SUB/model_output_pzap_v57.csv
P0=$CABASE/optimized_model_pzap/model_output_pzap.csv
for V in 57 58; do
  F=$HOME_DIR/NK92_fit_v$V/$SUB/model_output_pzap_v$V.csv
  if [ ! -f "$F" ]; then
    echo "exporting pZAP curves for v$V (10 reps, ~5 min)..."
    cp "$SRC/export_pzap.py" $HOME_DIR/NK92_fit_v$V/$SUB/
    (cd $HOME_DIR/NK92_fit_v$V/$SUB && python export_pzap.py 10 model_output_pzap_v$V.csv)
  fi
done

# ------------------------------------------------------------------------------------------------
# id  | pzap | t0 | params                 | fixed              | bounds override        | opt   | weighted | description
# ------------------------------------------------------------------------------------------------
VARIANTS="
c01 | P58 | 30 | C1,C2,g,k3,k4          | {}                 | {}                     | lbfgs | n | v58 pZAP, WIDE bounds, her optimizer (48 starts) - isolates the bound fix
c02 | P58 | 30 | C1,C2,g,k3,k4          | {}                 | {}                     | pso   | n | v58 pZAP, wide bounds, PSO+polish (MAIN)
c03 | P57 | 30 | C1,C2,g,k3,k4          | {}                 | {}                     | pso   | n | v57 pZAP, same as c02
c04 | P0  | 30 | C1,C2,g,k3,k4          | {}                 | {}                     | pso   | n | Indrani pZAP, wide bounds, PSO (control: does widening fix HER fit too?)
c05 | P58 | 0  | C1,C2,g,k3,k4          | {}                 | {}                     | pso   | n | c02 fitted from t=0 instead of 30
c06 | P58 | 30 | C1,C2,g,k3             | {\"k4\":0.0}         | {}                     | pso   | n | pure Hill (k4=0, her earlier form)
c07 | P58 | 30 | C1,C2,g,k3,k4,nn       | {}                 | {}                     | pso   | n | Hill coefficient fitted (was fixed 9)
c08 | P58 | 30 | C1,C2,g,k3,k4          | {}                 | {}                     | pso   | y | SE-weighted SSR (chi-square)
c09 | P58 | 30 | C1,C2,g,k3,k4,nn,be    | {}                 | {}                     | pso   | n | everything: nn + baseline influx be fitted
c10 | P58 | 30 | C1,C2,g,k3,k4          | {}                 | {\"C1\":[1,4],\"C2\":[-2,1],\"g\":[-5,-1],\"k3\":[-3,1],\"k4\":[-3,1]} | pso | n | HER bounds with PSO (is it bounds or optimizer?)
"

[ -d "$CABASE" ] || { echo "Ca base folder not found: $CABASE"; exit 1; }
for f in fit_ca.py check_ca.py; do [ -f "$SRC/$f" ] || { echo "missing $SRC/$f"; exit 1; }; done

echo "$VARIANTS" | grep -E '^c[0-9]' | while IFS='|' read -r ID PZ T0 PARAMS FIXED BOUNDS OPT WT DESC; do
  ID=$(echo $ID); PZ=$(echo $PZ); T0=$(echo $T0); PARAMS=$(echo $PARAMS); FIXED=$(echo $FIXED); BOUNDS=$(echo $BOUNDS); OPT=$(echo $OPT); WT=$(echo $WT); DESC=$(echo $DESC)
  case $PZ in P58) PZF=$P58;; P57) PZF=$P57;; P0) PZF=$P0;; esac
  [ -f "$PZF" ] || { echo "pZAP file missing: $PZF"; exit 1; }
  DEST=$HOME_DIR/Ca_fit_$ID
  echo "=== building $ID : $DESC"
  rm -rf "$DEST"; mkdir -p "$DEST/optimized_model_pzap"
  cp -r $CABASE/ca_data "$DEST/"
  cp "$PZF" "$DEST/optimized_model_pzap/model_output_pzap.csv"
  cp "$SRC/fit_ca.py" "$SRC/check_ca.py" "$DEST/"
  [ "$WT" == "y" ] && WTJ=true || WTJ=false
  python - <<PYEOF
import json
cfg = {"note": "$ID: $DESC", "pzap": "optimized_model_pzap/model_output_pzap.csv", "ca": "ca_data/Ca_NK92.csv",
       "t0": $T0, "params": "$PARAMS".split(","), "fixed": json.loads('$FIXED'), "bounds": json.loads('$BOUNDS'),
       "optimizer": "$OPT", "particles": 48, "iters": 200, "starts": 48, "seed": 7, "weighted": $WTJ}
json.dump(cfg, open("$DEST/ca_config.json", "w"), indent=2)
print("    config:", cfg["params"], "fixed", cfg["fixed"], cfg["optimizer"], "weighted" if cfg["weighted"] else "")
PYEOF
  cat > "$DEST/run.sh" <<WRAP
#!/bin/bash
#SBATCH --mail-user=$MAIL
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --job-name="Ca_$ID"
#SBATCH --output=$DEST/slurm_%j.out
#SBATCH --time=12:00:00
#SBATCH --chdir=$DEST

ml load Miniconda3/4.9.2
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2
python fit_ca.py ca_config.json
WRAP
done

echo; echo "Built:"; ls -d $HOME_DIR/Ca_fit_c*
if [ "$1" == "submit" ]; then
  cd ~; for D in $HOME_DIR/Ca_fit_c*; do sbatch $D/run.sh; done
  squeue -u $USER -o "%.10i %.12j %.2t %.10M %C"
else
  echo; echo "Not submitted. Smoke test:  cd $HOME_DIR/Ca_fit_c02 && python fit_ca.py ca_config.json   (edit iters to 3 first for a quick check)"
  echo "then: bash $SRC/setup_ca_batch.sh submit"
fi
