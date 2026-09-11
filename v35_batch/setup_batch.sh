#!/bin/bash
# ============================================================================
# setup_batch.sh  -  build + (optionally) submit the hedged V35-V44 batch
#
#   bash setup_batch.sh          -> build all folders + wrappers (no submit)
#   bash setup_batch.sh submit   -> build, then sbatch all versions
#
# Edit BASE / HOME_DIR / resources below if your paths differ.
# ============================================================================
set -e

HOME_DIR=/home/gddaslab/vxv016
BASE=$HOME_DIR/NK92_fit_v25                     # proven-good folder to clone from
SRC=$(cd "$(dirname "$0")" && pwd)              # folder containing this script + the .py files
SUB=estimate_params_pzap_cleaned_up
MAIL=Varun.Viswapriyan@nationwidechildrens.org
PARTICLES=32
ITERS=90

# ---------- the 10 variants ----------
# name | LB | UB | SHAPE_MODE | RATIO_W | ORDER_W | HET_W | hetero adaptor0 | PSO opts
# adaptor0 "keep" = leave .bngl untouched (70.0). Otherwise the number is written into JJ_mixed_0_July7.bngl
VERSIONS=(
"35|[0.5,-3.3,0.5,0.5]|[1.8,-2.6,3.2,2.5]|own|30|5|1|keep|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"36|[0.5,-3.3,0.5,0.5]|[1.8,-2.6,3.2,2.5]|own|100|5|1|keep|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"37|[0.5,-3.3,0.5,0.5]|[1.8,-2.6,3.2,2.5]|gamma|30|5|2|keep|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"38|[0.5,-3.3,0.5,0.5]|[2.4,-2.6,3.2,2.5]|own|30|5|2|keep|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"39|[0.5,-3.3,0.5,0.5]|[1.8,-2.6,3.2,2.5]|own|30|20|1|keep|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"40|[0.5,-3.3,0.5,0.5]|[1.8,-2.6,3.2,2.5]|own|30|5|1|105.0|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"41|[0.5,-3.3,0.5,0.5]|[2.4,-2.6,3.2,2.5]|own|30|5|1|105.0|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"42|[1.4,-4.5,0.5,0.5]|[2.4,-2.3,2.6,2.5]|gamma|0|5|1|keep|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"43|[0.5,-3.3,0.5,0.5]|[1.8,-2.6,3.2,2.5]|own|30|5|1|147.0|{\"c1\":1.5,\"c2\":1.5,\"w\":0.5}"
"44|[0.0,-3.3,0.5,0.5]|[2.6,-2.6,3.2,2.5]|own|50|10|2|keep|{\"c1\":2.0,\"c2\":1.5,\"w\":0.7}"
)
# 35 own-norm core             36 own-norm heavy ratio       37 ratio penalty under ORIGINAL norm
# 38 own-norm, PSO picks regime 39 own-norm strong ordering  40 own-norm + hetero density 105
# 41 density 105 + free regime  42 reference (=V14 settings)  43 own-norm + density 147 (=gamma)
# 44 wide + exploratory PSO

for entry in "${VERSIONS[@]}"; do
  IFS='|' read -r V LB UB MODE RW OW HW DENS OPTS <<< "$entry"
  DEST=$HOME_DIR/NK92_fit_v$V
  echo "=== building v$V ==="
  rm -rf "$DEST"
  cp -r "$BASE" "$DEST"
  cp "$SRC/pzap_param_estimation_NK92.py" "$DEST/$SUB/"
  cp "$SRC/check_fit.py"                  "$DEST/$SUB/"
  rm -f "$DEST/$SUB/analysis_param_residue.dat" "$DEST/$SUB/plot_N.png"

  cat > "$DEST/$SUB/v_config.json" <<EOF
{
  "LB": $LB,
  "UB": $UB,
  "N_REPS": 3,
  "HETERO_WEIGHT": $HW,
  "ORDER_WEIGHT": $OW,
  "RATIO_WEIGHT": $RW,
  "SHAPE_MODE": "$MODE",
  "PSO_OPTIONS": $OPTS,
  "FIT_TMAX": 300.0
}
EOF

  if [ "$DENS" != "keep" ]; then
    BNGL=$DEST/$SUB/JJ_mixed_0_July7.bngl
    if grep -q "^adaptor0 70.0" "$BNGL"; then
      sed -i "s/^adaptor0 70.0/adaptor0 $DENS/" "$BNGL"
      echo "    patched hetero adaptor0 -> $DENS"
    else
      echo "    WARNING: 'adaptor0 70.0' not found in $BNGL - density NOT patched"
    fi
  fi

  cat > $HOME_DIR/run_v$V.sh <<EOF
#!/bin/bash
#SBATCH --mail-user=$MAIL
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --job-name="Fit_v$V"
#SBATCH --output=$HOME_DIR/slurm_v${V}_%j.out
#SBATCH --time=2-00:00:00
#SBATCH --chdir=$DEST/$SUB

unset BNGPATH
ml load Miniconda3/4.9.2
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2

python $DEST/$SUB/pzap_param_estimation_NK92.py $PARTICLES $ITERS
EOF
done

echo
echo "Built folders:"; ls -d $HOME_DIR/NK92_fit_v{35..44}
echo "Built wrappers:"; ls $HOME_DIR/run_v{35..44}.sh

if [ "$1" == "submit" ]; then
  echo; echo "Submitting..."
  for V in 35 36 37 38 39 40 41 42 43 44; do sbatch $HOME_DIR/run_v$V.sh; done
  squeue -u $USER -o "%.10i %.20j %.2t %.10M %.6D %C %m"
else
  echo; echo "Not submitted. Smoke-test first:"
  echo "  cd $HOME_DIR/NK92_fit_v35/$SUB && python pzap_param_estimation_NK92.py 4 1"
  echo "then:  bash $SRC/setup_batch.sh submit"
fi
