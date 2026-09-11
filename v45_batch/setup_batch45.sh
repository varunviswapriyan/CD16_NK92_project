#!/bin/bash
# ============================================================================
# setup_batch45.sh  -  build + (optionally) submit the structural-fix batch V45-V54
#
#   bash setup_batch45.sh          -> build all folders + wrappers (no submit)
#   bash setup_batch45.sh submit   -> build, then sbatch all versions
#
# Base folder is V42 (best shape fit so far). Each variant changes the MODEL
# (.bngl) and/or the PSO search region, not just the cost function.
#
# Levers:
#   ZAPLB   - log10 lower bound on ZAP0. 2.0 forces PSO out of the ZAP-limited
#             regime where hetero/gamma is pinned at 1.0 (v42: Z0=39 -> ~980 ZAP
#             molecules vs thousands of phospho-ITAMs).
#   DENS    - hetero adaptor density in JJ_mixed (70 = untouched). At receptor
#             saturation hetero/gamma = DENS*4 / (147*2); 70 -> 0.95, 80 -> 1.09.
#   FRAC    - gamma-type ITAMs get ZAP70 on-rate kzb_g = kzb*FRAC ("keep" = no
#             patch). Applied to both Gamma ITAMs in JJ_gamma and to Hetero's
#             ITAM4 in JJ_mixed. Zeta-type ITAMs keep kzb. Biology: ZAP70 tandem
#             SH2 prefers zeta ITAMs over the gamma ITAM.
# ============================================================================
set -e
cd ~

HOME_DIR=/home/gddaslab/vxv016
BASE=$HOME_DIR/NK92_fit_v42                     # clone from here
SRC=$(cd "$(dirname "$0")" && pwd)              # folder containing this script + patch_affinity.py
SUB=estimate_params_pzap_cleaned_up
MAIL=Varun.Viswapriyan@nationwidechildrens.org
PARTICLES=32
ITERS=90

GAMMA_BNGL=JJ_gamma_0_July2.bngl
MIXED_BNGL=JJ_mixed_0_July7.bngl
GAMMA_RULES="r9_PZAP_binding_ITAM1 r10_PZAP_binding_ITAM2 r7_UZAP_binding_ITAM1 r8_UZAP_binding_ITAM2"
MIXED_RULES="rh14_PZAP_binding_ITAM4 rh18_UZAP_binding_ITAM4 rg9_PZAP_binding_ITAM1 rg10_PZAP_binding_ITAM2 rg7_UZAP_binding_ITAM1 rg8_UZAP_binding_ITAM2"

# v42 cost settings (kept for most variants): gamma-normalised, hetero weight 5, no penalties
# ----------------------------------------------------------------------------
# V  | ZAPLB | DENS | FRAC | SHAPE | RATIO_W | ORDER_W | HET_W | what it tests
# ----------------------------------------------------------------------------
VARIANTS="
45 | 2.0 | 70   | keep | gamma | 0  | 0 | 5 | regime only: force ZAP0 >= 100
46 | 2.0 | 80   | keep | gamma | 0  | 0 | 5 | regime + density 80
47 | 2.0 | 90   | keep | gamma | 0  | 0 | 5 | regime + density 90
48 | 2.0 | 70   | 0.5  | gamma | 0  | 0 | 5 | regime + gamma-ITAM ZAP affinity halved
49 | 2.0 | 70   | 0.25 | gamma | 0  | 0 | 5 | regime + gamma-ITAM ZAP affinity quartered
50 | 2.0 | 80   | 0.5  | gamma | 0  | 0 | 5 | regime + density 80 + affinity 0.5
51 | 0.5 | 70   | 0.5  | gamma | 0  | 0 | 5 | affinity 0.5 alone, v42 bounds (no regime push)
52 | 2.0 | 70   | 0.5  | own   | 30 | 5 | 1 | as v48 but own-norm + ratio penalty (locks ratio now that model can reach it)
53 | 2.3 | 70   | 0.5  | gamma | 0  | 0 | 5 | stronger regime push (ZAP0 >= 200) + affinity 0.5
54 | 2.0 | 80   | 0.5  | own   | 30 | 5 | 1 | full combo: density 80 + affinity 0.5 + own-norm + ratio penalty
"

[ -d "$BASE/$SUB" ] || { echo "BASE not found: $BASE/$SUB"; exit 1; }
[ -f "$SRC/patch_affinity.py" ] || { echo "patch_affinity.py not found next to this script"; exit 1; }

echo "$VARIANTS" | grep -E '^[0-9]' | while IFS='|' read -r V ZAPLB DENS FRAC SHAPE RATIOW ORDERW HETW DESC; do
  V=$(echo $V); ZAPLB=$(echo $ZAPLB); DENS=$(echo $DENS); FRAC=$(echo $FRAC)
  SHAPE=$(echo $SHAPE); RATIOW=$(echo $RATIOW); ORDERW=$(echo $ORDERW); HETW=$(echo $HETW)
  DEST=$HOME_DIR/NK92_fit_v$V
  echo "=== building v$V : $(echo $DESC) ==="

  rm -rf "$DEST"
  cp -r "$BASE" "$DEST"
  cd "$DEST/$SUB"
  rm -rf zeta_runs gamma_runs mixed_runs analysis_param_residue.dat plot_N.png v_config.json

  # --- PSO config (log10 bounds: [lig0, kd10, ZAP0, SYK0]) ---
  # ZAP0 upper bound raised to 3.2 so the pushed regime has room.
  cat > v_config.json <<EOF
{
  "LB": [1.4, -4.5, $ZAPLB, 0.5],
  "UB": [2.4, -2.3, 3.2, 2.5],
  "N_REPS": 3,
  "HETERO_WEIGHT": $HETW,
  "ORDER_WEIGHT": $ORDERW,
  "RATIO_WEIGHT": $RATIOW,
  "SHAPE_MODE": "$SHAPE",
  "PSO_OPTIONS": {"c1": 1.5, "c2": 1.5, "w": 0.5},
  "FIT_TMAX": 300.0,
  "NOTE": "v$V: ZAPLB=$ZAPLB DENS=$DENS FRAC=$FRAC"
}
EOF

  # --- hetero density patch ---
  if [ "$DENS" != "70" ]; then
    if grep -q "^adaptor0 70.0" "$MIXED_BNGL"; then
      sed -i "s/^adaptor0 70.0/adaptor0 $DENS.0/" "$MIXED_BNGL"
      echo "    patched hetero adaptor0 -> $DENS.0"
    else
      echo "    ERROR: 'adaptor0 70.0' not found in $MIXED_BNGL"; exit 1
    fi
  fi

  # --- zeta-vs-gamma ITAM ZAP70 affinity patch ---
  if [ "$FRAC" != "keep" ]; then
    python "$SRC/patch_affinity.py" "$GAMMA_BNGL" "$FRAC" $GAMMA_RULES
    python "$SRC/patch_affinity.py" "$MIXED_BNGL" "$FRAC" $MIXED_RULES
    # confirm pyBioNetGen still parses the patched models (needs CD16_v2 env active)
    python - <<EOF || { echo "    ERROR: patched .bngl failed to parse - is conda env CD16_v2 active?"; exit 1; }
import bionetgen
for f in ["$GAMMA_BNGL", "$MIXED_BNGL"]:
    m = bionetgen.bngmodel(f)
    assert "kzb_g" in str(m), f + " lost kzb_g after parse"
print("    parse check OK (kzb_g present in both patched models)")
EOF
  fi

  # --- SLURM wrapper ---
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
  cd ~
done

echo
echo "Built folders:"; ls -d $HOME_DIR/NK92_fit_v{45..54}
echo "Built wrappers:"; ls $HOME_DIR/run_v{45..54}.sh
echo
echo "Model patches applied:"
kzbg() { local x; x=$(grep -oE '^kzb_g .*' "$1" | awk '{print $NF}'); echo "${x:-none}"; }
for V in 45 46 47 48 49 50 51 52 53 54; do
  D=$HOME_DIR/NK92_fit_v$V/$SUB
  printf "  v%s  hetero adaptor0=%s  kzb_g(gamma)=%s  kzb_g(mixed)=%s\n" "$V" \
    "$(grep -oE '^adaptor0 [0-9.]+' $D/$MIXED_BNGL | cut -d' ' -f2)" \
    "$(kzbg $D/$GAMMA_BNGL)" "$(kzbg $D/$MIXED_BNGL)"
done

if [ "$1" == "submit" ]; then
  echo; echo "Submitting..."
  cd ~
  for V in 45 46 47 48 49 50 51 52 53 54; do sbatch $HOME_DIR/run_v$V.sh; done
  squeue -u $USER -o "%.10i %.20j %.2t %.10M %.6D %C %m"
else
  echo; echo "Not submitted. Smoke-test the patched model first:"
  echo "  cd $HOME_DIR/NK92_fit_v48/$SUB && python pzap_param_estimation_NK92.py 4 1"
  echo "then:  bash $SRC/setup_batch45.sh submit"
fi
