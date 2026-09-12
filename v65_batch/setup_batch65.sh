#!/bin/bash
# ============================================================================
# setup_batch65.sh  -  build + (optionally) submit the V65-V74 batch
#
#   bash setup_batch65.sh          -> build all folders + wrappers (no submit)
#   bash setup_batch65.sh submit   -> build, then sbatch all versions
#
# Base = V42 folder (clean Indrani .bngl files). Every version gets:
#   * the N-parameter fit script + check_fit (from this folder)
#   * calculate_SSR_*.py patched with the 'extra' hook
#   * .bngl files patched with fittable multipliers (KZP_MULT, KPR_MULT, KZU_MULT,
#     KP_MULT, KZD_MULT, KDL_MULT, KZBG_FRAC) - all default 1.0 = original model
# Per-version switches: gamma-ITAM affinity split (frac), per-receptor CBL degradation
# (recdeg), hetero density, and WHICH multipliers are fitted.
# ============================================================================
set -e
cd ~

HOME_DIR=/home/gddaslab/vxv016
BASE=$HOME_DIR/NK92_fit_v42
SRC=$(cd "$(dirname "$0")" && pwd)
SUB=estimate_params_pzap_cleaned_up
MAIL=Varun.Viswapriyan@nationwidechildrens.org
PARTICLES=32
ITERS=90
PSO_STD='{"c1":1.5,"c2":1.5,"w":0.5}'
PSO_EXP='{"c1":2.0,"c2":1.5,"w":0.7}'

# log10 bounds used for the multipliers
#   KZBG_FRAC [-1,0]  KZP_MULT [-0.3,1]  others [-1,1]  mixed:adaptor0 [1.7,2.1] (50-126)
# -----------------------------------------------------------------------------------------------------------
# V | frac | recdeg | dens | extra params | extra LB | extra UB | FIXED | shape | RW | OW | HW | ZAPLB | PSO | description
# -----------------------------------------------------------------------------------------------------------
VARIANTS="
65 | y | percond  | keep | KZBG_FRAC,KZP_MULT,gamma:KDL_MULT,zeta:KDL_MULT          | -1,-0.3,-1,-1      | 0,1,1,1      | {} | own   | 30 | 5 | 1 | 2.0 | STD | rule untouched: per-condition kdl (=v60 recipe, density 70)
66 | y | percond  | keep | KZBG_FRAC,KZP_MULT,gamma:KDL_MULT,zeta:KDL_MULT,KPR_MULT | -1,-0.3,-1,-1,-1   | 0,1,1,1,1    | {} | own   | 30 | 5 | 1 | 2.0 | STD | per-condition kdl + KPR (v58 recipe without rule change)
67 | y | itamnorm | keep | KZBG_FRAC,KZP_MULT                                       | -1,-0.3            | 0,1          | {} | own   | 30 | 5 | 1 | 2.0 | STD | rule untouched: kdl normalized by ITAM count (no extra params)
68 | y | itamnorm | keep | KZBG_FRAC,KZP_MULT,KPR_MULT                              | -1,-0.3,-1         | 0,1,1        | {} | own   | 30 | 5 | 1 | 2.0 | STD | ITAM-normalized kdl + KPR
69 | y | itamnorm | keep | KZBG_FRAC,KZP_MULT,KPR_MULT                              | -1,-0.3,-1         | 0,1,1        | {} | gamma | 0  | 0 | 5 | 2.0 | STD | as 68, v57-style cost (no penalties)
70 | y | percond  | keep | KZBG_FRAC,KZP_MULT,gamma:KDL_MULT,zeta:KDL_MULT,KPR_MULT,KZD_MULT | -1,-0.3,-1,-1,-1,-1 | 0,1,1,1,1,1 | {} | own | 30 | 5 | 1 | 2.0 | STD | per-condition kdl + KPR + SHP dephosph
71 | y | delmol   | keep | KZBG_FRAC,KZP_MULT                                       | -1,-0.3            | 0,1          | {} | own   | 30 | 5 | 1 | 2.0 | STD | DeleteMolecules: Cbl-b removes only the tagged kinase
72 | y | delmol   | keep | KZBG_FRAC,KZP_MULT,KPR_MULT                              | -1,-0.3,-1         | 0,1,1        | {} | own   | 30 | 5 | 1 | 2.0 | STD | DeleteMolecules + KPR
73 | y | itamnorm | keep | KZBG_FRAC,KZP_MULT,KPR_MULT                              | -1,-0.3,-1         | 0,1,1        | {} | own   | 30 | 5 | 1 | 2.0 | STD | replicate of 68 (PSO seed)
74 | y | percond  | keep | KZBG_FRAC,KZP_MULT,gamma:KDL_MULT,zeta:KDL_MULT,KPR_MULT | -1,-0.3,-1,-1,-1   | 0,1,1,1,1    | {} | own   | 30 | 5 | 1 | 2.0 | STD | replicate of 66 (PSO seed)
"

for f in pzap_param_estimation_NK92.py check_fit.py patch_model.py patch_ssr_extra.py make_config.py; do
  [ -f "$SRC/$f" ] || { echo "missing $SRC/$f"; exit 1; }
done
[ -d "$BASE/$SUB" ] || { echo "BASE not found: $BASE/$SUB"; exit 1; }

echo "$VARIANTS" | grep -E '^[0-9]' | while IFS='|' read -r V FRAC DECAY DENS EP ELB EUB FIXED SHAPE RW OW HW ZAPLB PSO DESC; do
  V=$(echo $V); FRAC=$(echo $FRAC); DECAY=$(echo $DECAY); DENS=$(echo $DENS)
  EP=$(echo $EP); ELB=$(echo $ELB); EUB=$(echo $EUB); FIXED=$(echo $FIXED)
  SHAPE=$(echo $SHAPE); RW=$(echo $RW); OW=$(echo $OW); HW=$(echo $HW); ZAPLB=$(echo $ZAPLB); PSO=$(echo $PSO)
  DESC=$(echo $DESC)
  if [ "$PSO" == "EXP" ]; then PSOJ=$PSO_EXP; else PSOJ=$PSO_STD; fi
  DEST=$HOME_DIR/NK92_fit_v$V
  echo "=== building v$V : $DESC ==="

  rm -rf "$DEST"
  cp -r "$BASE" "$DEST"
  cd "$DEST/$SUB"
  rm -rf zeta_runs gamma_runs mixed_runs analysis_param_residue.dat plot_N.png v_config.json
  cp "$SRC/pzap_param_estimation_NK92.py" "$SRC/check_fit.py" .

  for k in zeta gamma mixed; do python "$SRC/patch_ssr_extra.py" calculate_SSR_$k.py; done

  FL=""
  [ "$FRAC" == "y" ] && FL="$FL --frac"
  case "$DECAY" in
    percond)  ;;                          # rule untouched; per-condition KDL_MULT fitted via config
    itamnorm) FL="$FL --itamnorm" ;;      # rule untouched; kdl / ITAM count
    delmol)   FL="$FL --delmol" ;;        # rule kept, DeleteMolecules added
    recdeg)   FL="$FL --recdeg" ;;
    *) echo "    ERROR: unknown decay mode $DECAY"; exit 1 ;;
  esac
  python "$SRC/patch_model.py" JJ_zeta_0_July2.bngl  zeta  $FL
  python "$SRC/patch_model.py" JJ_gamma_0_July2.bngl gamma $FL
  python "$SRC/patch_model.py" JJ_mixed_0_July7.bngl mixed $FL
  if [ "$DENS" != "keep" ]; then
    grep -q "^adaptor0 70.0" JJ_mixed_0_July7.bngl || { echo "    ERROR: adaptor0 70.0 not found"; exit 1; }
    sed -i "s/^adaptor0 70.0/adaptor0 $DENS.0/" JJ_mixed_0_July7.bngl
    echo "    hetero adaptor0 -> $DENS.0"
  fi

  python "$SRC/make_config.py" v_config.json "$EP" "$ELB" "$EUB" "$FIXED" "$SHAPE" "$RW" "$OW" "$HW" "$ZAPLB" "$PSOJ" "v$V: $DESC"

  python - <<'PYEOF' || { echo "    ERROR: patched .bngl failed to parse - is conda env CD16_v2 active?"; exit 1; }
import bionetgen
for f in ["JJ_zeta_0_July2.bngl", "JJ_gamma_0_July2.bngl", "JJ_mixed_0_July7.bngl"]:
    m = bionetgen.bngmodel(f); s = str(m)
    for p in ["KZP_MULT", "KPR_MULT", "KZU_MULT", "KP_MULT", "KZD_MULT", "KDL_MULT", "KZBG_FRAC"]:
        assert p in s, f + " lost " + p
    raw = open(f).read()
    if "DeleteMolecules" in raw: assert "DeleteMolecules" in s, f + " lost DeleteMolecules after parse"
    if "ITAM_N" in raw: assert "ITAM_N" in s, f + " lost ITAM_N after parse"
print("    parse check OK (all multipliers present in 3 models)")
PYEOF

  cat > $HOME_DIR/run_v$V.sh <<WRAP_EOF
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
WRAP_EOF
  cd ~
done

echo
echo "Built folders:";  ls -d $HOME_DIR/NK92_fit_v{65..74}
echo "Built wrappers:"; ls $HOME_DIR/run_v{65..74}.sh
echo
echo "Summary:"
for V in 65 66 67 68 69 70 71 72 73 74; do
  D=$HOME_DIR/NK92_fit_v$V/$SUB
  M=$D/JJ_mixed_0_July7.bngl
  if grep -q DeleteMolecules $M; then DM=delmol; elif grep -q '^ITAM_N' $M; then DM=itamnorm; elif grep -q '^rdeg_hetero' $M; then DM=recdeg; else DM=original; fi
  printf "  v%s  dens=%-5s kzb_g=%-3s CBL=%-9s %s\n" "$V" \
    "$(grep -oE '^adaptor0 [0-9.]+' $M | cut -d' ' -f2)" \
    "$(grep -q '^kzb_g' $D/JJ_gamma_0_July2.bngl && echo yes || echo no)" "$DM" \
    "$(python -c "import json;c=json.load(open('$D/v_config.json'));print(c['PARAMS'][4:], c['SHAPE_MODE'])")"
done

if [ "$1" == "submit" ]; then
  echo; echo "Submitting..."; cd ~
  for V in 65 66 67 68 69 70 71 72 73 74; do sbatch $HOME_DIR/run_v$V.sh; done
  squeue -u $USER -o "%.10i %.20j %.2t %.10M %.6D %C %m"
else
  echo; echo "Not submitted. Smoke-test BOTH new mechanisms first (itamnorm=v67, delmol=v71):"
  echo "  cd $HOME_DIR/NK92_fit_v67/$SUB && python pzap_param_estimation_NK92.py 4 1"
  echo "  cd $HOME_DIR/NK92_fit_v71/$SUB && python pzap_param_estimation_NK92.py 4 1"
  echo "then:  bash $SRC/setup_batch65.sh submit"
fi
