#!/usr/bin/env bash
# run_noseq_prod.sh -- add the no-sequestration (AR-KO / IFNGR-blockade) growth run
# at PROD params (cap 10,000, recruitment 6.87204), so the prod growth figure shows
# male vs female vs male-no-seq. Uses in.full as-is with AR+ IFN-g uptake set to 0.
# Cluster: cd ~/CD16_NK92_project && git pull && cp thesis_repro/run_noseq_prod.sh . && bash run_noseq_prod.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
SEEDS="${SEEDS:-8377 8378}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
[ -x "$BIN" ] || { echo "build first"; exit 1; }

for seed in $SEEDS; do
  d=runs_sweep/prod/male_noseq_s$seed
  mkdir -p "$d/sex_bias_initialization_files"
  cp -f "$SRC/sex_bias_initialization_files/d14_M_4_SPPARKS_init" "$d/sex_bias_initialization_files/"
  # male slide, prod rates, AR+ IFN-g uptake (sequestration) set to 0
  sed -e "/^seed/ s/[0-9][0-9]*/$seed/" \
      -e "s|read_sites.*|read_sites sex_bias_initialization_files/d14_M_4_SPPARKS_init|" \
      -e "/^app_style/ s/0.0016667/0/" \
      "$SRC/in.full" > "$d/in.simulation"
  cat > "$d/job.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=noseqP_s$seed
#SBATCH --time=$TIME
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j.out
${PART:+#SBATCH --partition=$PART}
module load foss/2023a
cd "\$SLURM_SUBMIT_DIR"
mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
  ( cd "$d" && sbatch job.sh ) && echo "   submitted prod male_noseq seed $seed"
done

echo
echo "when done:  python3 ~/CD16_NK92_project/thesis_repro/parse_and_plot.py runs_sweep/prod"
echo "            (the prod growth figure will then show male, female AND male_noseq)"
