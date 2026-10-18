#!/usr/bin/env bash
# run_noseq.sh -- add the no-sequestration (AR-KO / IFNGR-blockade) run to the
# thesis-parameter set, so the figure shows male vs female vs male-no-seq.
# Cluster: cd ~/CD16_NK92_project && git pull && cp thesis_repro/run_noseq.sh . && bash run_noseq.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
SEEDS="${SEEDS:-8377 8378}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
[ -x "$BIN" ] || { echo "build first"; exit 1; }

# thesis-era parameters: capacity 1000, recruitment 0.1909 (same as the sweep's 'thesis')
sed -e '/^app_style/ s/10000/1000/' -e 's/6.87204/0.1909/g' -e 's/1.52712/0.0212/g' \
    "$SRC/in.full" > variant_thesis.in

for seed in $SEEDS; do
  d=runs_sweep/thesis/male_noseq_s$seed
  mkdir -p "$d/sex_bias_initialization_files"
  cp -f "$SRC/sex_bias_initialization_files/d14_M_4_SPPARKS_init" "$d/sex_bias_initialization_files/"
  # male slide, male rates, AR+ IFN-g uptake (sequestration) set to 0
  sed -e "/^seed/ s/[0-9][0-9]*/$seed/" \
      -e "s|read_sites.*|read_sites sex_bias_initialization_files/d14_M_4_SPPARKS_init|" \
      -e "/^app_style/ s/0.0016667/0/" \
      variant_thesis.in > "$d/in.simulation"
  cat > "$d/job.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=noseq_s$seed
#SBATCH --time=$TIME
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j.out
${PART:+#SBATCH --partition=$PART}
module load foss/2023a
cd "\$SLURM_SUBMIT_DIR"
mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
  ( cd "$d" && sbatch job.sh ) && echo "   submitted male_noseq seed $seed"
done

echo
echo "when done:  python3 ~/CD16_NK92_project/thesis_repro/parse_and_plot.py runs_sweep/thesis"
echo "            (the figure will now show male, female AND male_noseq)"
