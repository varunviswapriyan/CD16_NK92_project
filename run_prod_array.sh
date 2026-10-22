#!/usr/bin/env bash
# run_prod_array.sh -- launch NSEED replicate simulations per condition (male,
# female, male_noseq) at PROD params using a Slurm job array (all in parallel).
# Seeds = 9001..9000+NSEED, written under runs_big/prod/<cond>/s<seed>/.
# Run:  cd ~/CD16_NK92_project && git pull && NSEED=40 bash run_prod_array.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
NSEED="${NSEED:-40}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
[ -x "$BIN" ] || { echo "build first"; exit 1; }

for cond in male female male_noseq; do
  case $cond in
    male)       slide=d14_M_4_SPPARKS_init; fem=0; noseq=0;;
    female)     slide=d14_F_4_SPPARKS_init; fem=1; noseq=0;;
    male_noseq) slide=d14_M_4_SPPARKS_init; fem=0; noseq=1;;
  esac
  base=runs_big/prod/$cond
  mkdir -p "$base"
  cp -f "$SRC/sex_bias_initialization_files/$slide" "$base/$slide"

  # build the sed program as ONE semicolon-joined script (no quote-stripping trap);
  # \$seed stays literal so it expands inside the job at runtime.
  SEDPROG="/^seed/ s/[0-9][0-9]*/\$seed/ ; s|read_sites.*|read_sites sex_bias_initialization_files/$slide|"
  NOSEQ=""
  [ "$noseq" = 1 ] && NOSEQ=" ; /^app_style/ s/0.0016667/0/"

  cat > "$base/array.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${cond}
#SBATCH --array=1-$NSEED
#SBATCH --time=$TIME
#SBATCH --ntasks=1
#SBATCH --output=$base/slurm_%A_%a.out
${PART:+#SBATCH --partition=$PART}
module load foss/2023a
seed=\$((9000 + SLURM_ARRAY_TASK_ID))
d=$base/s\$seed
mkdir -p "\$d/sex_bias_initialization_files"
cp -f "$PWD/$base/$slide" "\$d/sex_bias_initialization_files/"
sed "${SEDPROG}${NOSEQ}" "$PWD/$SRC/in.full" > "\$d/in.simulation"
if [ "$fem" = 1 ]; then
  awk '{ if (\$0 ~ /recruitment/ && \$0 !~ /app_style/){ if(\$0 ~ /AR\\+/)\$6=0; else if(\$0 ~ /Tp recruitment/||\$0 ~ /Tm recruitment/)\$6=\$6/0.63 }
         if (\$0 ~ /Tp differentiation/ && \$0 !~ /app_style/){ if(\$0 ~ /AR\\+ Tm/)\$6=0; else \$6=\$6*2 } print }' \\
    "\$d/in.simulation" > "\$d/in.tmp" && mv "\$d/in.tmp" "\$d/in.simulation"
fi
cd "\$d"
mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
  sbatch "$base/array.sh" && echo "   submitted $cond  x$NSEED seeds"
done

echo
echo "watch:  squeue -u \$USER"
echo "when done:  bash aggregate.sh"
