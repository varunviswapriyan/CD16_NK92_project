#!/usr/bin/env bash
# run_thesis.sh -- set up and submit the thesis tumour-growth runs via Slurm.
# One line on the cluster:   cd ~/CD16_NK92_project && git pull && bash run_thesis.sh
#
# Submits male / female / male-no-sequestration, 3 seeds each (9 jobs), 1 week each.
# Override: SEEDS="8377 8378" STEPS=1008 TIME=24:00:00 PART=yourpartition bash run_thesis.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
SEEDS="${SEEDS:-8377 8378 8379}"
STEPS="${STEPS:-1008}"            # 1008 = day 14 -> 21
TIME="${TIME:-24:00:00}"
PART="${PART:-}"                  # set PART=... if your cluster needs a partition

[ -x ./spk_sexbias ] || { echo "build first: bash build_all.sh"; exit 1; }

sub() {  # cond  slide  sex  noseq(0/1)
  local cond=$1 slide=$2 sex=$3 noseq=$4 seed d
  for seed in $SEEDS; do
    d=runs/${cond}_s${seed}
    mkdir -p "$d/sex_bias_initialization_files"
    cp -f "$SRC/sex_bias_initialization_files/$slide" "$d/sex_bias_initialization_files/"
    # build in.simulation from Giuseppe's in.full
    sed -e "/^seed/ s/8377/$seed/" \
        -e "s|read_sites.*|read_sites sex_bias_initialization_files/$slide|" \
        -e "s|^run .*|run $STEPS|" \
        "$SRC/in.full" > "$d/in.simulation"
    # female rate edits (Giuseppe's: AR- recruit x1/0.63, AR+ recruit 0, Tp->Tm diff x2, Tp->AR+Tm diff 0)
    if [ "$sex" = female ]; then
      sed -i -e 's/6.87204/10.908/g' \
             -e 's/1.52712/0/' \
             -e '/differentiation to Tm/ s/0.00347/0.00694/' \
             -e '/differentiation to AR+ Tm/ s/0.00347/0/' "$d/in.simulation"
    fi
    # no-sequestration: zero the AR+ IFN-g uptake rate on the active app_style line
    [ "$noseq" = 1 ] && sed -i '/^app_style/ s/0.0016667/0/' "$d/in.simulation"
    # slurm job
    cat > "$d/job.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${cond}_s${seed}
#SBATCH --time=$TIME
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j.out
${PART:+#SBATCH --partition=$PART}
module load foss/2023a
cd "\$SLURM_SUBMIT_DIR"
mpirun -np 1 ../../spk_sexbias < in.simulation > log.spparks 2>&1
EOF
    ( cd "$d" && sbatch job.sh ) && echo "   submitted $cond seed $seed"
  done
}

sub male        d14_M_4_SPPARKS_init        male   0
sub female      d14_M_4_madeF_SPPARKS_init  female 0
sub male_noseq  d14_M_4_SPPARKS_init        male   1

echo
echo "watch jobs:   squeue -u \$USER"
echo "when done:    python3 thesis_repro/parse_and_plot.py runs"
