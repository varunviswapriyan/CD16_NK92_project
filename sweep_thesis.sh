#!/usr/bin/env bash
# sweep_thesis.sh -- submit a parameter sweep tonight: male vs REAL female,
# across 4 parameter sets, to bracket Giuseppe's thesis tumour-growth shape.
# One line on the cluster:  cd ~/CD16_NK92_project && git pull && bash sweep_thesis.sh
#
# Override: SEEDS="8377 8378" TIME=24:00:00 PART=yourpartition bash sweep_thesis.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
SEEDS="${SEEDS:-8377 8378}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
PSETS="${PSETS:-prod cap1k recr019 thesis}"

[ -x "$BIN" ] || { echo "build first: bash build_all.sh"; exit 1; }

# build a parameter-variant in.full
make_infull() {  # pset -> prints the variant file path
  local p=$1 f=$SRC/in.full o=variant_$1.in
  case $p in
    prod)    cp "$f" "$o";;
    cap1k)   sed '/^app_style/ s/10000/1000/' "$f" > "$o";;
    recr019) sed -e 's/6.87204/0.1909/g' -e 's/1.52712/0.0212/g' "$f" > "$o";;
    thesis)  sed -e '/^app_style/ s/10000/1000/' -e 's/6.87204/0.1909/g' -e 's/1.52712/0.0212/g' "$f" > "$o";;
    *) cp "$f" "$o";;
  esac
  echo "$o"
}

# build in.simulation (seed + slide, and Giuseppe's female transform if female)
make_sim() {  # out infull slide sex seed
  local out=$1 infull=$2 slide=$3 sex=$4 seed=$5
  sed -e "/^seed/ s/[0-9][0-9]*/$seed/" \
      -e "s|read_sites.*|read_sites sex_bias_initialization_files/$slide|" \
      "$infull" > "$out"
  if [ "$sex" = female ]; then
    awk '{
      if ($0 ~ /recruitment/ && $0 !~ /app_style/) {
        if ($0 ~ /AR\+/) $6=0;
        else if ($0 ~ /Tp recruitment/ || $0 ~ /Tm recruitment/) $6=$6/0.63;
      }
      if ($0 ~ /Tp differentiation/ && $0 !~ /app_style/) {
        if ($0 ~ /AR\+ Tm/) $6=0; else $6=$6*2;
      }
      print
    }' "$out" > "$out.tmp" && mv "$out.tmp" "$out"
  fi
}

n=0
for pset in $PSETS; do
  infull=$(make_infull "$pset")
  for seed in $SEEDS; do
    for spec in "male:d14_M_4_SPPARKS_init:male" "female:d14_F_4_SPPARKS_init:female"; do
      IFS=: read -r cond slide sex <<<"$spec"
      d=runs_sweep/$pset/${cond}_s${seed}
      mkdir -p "$d/sex_bias_initialization_files"
      cp -f "$SRC/sex_bias_initialization_files/$slide" "$d/sex_bias_initialization_files/"
      make_sim "$d/in.simulation" "$infull" "$slide" "$sex" "$seed"
      cat > "$d/job.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${pset}_${cond}_s${seed}
#SBATCH --time=$TIME
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j.out
${PART:+#SBATCH --partition=$PART}
module load foss/2023a
cd "\$SLURM_SUBMIT_DIR"
mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
      ( cd "$d" && sbatch job.sh ) && { echo "   submitted $pset $cond s$seed"; n=$((n+1)); }
    done
  done
done

echo
echo "== submitted $n jobs. watch: squeue -u \$USER"
echo "== when done, plot each set:"
for p in $PSETS; do echo "   python3 ~/CD16_NK92_project/thesis_repro/parse_and_plot.py runs_sweep/$p"; done
