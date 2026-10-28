#!/usr/bin/env bash
# run_overnight.sh -- submit the two SLOW figure runs overnight:
#   (1) no-sequestration male+female (n=NSEED)   -> IFNGR-blocking fold-change figure (Fig4 #4 / SI #2)
#   (2) spatial male+female (1 seed, with dumps)  -> spatial profiles figure (Fig4 #2)
#
# Submit:   cd ~/CD16_NK92_project && git pull && NSEED=40 bash run_overnight.sh
# Morning:  bash run_overnight.sh agg
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
NSEED="${NSEED:-40}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"

# ---- aggregate / plot mode ---------------------------------------------------
if [ "${1:-}" = "agg" ]; then
  source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
  echo "=== no-seq fold changes (IFNGR-blocking figure) ==="
  python3 - <<'PY'
import glob, math
def fold(f):
    col=first=last=None
    for line in open(f):
        s=line.split()
        if not s: continue
        if "sum(i1)" in s: col=s.index("sum(i1)"); continue
        if col is not None and len(s)>col:
            try:
                float(s[0]); v=float(s[col])
                if first is None: first=v
                last=v
            except ValueError: pass
    return last/first if first and last else None
for cond in ("male_noseq","female_noseq"):
    vs=[x for x in (fold(f) for f in sorted(glob.glob(f"runs_big/noseq/{cond}/s*/log.spparks"))) if x]
    if vs:
        m=sum(vs)/len(vs); sd=math.sqrt(sum((x-m)**2 for x in vs)/(len(vs)-1)) if len(vs)>1 else 0
        print(f"{cond:14s} n={len(vs):3d} mean fold={m:.4f} SD={sd:.4f}")
print("compare standard (sequestration ON):  male 2.3793  female 1.5456")
PY
  echo
  echo "=== spatial plots ==="
  python3 ~/CD16_NK92_project/thesis_repro/plot_spatial.py runs_spatial_prod && echo "-> PNGs in runs_spatial_prod/"
  exit 0
fi

# ---- submit mode -------------------------------------------------------------
[ -x "$BIN" ] || { echo "build first"; exit 1; }

# (1) no-seq m/f, NSEED seeds  (AR+ IFN-g internalization 0.0016667 -> 0)
for cond in male_noseq female_noseq; do
  case $cond in
    male_noseq)   slide=d14_M_4_SPPARKS_init; fem=0;;
    female_noseq) slide=d14_F_4_SPPARKS_init; fem=1;;
  esac
  base=runs_big/noseq/$cond
  mkdir -p "$base"; cp -f "$SRC/sex_bias_initialization_files/$slide" "$base/$slide"
  SEDPROG="/^seed/ s/[0-9][0-9]*/\$seed/ ; s|read_sites.*|read_sites sex_bias_initialization_files/$slide| ; /^app_style/ s/0.0016667/0/"
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
sed "${SEDPROG}" "$PWD/$SRC/in.full" > "\$d/in.simulation"
if [ "$fem" = 1 ]; then
  awk '{ if (\$0 ~ /recruitment/ && \$0 !~ /app_style/){ if(\$0 ~ /AR\\+/)\$6=0; else if(\$0 ~ /Tp recruitment/||\$0 ~ /Tm recruitment/)\$6=\$6/0.63 }
         if (\$0 ~ /Tp differentiation/ && \$0 !~ /app_style/){ if(\$0 ~ /AR\\+ Tm/)\$6=0; else \$6=\$6*2 } print }' \\
    "\$d/in.simulation" > "\$d/in.tmp" && mv "\$d/in.tmp" "\$d/in.simulation"
fi
cd "\$d"; mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
  sbatch "$base/array.sh" && echo "   submitted $cond x$NSEED"
done

# (2) spatial m/f, single seed, with dumps at day 14 / 17.5 / 21
SEED=8377
for spec in "male:d14_M_4_SPPARKS_init:male" "female:d14_F_4_SPPARKS_init:female"; do
  IFS=: read -r cond slide sex <<<"$spec"
  d=runs_spatial_prod/$cond
  mkdir -p "$d/sex_bias_initialization_files"
  cp -f "$SRC/sex_bias_initialization_files/$slide" "$d/sex_bias_initialization_files/"
  sed -e "/^seed/ s/[0-9][0-9]*/$SEED/" \
      -e "s|read_sites.*|read_sites sex_bias_initialization_files/$slide|" \
      "$SRC/in.full" > "$d/in.simulation"
  sed -i 's|^#dump .*|dump 1 sites 504.0 dump.2D_model.* id x y i1 i2 i3 i4 i5 i6 i7 d1|' "$d/in.simulation"
  if [ "$sex" = female ]; then
    awk '{ if ($0 ~ /recruitment/ && $0 !~ /app_style/){ if($0 ~ /AR\+/)$6=0; else if($0 ~ /Tp recruitment/||$0 ~ /Tm recruitment/)$6=$6/0.63 }
           if ($0 ~ /Tp differentiation/ && $0 !~ /app_style/){ if($0 ~ /AR\+ Tm/)$6=0; else $6=$6*2 } print }' \
      "$d/in.simulation" > "$d/in.tmp" && mv "$d/in.tmp" "$d/in.simulation"
  fi
  cat > "$d/job.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=spatial_$cond
#SBATCH --time=$TIME
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j.out
${PART:+#SBATCH --partition=$PART}
module load foss/2023a
cd "\$SLURM_SUBMIT_DIR"
mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
  ( cd "$d" && sbatch job.sh ) && echo "   submitted spatial $cond"
done

echo
echo "watch:    squeue -u \$USER   (wait until EMPTY)"
echo "morning:  bash run_overnight.sh agg"
