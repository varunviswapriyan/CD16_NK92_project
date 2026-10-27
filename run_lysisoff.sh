#!/usr/bin/env bash
# run_lysisoff.sh -- control: turn OFF direct T-cell lysis of cancer.
#   Sets the mature-CD8 (Tm) and AR+ mature (ARTm) cancer-lysis rate 0.06 -> 0
#   on both male and female. (Exhausted-T lysis is already 0.) All other
#   parameters = prod/manuscript. Tests whether IFN-g, not lysis, is the
#   dominant brake on tumor growth: if growth barely rises, lysis is minor.
#
# Submit:     cd ~/CD16_NK92_project && git pull && NSEED=40 bash run_lysisoff.sh
# Aggregate:  bash run_lysisoff.sh agg       (after jobs finish)
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
NSEED="${NSEED:-40}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
ROOT=runs_big/lysisoff

# ---- aggregate mode ----------------------------------------------------------
if [ "${1:-}" = "agg" ]; then
  source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
  python3 - "$ROOT" <<'PY'
import glob, math, sys
root=sys.argv[1]
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
for cond in ("male_nolysis","female_nolysis"):
    vs=[x for x in (fold(f) for f in sorted(glob.glob(f"{root}/{cond}/s*/log.spparks"))) if x]
    if vs:
        m=sum(vs)/len(vs)
        sd=math.sqrt(sum((x-m)**2 for x in vs)/(len(vs)-1)) if len(vs)>1 else 0
        se=sd/math.sqrt(len(vs)) if len(vs)>1 else 0
        print(f"{cond:16s} n={len(vs):3d}  mean fold={m:.4f}  SD={sd:.4f}  SE={se:.4f}")
print("\ncompare to standard (lysis ON):  male 2.3793   female 1.5456")
PY
  exit 0
fi

# ---- submit mode -------------------------------------------------------------
[ -x "$BIN" ] || { echo "build first"; exit 1; }
for cond in male_nolysis female_nolysis; do
  case $cond in
    male_nolysis)   slide=d14_M_4_SPPARKS_init; fem=0;;
    female_nolysis) slide=d14_F_4_SPPARKS_init; fem=1;;
  esac
  base=$ROOT/$cond
  mkdir -p "$base"
  cp -f "$SRC/sex_bias_initialization_files/$slide" "$base/$slide"
  SEDPROG="/^seed/ s/[0-9][0-9]*/\$seed/ ; s|read_sites.*|read_sites sex_bias_initialization_files/$slide|"
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
# turn off cancer lysis (0.06 -> 0) on every reaction whose comment says "lysis"
awk '/lysis/ && \$0 !~ /app_style/ { for(i=1;i<=NF;i++) if(\$i=="0.06") \$i="0" } {print}' \\
    "\$d/in.simulation" > "\$d/in.tmp" && mv "\$d/in.tmp" "\$d/in.simulation"
if [ "$fem" = 1 ]; then
  awk '{ if (\$0 ~ /recruitment/ && \$0 !~ /app_style/){ if(\$0 ~ /AR\\+/)\$6=0; else if(\$0 ~ /Tp recruitment/||\$0 ~ /Tm recruitment/)\$6=\$6/0.63 }
         if (\$0 ~ /Tp differentiation/ && \$0 !~ /app_style/){ if(\$0 ~ /AR\\+ Tm/)\$6=0; else \$6=\$6*2 } print }' \\
    "\$d/in.simulation" > "\$d/in.tmp" && mv "\$d/in.tmp" "\$d/in.simulation"
fi
cd "\$d"
mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
  sbatch "$base/array.sh" && echo "   submitted $cond x$NSEED"
done
echo
echo "watch:  squeue -u \$USER   (wait until EMPTY before aggregating)"
echo "when done:  bash run_lysisoff.sh agg"
