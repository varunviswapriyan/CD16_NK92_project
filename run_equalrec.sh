#!/usr/bin/env bash
# run_equalrec.sh -- control: make AR+ and AR- T-cell recruitment EQUAL.
#   Male:   AR+ Tm recruitment raised 1.52712 -> 6.87204 (= AR- rate).
#   Female: unchanged (females have no AR+ cells; AR+ recruitment stays 0).
#   All other parameters = prod/manuscript.  n = NSEED seeds per condition.
#
# Submit:     cd ~/CD16_NK92_project && git pull && NSEED=40 bash run_equalrec.sh
# Aggregate:  bash run_equalrec.sh agg      (after the jobs finish)
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
NSEED="${NSEED:-40}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
ROOT=runs_big/equalrec

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
S={}
for cond in ("male_eq","female"):
    vs=[x for x in (fold(f) for f in sorted(glob.glob(f"{root}/{cond}/s*/log.spparks"))) if x]
    if vs:
        m=sum(vs)/len(vs)
        sd=math.sqrt(sum((x-m)**2 for x in vs)/(len(vs)-1)) if len(vs)>1 else 0
        se=sd/math.sqrt(len(vs)) if len(vs)>1 else 0
        S[cond]=(m,se)
        print(f"{cond:10s} n={len(vs):3d}  mean fold={m:.4f}  SD={sd:.4f}  SE={se:.4f}")
if "male_eq" in S and "female" in S:
    (mm,ms),(fm,fs)=S["male_eq"],S["female"]
    r=mm/fm; rse=r*math.sqrt((ms/mm)**2+(fs/fm)**2) if mm and fm else 0
    print(f"\nmale(equal-recruit)/female ratio = {r:.4f} +/- {rse:.4f}")
    print("(compare: standard prod ratio = 1.5394)")
PY
  exit 0
fi

# ---- submit mode -------------------------------------------------------------
[ -x "$BIN" ] || { echo "build first"; exit 1; }
for cond in male_eq female; do
  case $cond in
    male_eq) slide=d14_M_4_SPPARKS_init; fem=0; eq=1;;
    female)  slide=d14_F_4_SPPARKS_init; fem=1; eq=0;;
  esac
  base=$ROOT/$cond
  mkdir -p "$base"
  cp -f "$SRC/sex_bias_initialization_files/$slide" "$base/$slide"
  SEDPROG="/^seed/ s/[0-9][0-9]*/\$seed/ ; s|read_sites.*|read_sites sex_bias_initialization_files/$slide|"
  EQ=""
  [ "$eq" = 1 ] && EQ=" ; s/1.52712/6.87204/g"
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
sed "${SEDPROG}${EQ}" "$PWD/$SRC/in.full" > "\$d/in.simulation"
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
echo "watch:  squeue -u \$USER"
echo "when done:  bash run_equalrec.sh agg"
