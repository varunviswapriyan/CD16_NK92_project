#!/usr/bin/env bash
# run_cd4off.sh -- SI control: remove CD4+ helper T cells (SI #3).
#   CD4 helpers are columns i8 (Th, AR-) and i9 (ARTh, AR+) in the init slide.
#   We zero those columns -> no CD4 cells -> no CD4 IFN-g secretion. (CD4's
#   only functional role is IFN-g secretion, so removing them = silencing them
#   for the fold-change question.) All other params = prod/manuscript.
#
# Submit:     cd ~/CD16_NK92_project && git pull && NSEED=40 bash run_cd4off.sh
# Aggregate:  bash run_cd4off.sh agg
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
NSEED="${NSEED:-40}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
ROOT=runs_big/cd4off

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
for cond in ("male_nocd4","female_nocd4"):
    vs=[x for x in (fold(f) for f in sorted(glob.glob(f"{root}/{cond}/s*/log.spparks"))) if x]
    if vs:
        m=sum(vs)/len(vs); sd=math.sqrt(sum((x-m)**2 for x in vs)/(len(vs)-1)) if len(vs)>1 else 0
        print(f"{cond:14s} n={len(vs):3d} mean fold={m:.4f} SD={sd:.4f}")
print("compare standard (CD4 present):  male 2.3793  female 1.5456")
PY
  exit 0
fi

[ -x "$BIN" ] || { echo "build first"; exit 1; }
for cond in male_nocd4 female_nocd4; do
  case $cond in
    male_nocd4)   slide=d14_M_4_SPPARKS_init; fem=0;;
    female_nocd4) slide=d14_F_4_SPPARKS_init; fem=1;;
  esac
  base=$ROOT/$cond
  mkdir -p "$base/sex_bias_initialization_files"
  # make a no-CD4 slide once: zero i8 (field 10) and i9 (field 11) in Values rows
  nocd4=$base/sex_bias_initialization_files/${slide}_noCD4
  before=$(awk 'f&&NF>15{c+=($10!=0)+($11!=0)} /^Values/{f=1} END{print c+0}' "$SRC/sex_bias_initialization_files/$slide")
  awk 'f && NF>15 {$10=0; $11=0} /^Values/{f=1} {print}' \
      "$SRC/sex_bias_initialization_files/$slide" > "$nocd4"
  after=$(awk 'f&&NF>15{c+=($10!=0)+($11!=0)} /^Values/{f=1} END{print c+0}' "$nocd4")
  echo "   $cond: CD4 entries before=$before after=$after (should drop to 0)"
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
cp -f "$PWD/$nocd4" "\$d/sex_bias_initialization_files/${slide}_noCD4"
sed "/^seed/ s/[0-9][0-9]*/\$seed/ ; s|read_sites.*|read_sites sex_bias_initialization_files/${slide}_noCD4|" \\
    "$PWD/$SRC/in.full" > "\$d/in.simulation"
if [ "$fem" = 1 ]; then
  awk '{ if (\$0 ~ /recruitment/ && \$0 !~ /app_style/){ if(\$0 ~ /AR\\+/)\$6=0; else if(\$0 ~ /Tp recruitment/||\$0 ~ /Tm recruitment/)\$6=\$6/0.63 }
         if (\$0 ~ /Tp differentiation/ && \$0 !~ /app_style/){ if(\$0 ~ /AR\\+ Tm/)\$6=0; else \$6=\$6*2 } print }' \\
    "\$d/in.simulation" > "\$d/in.tmp" && mv "\$d/in.tmp" "\$d/in.simulation"
fi
cd "\$d"; mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
  sbatch "$base/array.sh" && echo "   submitted $cond x$NSEED"
done
echo
echo "watch:  squeue -u \$USER"
echo "when done:  bash run_cd4off.sh agg"
