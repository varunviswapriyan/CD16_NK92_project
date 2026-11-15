#!/usr/bin/env bash
# zihai_overnight.sh -- everything with Giuseppe's Zenodo mIF data (Processed_Zihai_Data.zip).
#
#   tonight:   cd ~/CD16_NK92_project && git pull --no-edit && bash zihai_overnight.sh
#   morning:   cd ~/CD16_NK92_project && bash zihai_overnight.sh analyze
#
# Tonight (a few minutes, then jobs run overnight):
#   1. unzips the data (looks for Processed_Zihai_Data.zip in SPPARKS_Sex_Bias/, the repo, ~ or ~/Downloads)
#   2. CHECKS the converter: rebuilds Giuseppe's d14_M_4 and d14_F_4 starting files from the raw slides
#      and compares column by column. Nothing is submitted unless they match (within 3%).
#   3. builds day-7 starting files WITH CD4 for the same slides/windows Giuseppe used (d7_M_1, d7_F_3)
#   4. cancels the old no-CD4 day-7 jobs (d7_m, d7_f) and submits day 7 -> 14 runs, n=20 per sex
#   5. per-slide numbers from all 29 slides (AR+ fraction, CD4/CD8, densities, infiltration profiles)
# Morning: model-vs-data figures + day 7->14 fold change and male/female ratio (zihai_summary.txt).
set -uo pipefail
P=~/CD16_NK92_project
cd $P/SPPARKS_Sex_Bias || exit 1
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
export PYTHONPATH=$P
SRC=rsc_out_2026-07-03; DECK="$PWD/$SRC/in.full"; BIN="$PWD/spk_sexbias"; JOB=$PWD/runs_big/_job_generic.sh
ND7="${ND7:-20}"; TIME="${TIME:-30:00:00}"
DATA=$PWD/zihai_data/Processed_Zihai_Data
REF14=$SRC/sex_bias_initialization_files
REF7=~/giuseppe_ref/sex_bias_initialization_files
OUTI=$PWD/zihai_init
mkdir -p figures $OUTI

if [ "${1:-}" = "analyze" ]; then
  echo "jobs still queued/running: $(squeue -u $USER -h -n d7c_m,d7c_f 2>/dev/null | wc -l)"
  python3 $P/zihai_analysis.py data "$DATA" > /dev/null && python3 $P/zihai_analysis.py figures
  exit 0
fi

# ---------- 1. data ----------
if [ ! -d "$DATA" ]; then
  Z=$(ls -t $PWD/Processed_Zihai_Data.zip $P/Processed_Zihai_Data.zip ~/Processed_Zihai_Data.zip ~/Downloads/Processed_Zihai_Data.zip 2>/dev/null | head -1)
  [ -n "$Z" ] || { echo "!! Processed_Zihai_Data.zip not found. Put it in ~/CD16_NK92_project/SPPARKS_Sex_Bias/ and re-run."; exit 1; }
  echo "unzipping $Z"
  [ "$(md5sum "$Z" | cut -d' ' -f1)" = "287e9c3df397de316db6583bebd896ad" ] && echo "  checksum matches Zenodo" || echo "  (checksum differs from Zenodo's; continuing, the file check below will catch a bad copy)"
  mkdir -p zihai_data && unzip -q -o "$Z" -d zihai_data
fi
n=$(ls $DATA/d*_*.csv 2>/dev/null | wc -l); echo "slides found: $n (expected 29)"; [ "$n" -ge 20 ] || exit 1

# ---------- 2+3. check converter, build day-7 files ----------
python3 - "$DATA" "$REF14" "$REF7" "$OUTI" <<'PY'
import sys, os
from convert_zihai import read_slide, read_init, find_window, build, write, compare
data, ref14, ref7, outi = sys.argv[1:]
ok = True
print("\n=== CHECK: rebuild Giuseppe's day-14 starting files from the raw slides ===")
for s in ("d14_M_4", "d14_F_4"):
    ref = read_init(f"{ref14}/{s}_SPPARKS_init"); sl = read_slide(f"{data}/{s}.csv")
    c0, r0, fx, fy = find_window(sl, ref); print(f"-- {s}: window at ({c0*10}, {r0*10}) um")
    r = compare(build(sl, c0, r0, fx, fy), ref)
    good = abs(r["cd8_total"]) < .03 and abs(r["cd4_total"]) < .03 and abs(r["cancer_total"]) < .03 and r["cd8_site_match"] > .85
    print("   ->", "MATCHES" if good else "DOES NOT MATCH"); ok &= good
print("\n=== BUILD day-7 starting files with CD4 (same windows as Giuseppe's day-7 files) ===")
for s in ("d7_M_1", "d7_F_3"):
    sl = read_slide(f"{data}/{s}.csv"); rp = f"{ref7}/{s}_SPPARKS_init"
    ref = read_init(rp) if os.path.exists(rp) else None
    c0, r0, fx, fy = find_window(sl, ref)
    out = build(sl, c0, r0, fx, fy)
    if ref is not None:
        print(f"-- {s}: window at ({c0*10}, {r0*10}) um, compared with Giuseppe's old (no-CD4) file:")
        r = compare(out, ref)
        if r["cd8_site_match"] < .6:
            c0, r0, fx, fy = find_window(sl, None); out = build(sl, c0, r0, fx, fy)
            print(f"   !! poor match to his window; using the 6x6 mm window with most cells instead ({c0*10}, {r0*10}) um")
    else:
        print(f"-- {s}: no old file found; using the 6x6 mm window with most cells ({c0*10}, {r0*10}) um")
    print(f"   cells: cancer {out['i1'].sum():.0f}, AR- CD8 {sum(out[k].sum() for k in ('i2','i3','i4')):.0f}, "
          f"AR+ CD8 {out['i6'].sum()+out['i7'].sum():.0f}, CD4 {out['i13'].sum()+out['i14'].sum():.0f}")
    write(out, f"{outi}/{s}_CD4_SPPARKS_init"); print(f"   wrote zihai_init/{s}_CD4_SPPARKS_init")
sys.exit(0 if ok else 2)
PY
rc=$?
if [ $rc -ne 0 ]; then echo; echo "!! Converter check failed: NOT submitting runs. Send Claude this output."; exit 1; fi

# ---------- 4. day 7 -> 14 runs ----------
echo; echo "=== day 7 -> 14 runs (n=$ND7 per sex) ==="
[ -x "$BIN" ] && [ -f "$JOB" ] || { echo "missing $BIN or $JOB"; exit 1; }
scancel -n d7_m,d7_f 2>/dev/null && echo "  cancelled old no-CD4 day-7 jobs (if any)"
submit(){ # name base n slide fem
  local name=$1 base=$2 n=$3 sl=$4 fem=$5
  [ -f $base/cfg.sh ] && { echo "  $base already submitted"; return; }
  mkdir -p $base; : > $base/extra.sed
  cat > $base/cfg.sh <<EOF
BASE=$PWD/$base
SEEDBASE=9000
SLIDEPATH=$(readlink -f $sl)
SLIDENAME=$(basename $sl)
FEM=$fem
DECK=$DECK
BIN=$BIN
EXTRASED=$PWD/$base/extra.sed
EOF
  sbatch --job-name=$name --array=1-$n --time=$TIME --ntasks=1 --output=$PWD/$base/slurm_%A_%a.out \
         --export=ALL,CFG=$PWD/$base/cfg.sh "$JOB" >/dev/null && echo "  submitted $name x$n -> $base"
}
submit d7c_m runs_d7cd4/male   $ND7 $OUTI/d7_M_1_CD4_SPPARKS_init 0
submit d7c_f runs_d7cd4/female $ND7 $OUTI/d7_F_3_CD4_SPPARKS_init 1

# ---------- 5. slide numbers ----------
echo; echo "=== per-slide numbers from all slides ==="
python3 $P/zihai_analysis.py data "$DATA"
echo; echo "Done for tonight. Tomorrow: cd ~/CD16_NK92_project && bash zihai_overnight.sh analyze"
