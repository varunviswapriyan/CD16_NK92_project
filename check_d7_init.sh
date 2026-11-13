#!/usr/bin/env bash
# check_d7_init.sh -- are the day-7 init files in the same format as the day-14 file we use?
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash check_d7_init.sh
G=~/giuseppe_ref
TOP=$G/sex_bias_initialization_files
CD4=$G/src_2024-11-18_post-Zihai_adding_CD4s/sex_bias_initialization_files
D14=~/CD16_NK92_project/SPPARKS_Sex_Bias/rsc_out_2026-07-03/sex_bias_initialization_files/d14_M_4_SPPARKS_init
for f in $TOP/d7_M_1_SPPARKS_init $TOP/d7_F_3_SPPARKS_init $CD4/d7_M_1_SPPARKS_init $CD4/d7_F_3_SPPARKS_init $D14; do
  echo "== $f"
  [ -f "$f" ] || { echo "   (missing)"; continue; }
  awk '/^Values/{v=1;next} v&&NF>5{print "   columns:",NF; exit}' "$f"
  grep -m3 -iE "sites|values|i[0-9]+ " "$f" | sed 's/^/   /'
done
for s in d7_M_1 d7_F_3; do
  cmp -s $TOP/${s}_SPPARKS_init $CD4/${s}_SPPARKS_init && echo "$s: top-level and CD4 copies identical" || echo "$s: top-level and CD4 copies DIFFER"
done
