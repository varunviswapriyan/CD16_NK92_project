#!/usr/bin/env bash
# show_init_builder.sh -- print only the relevant parts of init_builder_report.txt
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash show_init_builder.sh
R=~/CD16_NK92_project/init_builder_report.txt
[ -f "$R" ] || { echo "no $R -- run find_init_builder.sh first"; exit 1; }
NOISE="melanoIMC|Aaron|_100_samples|dump|\.git/"
echo "################ SECTION 2: what each builder script reads ################"
sed -n '/=== 2\./,/=== 3\./p' "$R" | head -60
echo
echo "################ SECTION 3: mouse-looking file names (noise removed) ################"
sed -n '/=== 3\./,/=== 4\./p' "$R" | grep -v -i -E "$NOISE" | head -40
echo
echo "################ SECTION 4: large data files (noise removed) ################"
sed -n '/=== 4\./,$p' "$R" | grep -v -i -E "$NOISE" | head -40
echo
echo "################ converter script: inputs and markers ################"
F=~/rsgxg003_mount/enableMedicine_Data_Analysis/EnableMedicine_Codex_Analysis_Code/Data_To_SPPARKS_Initialization_File.py
if [ ! -f "$F" ]; then
  mountpoint -q ~/rsgxg003_mount 2>/dev/null || squashfuse /gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh ~/rsgxg003_mount
fi
grep -n -E "read_csv|read_excel|open\(|glob|listdir|_path|_dir|\.csv|\.txt|\.xlsx|CD8|CD3|TIM3|TCF1|AR|0\.7" "$F" 2>/dev/null | head -50
