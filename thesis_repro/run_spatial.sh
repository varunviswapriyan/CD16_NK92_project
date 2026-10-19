#!/usr/bin/env bash
# run_spatial.sh -- re-run male & female (thesis params) WITH spatial dumps,
# so we can draw the cell-position and IFN-g maps at day 14, 17.5, 21.
# Cluster: cd ~/CD16_NK92_project && git pull && cp thesis_repro/run_spatial.sh . && bash run_spatial.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
SRC=rsc_out_2026-07-03
BIN="$PWD/spk_sexbias"
SEED="${SEED:-8377}"
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
[ -x "$BIN" ] || { echo "build first"; exit 1; }

# thesis-era parameters (same as the validated growth figure)
sed -e '/^app_style/ s/10000/1000/' -e 's/6.87204/0.1909/g' -e 's/1.52712/0.0212/g' \
    "$SRC/in.full" > variant_thesis.in

for spec in "male:d14_M_4_SPPARKS_init:male" "female:d14_F_4_SPPARKS_init:female"; do
  IFS=: read -r cond slide sex <<<"$spec"
  d=runs_spatial/$cond
  mkdir -p "$d/sex_bias_initialization_files"
  cp -f "$SRC/sex_bias_initialization_files/$slide" "$d/sex_bias_initialization_files/"
  sed -e "/^seed/ s/[0-9][0-9]*/$SEED/" \
      -e "s|read_sites.*|read_sites sex_bias_initialization_files/$slide|" \
      variant_thesis.in > "$d/in.simulation"
  # turn ON dumps at t = 0, 504, 1008 (day 14, 17.5, 21); fields: coords, cancer(i1),
  # CD8 subsets (i2..i7) and free IFN-g (d1)
  sed -i 's|^#dump .*|dump 1 sites 504.0 dump.2D_model.* id x y i1 i2 i3 i4 i5 i6 i7 d1|' "$d/in.simulation"
  if [ "$sex" = female ]; then
    awk '{
      if ($0 ~ /recruitment/ && $0 !~ /app_style/) {
        if ($0 ~ /AR\+/) $6=0; else if ($0 ~ /Tp recruitment/ || $0 ~ /Tm recruitment/) $6=$6/0.63; }
      if ($0 ~ /Tp differentiation/ && $0 !~ /app_style/) {
        if ($0 ~ /AR\+ Tm/) $6=0; else $6=$6*2; }
      print }' "$d/in.simulation" > "$d/in.tmp" && mv "$d/in.tmp" "$d/in.simulation"
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
echo "when done:  python3 ~/CD16_NK92_project/thesis_repro/plot_spatial.py runs_spatial"
