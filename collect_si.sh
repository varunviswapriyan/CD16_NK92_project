#!/usr/bin/env bash
# collect_si.sh -- copy all SI figures into SPPARKS_Sex_Bias/SI_figures/ and zip them,
# remaking any fold-change plot that is missing. Download the zip from the browser.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash collect_si.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
F=figures
[ -f $F/Fold_no_ARp_internalization.png ] && [ -f $F/Fold_CD4_removed.png ] || bash ~/CD16_NK92_project/plot_ifngr_fold.sh >/dev/null
[ -f $F/Fold_feedback_off.png ] || bash ~/CD16_NK92_project/plot_fboff_fold.sh >/dev/null
mkdir -p SI_figures && rm -f SI_figures/*
n=0
for f in Fold_no_ARp_internalization.png Fold_CD4_removed.png Fold_feedback_off.png \
         Diffusion_male_vs_female.png Kinetics_free_IFNg_by_condition.png ARp_fraction_CD8.png \
         CD8_subtypes_standard.png Kinetics_standard.png Kinetics_ARp_internalization_off.png \
         Kinetics_IFNGR_blocked_all_T_cells.png Kinetics_CD4_helpers_removed_corrected.png; do
  if [ -f $F/$f ]; then cp $F/$f SI_figures/; n=$((n+1)); echo "  ok      $f"; else echo "  MISSING $f"; fi
done
rm -f SI_figures.zip && zip -qj SI_figures.zip SI_figures/*.png
echo "$n figures -> SPPARKS_Sex_Bias/SI_figures/  and  SI_figures.zip"
echo "Download: python -m http.server 8001, then open http://<node-ip>:8001/SI_figures.zip"
