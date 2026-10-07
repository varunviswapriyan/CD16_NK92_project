# Sex-bias tumour model: manuscript figures

This folder has the code that makes the four figures agreed in the meeting:

| Figure | Shows |
|---|---|
| `Fig_cancer_timecourse.png` | Cancer fold change, male vs female (mean ± SD over seeds), plus male without AR+ internalization. Second panel: total free IFN-γ. |
| `Fig_spatial_cells.png` | Cell positions at day 14, 17.5 and 21, male and female. |
| `Fig_spatial_ifng.png` | Free IFN-γ maps at the same times. |
| `Fig_effective_diffusion.png` | Remake of thesis Fig 4.8: IFN-γ profiles, front² vs time, and D_eff against D/(1+R₀/K_D). |

## Files

- `sex_bias_sim_v7.cpp` is the simulator. It is your v6 port with the `in.full` model, with these changes:
  - the compile bug is fixed;
  - the tuning options are removed;
  - `--no-sequestration` is added;
  - time is given in days;
  - it runs about 3× faster.

  On seed 8377 it gave the same trajectory as v6, number for number.
- `ifng_effective_diffusion.cpp` runs the simplified model from thesis §4.3: one fixed IFN-γ source and one cancer cell per site. It uses the same IFN-γ physics as the simulator.
- `make_figures.py` draws the figures.
- `run_all.sh` builds everything, runs all simulations and draws the figures.

## Setup (once)

```bash
mkdir -p init
cp /path/to/rsc_out_2026-07-03/sex_bias_initialization_files/* init/
# needs: g++ (or clang++), python3 with numpy + matplotlib
```

## Run everything

```bash
./run_all.sh                         # 5 seeds x 3 conditions, 4 at a time
SEEDS="8377 8378 8379" JOBS=3 ./run_all.sh
```

Results go in `results/`:
- `results/runs/` holds the time courses and snapshots;
- `results/figures/` holds the PNGs, plus `summary_day21.txt` and `effective_diffusion.txt`, which have the numbers for the text.

If you re-run, finished runs are skipped.

**Run time.** Each 7-day run takes about 2–3 h on one core (roughly 10 s of CPU per 10-minute unit). Set `JOBS` to your number of physical cores. The diffusion model takes about 15 s.

**Female condition.** By default, "female" means the `d14_M_4_madeF` slide (the male tumour with its AR+ cells relabelled AR−) with the female rates. That isolates the effect of AR. To use the real female slide instead:

```bash
FEMALE_INIT=init/d14_F_4_SPPARKS_init ./run_all.sh
```

That slide starts with more cancer cells, so compare fold changes rather than raw counts.

## Running one simulation by hand

```bash
g++ -O3 -std=c++17 sex_bias_sim_v7.cpp -o sex_bias_sim_v7
./sex_bias_sim_v7 --init init/d14_M_4_SPPARKS_init --sex male --days 7 \
    --seed 8377 --report-hours 6 --snapshot-days 0,3.5,7 --out male_s8377
```

| Option | Meaning |
|---|---|
| `--days D` | Length of the run. 7 days = day 14 → 21 = 1008 native units. |
| `--sex female` | Applies the female rates: no AR+ recruitment, no Tp→AR+ PE, AR− recruitment ×1/0.63, Tp→Tm ×2. |
| `--no-sequestration` | AR+ cells keep their receptors bound instead of internalizing IFN-γ. Nothing else changes. |
| `--te-death-weekly` | Exhausted-cell death of 1/week (rates table) instead of 1/day (`in.full`). |

**Outputs.**
- `PREFIX_timecourse.csv` has one row per report time, with all cell counts, free and bound IFN-γ, and mean bound fractions.
- `PREFIX_snap_tNNNN.txt` holds a spatial snapshot. NNNN is the time in native units.

## Notes

- **Time units.** One native unit is 10 minutes, so 144 units is one day. Older runs labelled "168h" were really 168 units, about 28 hours.
- **The model is `in.full` unchanged.** Nothing is tuned to get a particular result. Whatever the figures show is what this model predicts.
- **`D_eff` calibration.** The diffusion figure gets `D_eff` by comparing the front's x² slope against a bare-diffusion run with the same source. That avoids guessing the geometric prefactor. The analytic estimate D/(1+R₀/K_D) is shown next to it.
- **Open questions for Giuseppe and Jayajit.** These don't change the code above:
  - Exhausted-cell death is 1/day in `in.full` but 1/week in the rates table.
  - AR+ progenitor-exhausted cells kill cancer at the full mature rate.
