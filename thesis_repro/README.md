# Reproducing Giuseppe's thesis tumour-growth figure

This runs the **with-CD4s thesis model** — the branch that still has the CD4
helper T cells and Tregs your current `main` code deleted. That deletion is
why your figures don't match his shape.

Put this `thesis_repro/` folder at the **root of the repo**, on the thesis
branch:

```bash
cd ~/SPPARKS_Sex_Bias
git checkout Sex_bias_post-Zihai_with_CD4s
# copy the four files of this kit into ~/SPPARKS_Sex_Bias/thesis_repro/
chmod +x thesis_repro/*.sh
```

## 1. Build (loads MPI, compiles all sources)

```bash
module spider openmpi        # find your cluster's module names
module load gcc openmpi      # adjust to what's available
bash thesis_repro/build.sh   # -> ./spk_sexbias
```

## 2. Smoke test first (≈1 min) — confirms it builds, reads input, emits stats

```bash
SMOKE=1 bash thesis_repro/run_thesis.sh
sed -n '1,40p' thesis_runs/male_s8377/log.spparks
```

Look for a line starting with `Time` followed by numeric rows. Send me the
top of that log so I can lock the cancer column if `i1` isn't the label.

## 3. Full runs (male vs female, + male no-sequestration)

```bash
JOBS=4 bash thesis_repro/run_thesis.sh           # 5 seeds x 3 conditions, run 1008
python3 thesis_repro/parse_and_plot.py thesis_runs
```

Output: `thesis_runs/Fig_cancer_timecourse.png` and a printed day-21 summary
with the male/female fold ratio (his ~1.29).

### What each condition is
- **male** — `d14_M_4` slide, `in.full` unchanged.
- **female** — `d14_M_4_madeF` slide (AR+ relabelled AR−, isolates AR) with his
  exact female rate edits (factor 0.9×0.7 = 0.63). For the **real** female
  slide instead, run `CONDS="male female_real" ...`.
- **male_noseq** — male with AR+ IFN-γ uptake (sequestration) set to 0.

### Cluster (sbatch) instead of local parallel
Giuseppe's `auto_run.sh` isn't in the repo, but each run is just
`mpirun -np 1 ./spk_sexbias < in.simulation > log.spparks`. Wrap the body of
`run_one` in a one-line sbatch script if you want them queued instead of run
with `xargs -P`.
