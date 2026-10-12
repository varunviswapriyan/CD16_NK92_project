# nk92-ca-data

Turn Oscar's exported FCS files into the `Ca_NK92.csv` that `bootstrap_ca2.py` fits.

This is the step Indrani handed off on Oct 8:

> You should be able to calculate the mean and SE of the Ca signal using the
> `Ca_flow_analysis_mean_se.ipynb` notebook located in the `ca_data` folder on the cluster.

---

## Cross-check of her notebook — what I found

Her method is sound and this repo reproduces it to **1e-14** (verified against her exact
`np.digitize` + `groupby.std(ddof=1)` code on Oscar's real files). Three things break if
you run her notebook on the new folders as-is:

**1. Hard-coded file names (will crash).** Cell 1 expects

```
export_liveNK92_344663_NK92_zeta.fcs
```

Oscar's actual exports are

```
export_liveNK92_344663_NK92 (CD3z).fcs          Exp 1
export_20260618NK_361774_NK92 (CD3z).fcs        Exp 2
```

`flowkit.Sample()` raises on the missing path. Fixed by matching the tag in parentheses.

**2. The time scale depends on undocumented FlowKit behaviour.** Her code does
`raw_time * $TIMESTEP` = `raw * 0.01`. On its own that gives **33,133 s (9.2 hours)** for a
5-minute acquisition. It comes out right only because `get_events(source="raw")` already
applied the `$PnG` gain (0.01) inside FlowKit, so the net factor is 1e-4 → **331 s**.

The decisive evidence is `$P7R`, the declared range of the Time parameter: Oscar's files
say **331**, and `raw_max × 1e-4 = 331.3`. Exact match. So her numbers are correct — but
nothing in the notebook checks it, and if FlowKit ever changes that behaviour the time axis
silently becomes 100× too long and every fit is garbage with no error. Both the fixed
notebook and `parse_ca_fcs.py` now assert converted-max ≈ `$PnR` and refuse to run otherwise.

**3. Experiment 2 produces ~80 NaN rows.** She puts 250 bins across the max time over *all
four* tubes. In Exp 2 the CD3z tube stops at 224.8 s and FcRg at 221.3 s, while CD3z+FcRg
runs to 328.8 s. Bins past ~222 s are empty for zeta and gamma → `mean` NaN, `SE` NaN. Those
rows go straight into the CSV and into the fit. Fixed with `TRIM = True` / `--trim`, which
keeps the window all three share (0.7–221.6 s, 169 of 250 bins).

Everything else is unchanged: **250 bins**, `np.linspace(0, global_max, 251)`,
`np.digitize`, `std(ddof=1)`, `SE = std/√n`, and the export columns
`time_seconds, mean_zeta, mean_gamma, mean_hetero, SE_zeta, SE_gamma, SE_hetero`.

## Heads-up before you reply to Indrani

The `20260407NFH` folder Oscar reshared is **Experiment 1 Replicate 1** — the replicate we
already had. Parsed output is byte-for-byte identical to the earlier
`Ca_NK92_exp1rep1_from_fcs.csv` (max abs diff `0.0`).

Indrani asked for **Exp 1 Replicate 2** and **Exp 2 Replicate 1**. We got Exp 1 Rep 1 and
Exp 2 Rep 1. **Exp 1 Rep 2 is still outstanding** — worth one line back to Oscar.

---

## What is in Oscar's two folders

| zip | experiment | cells/tube | usable window |
|---|---|---|---|
| `Exported_FCS_Files_20260407NFH.zip` | Exp 1 (7 Apr) | 34k–79k | 0–331 s |
| `exported_NK_FCS_files_20260518NFH_2.zip` | Exp 2 (18 May) | 120k each | 0–222 s * |

\* truncated by the CD3z/FcRg tubes, see point 3 above.

| Oscar's tube | condition | ITAMs |
|---|---|---|
| `NK92` (or `Parental`) | `parental` | — |
| `CD3z` (or `KI-1 CD3z_CD16`) | `zeta` | 6 |
| `FcRg` (or `KI-2 Fcg_CD16`) | `gamma` | 2 |
| `CD3z+FcRg` (or `KI-6 Fcg_CD3z_CD16`) | `hetero` | 4 |

Both naming conventions are handled.

---

## 0. The push did not work — fix this first

Your `git pull` printed:

```
create mode 160000 nk92-ca-data
```

Mode `160000` is a **gitlink**, not a directory. `nk92-ca-data` had its own `.git` inside it
(I ran `git init` when I built it), so git recorded it as a submodule and pushed a single
pointer instead of the files — same thing that happened with `sexbias_paper_code`. On the
cluster that path is an empty stub.

On your Mac:

```bash
rm -rf ~/Downloads/indrani_code/nk92-ca-data/.git
cd ~/Downloads/indrani_code
git rm --cached nk92-ca-data
git add nk92-ca-data
git commit -m "add Ca FCS pipeline (real files, not a gitlink)"
git push
```

Confirm before pushing: `git ls-files nk92-ca-data | head` should list several files,
not one line.

Then on the cluster:

```bash
cd ~/CD16_NK92_project && git pull
ls nk92-ca-data        # should show ca_go.sh, parse_ca_fcs.py, fcs_zips/, ...
```

## 1. Then run one command

```bash
bash ~/CD16_NK92_project/nk92-ca-data/ca_go.sh
```

That is the whole of the old steps 1-3. It finds the FCS zips, Indrani's notebook, and
every `Ca_NK92.csv` on the system by itself; unzips and parses both experiments; runs the
fixed notebook if flowkit is importable and diffs it against the CLI parser; then runs the
gate — new Exp 1 vs whatever the current fits were actually run on — and prints a verdict.
It writes everything to `~/ca_check/` and changes nothing else.

Two more:

```bash
bash ~/CD16_NK92_project/nk92-ca-data/ca_go.sh serve   # http.server, to see the PNGs
bash ~/CD16_NK92_project/nk92-ca-data/ca_go.sh swap    # back up + install + print the sbatch line
```

`swap` is the only one that writes outside `~/ca_check/`, and it backs the old file up first.
If more than one `Ca_NK92.csv` exists it stops and lists them; pick with
`TARGET=<path> bash ca_go.sh swap`.

---

## Appendix: the manual steps, if you ever need them

### Parse the zips on your Mac (~10 s)

```bash
cd ~/Downloads
bash nk92-ca-data/run_local.sh \
  Exported_FCS_Files_20260407NFH.zip \
  exported_NK_FCS_files_20260518NFH_2.zip
```

Writes `./ca_out/`:

```
Ca_NK92_exp1_20260407.csv        250 rows, drop-in for ca_data/Ca_NK92.csv
Ca_NK92_exp2_20260518.csv        169 rows (trimmed to 221.6 s)
*_full.csv                       adds n_* and the parental tube, for QC
ca_data.png                      mean +/- SE, both experiments
```

Open `ca_data.png` and check against Oscar's slide panels before going further.

Knobs: `NBINS=250`, or call `parse_ca_fcs.py` directly for `--grid fixed --bin 1.0`.

### Onto the cluster — via git, not scp

`scp` to `10.73.170.128` times out; that is an internal address. Push and pull instead.

```bash
# on your Mac
cp -r ~/Downloads/nk92-ca-data ~/Downloads/indrani_code/
mkdir -p ~/Downloads/indrani_code/ca_data_new
cp ~/Downloads/ca_out/Ca_NK92_exp*.csv ~/Downloads/indrani_code/ca_data_new/

cd ~/Downloads/indrani_code
git add nk92-ca-data ca_data_new
git commit -m "Ca FCS -> mean/SE pipeline + parsed CSVs for exp1 and exp2" \
           -m "Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
git push
```

```bash
# on the cluster
ssh vxv016@10.73.170.129
cd ~/CD16_NK92_project && git pull
```

### Cross-check against her notebook (what ca_go.sh automates)

`Ca_flow_analysis_mean_se_fixed.ipynb` in this repo is her notebook with the three fixes.
Run it on an experiment folder, then let its last cell diff itself against the CLI parser.

```bash
cd ~/CD16_NK92_project/nk92-ca-data
cp Ca_flow_analysis_mean_se_fixed.ipynb <folder with the 4 .fcs files>/
cd <that folder>
jupyter nbconvert --to notebook --execute Ca_flow_analysis_mean_se_fixed.ipynb \
        --output executed.ipynb --ExecutePreprocessor.timeout=1200
```

Cell 2 prints the time-scale check and raises if the gain assumption ever breaks. Cell 6
runs `parse_ca_fcs.py` on the same folder and `compare_ca_csv.py` on the two outputs;
expect `IDENTICAL`.

Then the gate that matters — diff the new CSV against the one the current fits were run on:

```bash
cd /home/gddaslab/vxv016/Ca_fit_c02
python3 ~/CD16_NK92_project/nk92-ca-data/compare_ca_csv.py \
        ca_data/Ca_NK92.csv \
        ~/CD16_NK92_project/ca_data_new/Ca_NK92_exp1_20260407.csv
```

**If that prints `DIFFERENT`, every SSR in the draft email (orig = 5.206e4 and friends) is
stale.** Re-run the fits before that email goes out.

### Point the fits at the data and submit

```bash
cd /home/gddaslab/vxv016/Ca_fit_c02
cp ca_data/Ca_NK92.csv ca_data/Ca_NK92.csv.bak.$(date +%F)
cp ~/CD16_NK92_project/ca_data_new/Ca_NK92_exp1_20260407.csv ca_data/Ca_NK92.csv

bash ~/CD16_NK92_project/filesCC/finish.sh ca
squeue -u vxv016
```

Use `finish.sh ca`. **Not** `finish.sh submit` — that one `rmtree`s the 60 finished grid
cells and you lose the ZAP0 × kzp map.

### Look at the results

```bash
cd ~/nk92-fits
bash show_results.sh
python3 plot_ca_fit.py orig
```

To see the PNGs, serve your home directory and open it in the cluster VM's Chrome:

```bash
cd ~ && python -m http.server 8000 &
```

then, depending on which login node you are on:

```
http://10.73.170.128:8000/...      # log01
http://10.73.170.129:8000/...      # log02
```

e.g. `http://10.73.170.129:8000/CD16_NK92_project/nk92-ca-data/ca_data.png`.

---

## Files

| file | what it does |
|---|---|
| `parse_ca_fcs.py` | FCS 3.0 reader → binned mean/SE → CSV. No flowkit. Asserts the time scale. |
| `Ca_flow_analysis_mean_se_fixed.ipynb` | Indrani's notebook, three fixes, self-cross-checking |
| `ca_go.sh` | **the one command** — discover, parse, cross-check, gate; `swap` / `serve` |
| `run_local.sh` | unzip + parse both experiments + plot, one command |
| `compare_ca_csv.py` | diff two Ca CSVs; exit 1 if they differ |
| `plot_ca_data.py` | mean ± SE traces, to check against the slides |

### How the parsing works

Oscar's files are FCS 3.0, list mode, `$DATATYPE F` (float32), `$BYTEORD 4,3,2,1`
(big-endian), 11 parameters, with the FlowJo gate already applied — every event in the file
is a live NK92, so no gating is needed here.

`parse_ca_fcs.py` reads the header directly, so it has no flowkit dependency and the gain
handling is explicit rather than inherited from a library default. Binning uses
`np.bincount` on running sums of `r` and `r²`; a 120k-event file takes well under a second,
and the result matches the notebook's `groupby` to 1e-14.

### Open items this does not settle

* Exp 1 Replicate 2 was never delivered — ask Oscar.
* How to combine replicates (average the means? fit each separately and compare parameters?)
  is a modeling call for Indrani/Das.
* Whether to fit Exp 2 on the trimmed 0–222 s window or ask Oscar why those two tubes stop
  early.
