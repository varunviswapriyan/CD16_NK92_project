#!/bin/bash
# setup_batch75.sh -- v75..v90: fit pZAP WITHOUT the gamma-affinity fraction.
#
# Every variant sets KZBG_FRAC = 1.0 (FIXED, not fitted), per Indrani's point #3.
# The zeta/gamma amplitude ordering must then come from something else. Four
# mechanisms are tested, all evidence-based or already-fitted quantities:
#
#   A) nothing extra (controls: does the fit survive FRAC removal at all?)
#   B) zeta-3 ITAM lower ZAP70 affinity -- MEASURED, Bu/Shaw/Chan 1995 (her
#      Table S1 ref 7): zeta3 binds ~2.5x less well than zeta1/zeta2/gamma.
#      Implemented as kzb_z3 = kzb * Z3_FRAC on ITAM3 rules (zeta + hetero).
#   C) limited ZAP70 / higher SYK competition (ZAP0, SYK0 already fitted --
#      just widened bounds; 6 zeta ITAMs compete harder than 2 gamma ITAMs)
#   D) other existing rate multipliers (kzu, kdl) added to the fitted set
#
# Warm start = bounds narrowed around v69's optimum, so 30 iterations suffice
# (v69 needed 90 from cold). No change to her fitting script.
#
# Usage:
#   bash ~/CD16_NK92_project/v75_batch/setup_batch75.sh          # build only
#   bash ~/CD16_NK92_project/v75_batch/setup_batch75.sh submit   # build+submit
#   bash ~/CD16_NK92_project/v75_batch/setup_batch75.sh report   # summary

set -e
SRC="$HOME/NK92_fit_v69"
SUB="estimate_params_pzap_cleaned_up"
HOME_DIR="$HOME"
PARTICLES=24
ITERS=30
MODE="${1:-build}"

VARIANTS="75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96"

if [ "$MODE" = "report" ]; then
  printf '%-5s %-34s %12s\n' V description cost
  printf -- '---------------------------------------------------------------\n'
  for V in $VARIANTS; do
    D="$HOME_DIR/NK92_fit_v$V/$SUB"
    NOTE=$(python3 -c "import json;print(json.load(open('$D/v_config.json'))['NOTE'])" 2>/dev/null || echo '?')
    C=$(grep -h "best cost" $HOME_DIR/NK92_fit_v$V/*.log $D/slurm*.out 2>/dev/null | tail -1 | sed 's/.*best cost: *//' | cut -c1-10)
    R=$(grep -h "residue \*\*" $D/slurm*.out $HOME_DIR/NK92_fit_v$V/*.log 2>/dev/null | tail -1 | sed 's/.*= *//' | cut -c1-8)
    printf '%-5s %-34s %12s %s\n' "v$V" "${NOTE:0:34}" "${C:--}" "${R:+res=$R}"
  done
  printf -- '---------------------------------------------------------------\n'
  echo "v69 benchmark: worst single-point miss 0.043 (WITH the affinity fraction)"
  exit 0
fi

[ -d "$SRC/$SUB" ] || { echo "ERROR: $SRC/$SUB not found"; exit 1; }

# ---- v69 optimum (log10), order: lig0 kdl0 ZAP0 SYK0 KZBG_FRAC KZP_MULT KPR_MULT
BEST="1.97853984 -2.54048026 2.19862004 1.71998766 -0.25851195 0.34152821 -0.26705681"

build_variant () {   # $1=V  $2=note  $3=python-config-edits  $4=z3flag
  local V="$1" NOTE="$2" EDITS="$3" Z3="$4"
  local D="$HOME_DIR/NK92_fit_v$V"
  rm -rf "$D"; cp -r "$SRC" "$D"
  rm -rf "$D/$SUB"/gamma_runs "$D/$SUB"/zeta_runs "$D/$SUB"/mixed_runs \
         "$D/$SUB"/*.png "$D/$SUB"/slurm*.out 2>/dev/null || true

  # --- zeta-3 ITAM affinity patch (zeta file + hetero/mixed file) ---
  if [ "$Z3" = "yes" ]; then
    python3 - "$D/$SUB" <<'PY'
import sys, re, glob, os
d = sys.argv[1]
files = [f for f in glob.glob(os.path.join(d, 'JJ_*.bngl'))
         if 'gamma' not in os.path.basename(f)]   # zeta + mixed only
total = 0
for f in files:
    L = open(f).read().split('\n')
    hits = [i for i, l in enumerate(L) if re.match(r'^\s*kzb\s', l)]
    if len(hits) != 1:
        sys.exit(f'ERROR {f}: expected 1 kzb param line, found {len(hits)}')
    if not any(l.startswith('kzb_z3') for l in L):
        L.insert(hits[0] + 1, 'kzb_z3 kzb * Z3_FRAC')
    n = 0
    for i, l in enumerate(L):
        # zeta-derived ITAM3 binding rules only
        if re.search(r'\bITAM3\b', l) and 'binding' in l.lower():
            for j in range(i, min(i + 4, len(L))):
                if re.search(r'\bkzb\s*,', L[j]):
                    L[j] = re.sub(r'\bkzb\s*,', 'kzb_z3,', L[j], count=1); n += 1; break
    if n == 0:
        sys.exit(f'ERROR {f}: no ITAM3 binding rules patched')
    open(f, 'w').write('\n'.join(L))
    print(f'    {os.path.basename(f)}: patched {n} ITAM3 rule(s)')
    total += n
if total == 0: sys.exit('ERROR: no rules patched anywhere')
PY
    # ensure Z3_FRAC parameter exists in every bngl that uses it
    python3 - "$D/$SUB" <<'PY'
import sys, glob, os, re
for f in glob.glob(os.path.join(sys.argv[1], 'JJ_*.bngl')):
    s = open(f).read()
    if 'kzb_z3' in s and not re.search(r'^\s*Z3_FRAC\s', s, re.M):
        L = s.split('\n')
        i = next(i for i, l in enumerate(L) if re.match(r'^\s*KZBG_FRAC\s', l))
        L.insert(i, 'Z3_FRAC 1.0')
        open(f, 'w').write('\n'.join(L))
PY
  fi

  # --- config ---
  BEST="$BEST" NOTE="$NOTE" EDITS="$EDITS" python3 - "$D/$SUB/v_config.json" <<'PY'
import json, os, sys
p = sys.argv[1]
c = json.load(open(p))
best = [float(x) for x in os.environ['BEST'].split()]
names = ["lig0","kdl0","ZAP0","SYK0","KZBG_FRAC","KZP_MULT","KPR_MULT"]
opt = dict(zip(names, best))
orig_lb = dict(zip(c["PARAMS"], c["LB"]))
orig_ub = dict(zip(c["PARAMS"], c["UB"]))

def drop(n):
    if n in c["PARAMS"]:
        i = c["PARAMS"].index(n)
        for k in ("PARAMS","LB","UB"): c[k].pop(i)

def add(n, lb, ub):
    if n not in c["PARAMS"]:
        c["PARAMS"].append(n); c["LB"].append(lb); c["UB"].append(ub)

def setb(n, lb, ub):
    if n in c["PARAMS"]:
        i = c["PARAMS"].index(n); c["LB"][i] = lb; c["UB"][i] = ub

def warm(half=0.25):
    for n in list(c["PARAMS"]):
        if n in opt:
            i = c["PARAMS"].index(n)
            lo = max(orig_lb.get(n, -9), opt[n] - half)
            hi = min(orig_ub.get(n, 9), opt[n] + half)
            c["LB"][i], c["UB"][i] = round(lo, 4), round(hi, 4)

# always: remove the gamma affinity fraction from fitting, hold it at 1.0
drop("KZBG_FRAC")
c.setdefault("FIXED", {})["KZBG_FRAC"] = 1.0

exec(os.environ["EDITS"])           # per-variant edits
c["NOTE"] = os.environ["NOTE"]
json.dump(c, open(p, "w"), indent=2)
print(f'    params: {c["PARAMS"]}')
print(f'    fixed : {c["FIXED"]}')
PY

  cat > "$HOME_DIR/run_v$V.sh" <<EOF
#!/bin/bash
#SBATCH --job-name=v$V
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --output=$D/v$V.log
module load Miniconda3/4.9.2
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2
cd $D/$SUB
python pzap_param_estimation_NK92.py $PARTICLES $ITERS
EOF
  chmod +x "$HOME_DIR/run_v$V.sh"
  echo "  v$V built: $NOTE"
}

echo "Building v75..v90 (KZBG_FRAC removed everywhere)..."

# ---- A: controls, no extra mechanism -------------------------------------
build_variant 75 "FRAC=1, warm, kzp+KPR fitted"          'warm()'                                        no
build_variant 76 "FRAC=1, cold (full bounds)"            'pass'                                          no
build_variant 77 "FRAC=1, warm, kzp+KPR FIXED at approved" '
warm()
for n,v in (("KZP_MULT",2.2),("KPR_MULT",0.54)):
    drop(n); c["FIXED"][n]=v'                                                                            no
build_variant 78 "FRAC=1, warm, replicate seed"          'warm()'                                        no

# ---- C: competition for limited ZAP70 / more SYK -------------------------
build_variant 79 "FRAC=1, low-ZAP0 competition"          'warm(); setb("ZAP0", 1.0, 2.6)'                no
build_variant 80 "FRAC=1, high-SYK0 competition"         'warm(); setb("SYK0", 1.5, 3.5)'                no
build_variant 81 "FRAC=1, low ZAP0 + high SYK0"          'warm(); setb("ZAP0",1.0,2.6); setb("SYK0",1.5,3.5)' no

# ---- D: other existing multipliers into the fitted set -------------------
build_variant 82 "FRAC=1, + kzu fitted"                  'warm(); add("KZU_MULT",-1.0,1.0)'              no
build_variant 83 "FRAC=1, + kdl fitted"                  'warm(); add("KDL_MULT",-1.0,1.0)'              no
build_variant 84 "FRAC=1, + kzu + kdl fitted"            'warm(); add("KZU_MULT",-1.0,1.0); add("KDL_MULT",-1.0,1.0)' no

# ---- B: zeta-3 ITAM affinity (Bu et al. 1995) ---------------------------
build_variant 85 "z3 kzb/2.5 FIXED (Bu 1995), warm"      'warm(); c["FIXED"]["Z3_FRAC"]=0.4'             yes
build_variant 86 "z3 kzb/2.5 FIXED, cold"                'c["FIXED"]["Z3_FRAC"]=0.4'                     yes
build_variant 87 "z3 FRAC fitted, warm"                  'warm(); add("Z3_FRAC",-1.0,0.0)'               yes
build_variant 88 "z3 fitted + kzu fitted"                'warm(); add("Z3_FRAC",-1.0,0.0); add("KZU_MULT",-1.0,1.0)' yes
build_variant 89 "z3 fixed 0.4 + low ZAP0"               'warm(); c["FIXED"]["Z3_FRAC"]=0.4; setb("ZAP0",1.0,2.6)' yes
build_variant 90 "z3 fitted, ratio-weighted cost"        '
warm(); add("Z3_FRAC",-1.0,0.0)
c["RATIO_WEIGHT"]=1.0'                                                                                   yes

# ---- E: partial phosphorylation / Lck limitation -------------------------
#   ZAP70 needs DOUBLY phosphorylated ITAMs. Lck is fixed and shared, so 6
#   zeta ITAMs per receptor are phosphorylated less completely than 2 gamma
#   ITAMs -> fewer ZAP70-competent sites. Pure ITAM-number consequence.
build_variant 91 "FRAC=1, + kp fitted (Lck phos rate)"   'warm(); add("KP_MULT",-1.0,1.0)'               no
build_variant 92 "FRAC=1, + kzd fitted (dephos)"         'warm(); add("KZD_MULT",-1.0,1.0)'              no
build_variant 93 "FRAC=1, + kp + kzd (phospho balance)"  'warm(); add("KP_MULT",-1.0,1.0); add("KZD_MULT",-1.0,1.0)' no
build_variant 94 "z3 fixed 0.4 + kp + kzd"               'warm(); c["FIXED"]["Z3_FRAC"]=0.4; add("KP_MULT",-1.0,1.0); add("KZD_MULT",-1.0,1.0)' yes

# ---- F: objective variants + kitchen sink --------------------------------
build_variant 95 "FRAC=1, order-weighted cost"           '
warm(); add("KP_MULT",-1.0,1.0)
c["ORDER_WEIGHT"]=1.0; c["RATIO_WEIGHT"]=0.5'                                                            no
#   everything open at once: NOT for presenting -- if even this cannot fit,
#   that is hard evidence no parameter set works without an affinity difference
build_variant 96 "KITCHEN SINK: z3+kzu+kdl+kp+kzd"       '
warm(0.4)
add("Z3_FRAC",-1.0,0.0); add("KZU_MULT",-1.0,1.0); add("KDL_MULT",-1.0,1.0)
add("KP_MULT",-1.0,1.0); add("KZD_MULT",-1.0,1.0)'                                                       yes

echo
echo "Checking for phosphatase / Lck parameters available for future variants:"
grep -iE "^\s*(LCK|SHP|ksb|ksu|kd?ph)" "$SRC/$SUB"/JJ_zeta*.bngl | head -8 || echo "  (none matched -- SHP-1/Lck are concentrations set elsewhere)"

echo
if [ "$MODE" = "submit" ]; then
  for V in $VARIANTS; do sbatch "$HOME_DIR/run_v$V.sh"; done
  echo "Submitted. squeue -u \$USER"
else
  echo "Built but NOT submitted. Smoke-test one zeta3 variant first (2 min):"
  echo "  cd $HOME_DIR/NK92_fit_v85/$SUB && python pzap_param_estimation_NK92.py 4 1"
  echo "Then: bash \$0 submit"
fi
echo "Later: bash \$0 report"
