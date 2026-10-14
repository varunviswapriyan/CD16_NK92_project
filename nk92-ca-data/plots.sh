#!/usr/bin/env bash
# plots.sh -- build every figure that exists right now, put them on one page,
#             serve it, and print ONE url.
#
#   bash plots.sh
#
# Anything whose bootstrap has not finished is skipped with a note, not an error.
set -uo pipefail

PY="${PY:-python3}"
OUT="${OUT:-$HOME/ca_plots}"
PORT="${PORT:-8000}"
F="${F:-$HOME/CD16_NK92_project/filesCC}"
FITS="${FITS:-$HOME/nk92-fits}"
EXP2="${EXP2:-$HOME/Ca_fit_exp2}"
CHECK="${CHECK:-$HOME/ca_check}"
BOOT="$F/bootstrap_ca2.py"

ok()   { printf '   \033[32m+\033[0m %s\n' "$*"; }
skip() { printf '   \033[33m-\033[0m %s\n' "$*"; }

# ------------------------------------------------------------------ python env
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
[ -f "$CONDA_SH" ] && { . "$CONDA_SH" >/dev/null 2>&1; conda activate CD16_v2 >/dev/null 2>&1; }
if [ -z "${PZ_PY:-}" ]; then
  for c in "${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}" python python3; do
    [ -n "$c" ] || continue
    q=$(command -v "$c" 2>/dev/null || { [ -x "$c" ] && echo "$c"; })
    [ -n "$q" ] || continue
    "$q" -c "import numpy,pandas,matplotlib" >/dev/null 2>&1 && { PY="$q"; break; }
  done
else
  PY="$PZ_PY"
fi

rm -rf "$OUT"; mkdir -p "$OUT"
echo "building figures with $PY"

# collect <src.png> <name.png> <caption>
CAPS="$OUT/.caps"; : > "$CAPS"
grab() {  # grab <src> <newname> <caption>
  if [ -f "$1" ]; then
    cp "$1" "$OUT/$2"
    printf '%s\t%s\n' "$2" "$3" >> "$CAPS"
    ok "$3"
  else
    skip "$3  (not built yet)"
  fi
}

# ------------------------------------------------------- 1. raw data, both exps
grab "$CHECK/ca_data.png" "01_data.png" \
     "Ca data: Indo-1 ratio mean +/- SE, Experiment 1 and Experiment 2"

# ------------------------------------------------------- 2. Exp 1 fits
if [ -d "$FITS" ]; then
  ( cd "$FITS" && "$PY" plot_ca_fit.py orig   >/dev/null 2>&1 )
  ( cd "$FITS" && "$PY" plot_pzap_fit.py      >/dev/null 2>&1 )
  for p in "$FITS"/*.png; do
    [ -f "$p" ] || continue
    b="$(basename "$p")"
    case "$b" in
      *ca*)   grab "$p" "02_$b" "Experiment 1 Ca fit, variant orig (her 5 parameters)" ;;
      *pzap*) grab "$p" "03_$b" "pZAP70 fit" ;;
      *)      grab "$p" "04_$b" "$b" ;;
    esac
  done
else
  skip "nk92-fits repo not at $FITS"
fi

# ------------------------------------------------------- 3. Exp 1 bands
if [ -d "$HOME/Ca_fit_c02/out_boot_ca2" ]; then
  for v in orig; do
    n=$(ls "$HOME/Ca_fit_c02/out_boot_ca2/$v"/boot_*.json 2>/dev/null | wc -l)
    if [ "$n" -gt 2 ]; then
      ( cd "$HOME/Ca_fit_c02" && "$PY" "$BOOT" band "$v" >/dev/null 2>&1 )
      grab "$HOME/Ca_fit_c02/out_boot_ca2/band_$v.png" "05_exp1_band_$v.png" \
           "Experiment 1, $v: bootstrap prediction envelope ($n samples)"
    else
      skip "Experiment 1 band_$v  ($n samples so far)"
    fi
  done
fi

# ------------------------------------------------------- 4. Exp 2 bands
if [ -d "$EXP2/out_boot_ca2" ]; then
  for v in orig orig1 lin1 min4; do
    n=$(ls "$EXP2/out_boot_ca2/$v"/boot_*.json 2>/dev/null | wc -l)
    if [ "$n" -gt 2 ]; then
      ( cd "$EXP2" && "$PY" "$BOOT" band "$v" >/dev/null 2>&1 )
      grab "$EXP2/out_boot_ca2/band_$v.png" "06_exp2_band_$v.png" \
           "Experiment 2, $v: bootstrap prediction envelope ($n samples)"
    else
      skip "Experiment 2 band_$v  ($n samples so far)"
    fi
  done
fi

# ------------------------------------------------------------------ index page
"$PY" - "$OUT" <<'PYIDX'
import os, sys, html, datetime
out = sys.argv[1]
caps = {}
order = []
cf = os.path.join(out, '.caps')
if os.path.exists(cf):
    for line in open(cf):
        if '\t' in line:
            n, c = line.rstrip('\n').split('\t', 1)
            caps[n] = c; order.append(n)
order = sorted(set(order))
rows = []
for n in order:
    rows.append(
        f'<figure><img src="{html.escape(n)}" alt="{html.escape(caps[n])}">'
        f'<figcaption>{html.escape(caps[n])}</figcaption></figure>')
body = '\n'.join(rows) or '<p class="empty">No figures built yet. Jobs still running?</p>'
open(os.path.join(out, 'index.html'), 'w').write(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NK92 CD16 fits</title>
<style>
 body {{ background:#fff; color:#1a1a1a; font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;
        margin:0; padding:24px 16px 60px; }}
 .wrap {{ max-width:1100px; margin:0 auto; }}
 h1 {{ font-size:20px; margin:0 0 4px; }}
 .sub {{ color:#666; font-size:13px; margin:0 0 28px; }}
 figure {{ margin:0 0 34px; border:1px solid #e3e3e3; border-radius:8px; overflow:hidden; }}
 img {{ display:block; width:100%; height:auto; background:#fafafa; }}
 figcaption {{ padding:10px 14px; font-size:13px; color:#444; border-top:1px solid #e3e3e3; }}
 .empty {{ color:#666; }}
</style></head><body><div class="wrap">
<h1>NK92 CD16 &mdash; current fits</h1>
<p class="sub">built {datetime.datetime.now():%Y-%m-%d %H:%M} &middot; {len(order)} figure(s)</p>
{body}
</div></body></html>""")
print(f'index: {len(order)} figure(s)')
PYIDX
rm -f "$CAPS"

# ---------------------------------------------------------------------- serve
pkill -f "http.server $PORT" 2>/dev/null
( cd "$HOME" && setsid nohup "$PY" -m http.server "$PORT" >/dev/null 2>&1 </dev/null & ) 2>/dev/null \
  || ( cd "$HOME" && nohup "$PY" -m http.server "$PORT" >/dev/null 2>&1 </dev/null & )
sleep 1
IP="$(hostname -I 2>/dev/null | awk '{print $1}')"
echo
echo "   open this:"
echo "   http://${IP:-10.73.170.129}:$PORT/$(basename "$OUT")/"
echo
