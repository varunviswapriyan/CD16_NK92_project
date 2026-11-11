#!/usr/bin/env bash
# check_fboff.sh -- early sanity check of the feedback-off runs (run anytime).
# 1) the deck edit is really in the simulation input (max_ifng_production 2.0 -> 0.5)
# 2) free IFN-g in feedback-off vs the standard run at the SAME time, seed s9001:
#    feedback off should be LOWER (no boost), most clearly in females.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash check_fboff.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
for sex in male female; do
  f=runs_big/fboff/$sex/s9001
  echo "=== $sex ==="
  if grep -q "^app_style.* 4 4 0.5 1000 1000 1000 " $f/in.simulation 2>/dev/null; then
    echo "  deck edit present: feedback max 0.5  OK"
  else
    echo "  !! deck edit NOT found in $f/in.simulation -- stop and tell Claude"; continue
  fi
  python3 - "$f/log.spparks" "runs_big/prod/$sex/s9001/log.spparks" <<'PY'
import sys
def read(p):
    hdr=None; rows=[]
    for line in open(p):
        s=line.split()
        if "sum(i1)" in s:
            if rows: break
            hdr=s; continue
        if hdr and len(s)==len(hdr):
            try: rows.append([float(x) for x in s])
            except ValueError: pass
    return hdr,rows
h,fb=read(sys.argv[1]); h2,st=read(sys.argv[2])
if not fb: print("  no data yet"); sys.exit()
ci=h.index("sum(i1)"); di=h.index("sum(d1)"); pi=h.index("Propnsty") if "Propnsty" in h else None
t=fb[-1][0]
ref=min(st,key=lambda r:abs(r[0]-t))
print(f"  time {t:.0f} (of 1008)   propensity {fb[-1][pi] if pi is not None else '?'}")
print(f"  free IFN-g : feedback off {fb[-1][di]:.3g}   standard {ref[di]:.3g}   ratio {fb[-1][di]/max(ref[di],1e-9):.2f}")
print(f"  cancer     : feedback off {fb[-1][ci]:.0f}   standard {ref[ci]:.0f}")
PY
done
echo
echo "Looks right if: deck edit OK, no 'nan', time increasing, free IFN-g ratio < 1 (clearest in female)."
