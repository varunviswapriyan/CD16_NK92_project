#!/usr/bin/env bash
# check_running.sh -- sanity check on the IFNGR-blockade and corrected-CD4 runs
# while they are still going. For seed s9001 of each: latest Time (should climb to
# 1008), Propnsty (must be a number, not nan), cancer, AR- CD8, CD4 helpers, free IFN-g.
# Also: how many seeds are done, and any nan anywhere.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash check_running.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
echo "=== queue ==="; squeue -u "$USER" -h -o "%j %T" | sort | uniq -c
for d in runs_big/blockall/male runs_big/blockall/female runs_big/cd4off_v2/male runs_big/cd4off_v2/female; do
  echo; echo "=== $d ==="
  logs=$(ls $d/s*/log.spparks 2>/dev/null)
  [ -z "$logs" ] && { echo "  no logs yet"; continue; }
  done_n=0; nan_n=0
  for f in $logs; do
    grep -qi nan "$f" && nan_n=$((nan_n+1))
    t=$(awk '/sum\(i1\)/{h=1;next} h&&NF>20&&$1+0==$1{t=$1} END{print t+0}' "$f")
    awk -v t="$t" 'BEGIN{exit !(t>=1000)}' && done_n=$((done_n+1))
  done
  echo "  seeds with logs: $(echo $logs | wc -w)   finished (Time>=1000): $done_n   with nan: $nan_n"
  f=$(echo $logs | tr ' ' '\n' | head -1)
  echo "  $f  (first row, latest row):"
  awk '/sum\(i1\)/{for(i=1;i<=NF;i++)c[$i]=i; h=1; next}
       h&&NF>20&&$1+0==$1{ if(!n++) first=$0; last=$0 }
       END{ split(first,a); split(last,b);
            printf "     %-8s %-10s %-9s %-8s %-8s %-8s %-12s\n","Time","Propnsty","cancer","AR-CD8","CD4AR-","CD4AR+","freeIFNg";
            printf "     %-8s %-10s %-9s %-8d %-8s %-8s %-12s\n",a[c["Time"]],a[c["Propnsty"]],a[c["sum(i1)"]],a[c["sum(i2)"]]+a[c["sum(i3)"]]+a[c["sum(i4)"]],a[c["sum(i13)"]],a[c["sum(i14)"]],a[c["sum(d1)"]];
            printf "     %-8s %-10s %-9s %-8d %-8s %-8s %-12s\n",b[c["Time"]],b[c["Propnsty"]],b[c["sum(i1)"]],b[c["sum(i2)"]]+b[c["sum(i3)"]]+b[c["sum(i4)"]],b[c["sum(i13)"]],b[c["sum(i14)"]],b[c["sum(d1)"]] }' "$f"
done
cat <<'TXT'

What "working" looks like:
  * Propnsty is a number (never nan); nan count = 0
  * Time on the latest row increases each time you rerun this
  * cancer changes from the first row (not stuck at 136840 / 186840)
  * cd4off_v2: CD4AR- and CD4AR+ = 0, AR-CD8 still present and recovering (recruitment intact)
  * blockall: free IFN-g rising (T cells no longer bind it)
TXT
