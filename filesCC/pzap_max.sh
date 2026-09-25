#!/bin/bash
# pzap_max.sh -- the two runs that between them attack every fixable CI, based
#                on last night's diagnostic.  Nothing already computed is lost.
#
#   bash ~/CD16_NK92_project/filesCC/pzap_max.sh          # submit both (clears queue first)
#   bash ~/CD16_NK92_project/filesCC/pzap_max.sh add      # submit only what is NOT
#                                                        # already running; cancels nothing
#   bash ~/CD16_NK92_project/filesCC/pzap_max.sh report   # both tables
#
# WHAT LAST NIGHT SHOWED, AND WHAT EACH ARM DOES ABOUT IT
#
#   lig0   w 1.39   only 10% simulator noise  -> DATA-limited  -> arm B
#   kd10   w 1.29   only 10% simulator noise  -> DATA-limited  -> arm B
#   KPR    w 1.26   128% simulator noise      -> NOISE-limited -> arm A2
#   SYK0   w 0.77   211% simulator noise      -> NOISE-limited -> arm A2 (partly)
#   ZAP0   w 0.65   replicate sat on a bound  -> BOX-limited   -> arm A2
#   kzp    w 0.50   fine already
#
#   arm A2  0-300 s, Indrani's approved window.  N_REPS 6 -> 16 (noise SD x 0.61)
#           and the search box widened 0.25 -> 0.35 log10 so ZAP0 stops hitting
#           a wall.  ~8 h.   -> ~/boot_pzap_n16
#   arm B   0-600 s, adds the 10-minute point: 12 informative points instead of
#           9.  The only lever that exists for lig0 and kd10.  ~4 h.
#           -> ~/boot_pzap_t600
#
# Last night's results stay in ~/boot_pzap, untouched, as the fallback.
set -e
F=~/CD16_NK92_project/filesCC

# ---- which roots have jobs on the nodes right now --------------------------
running_roots() {
  for id in $(squeue -u "$USER" -h -o "%i" 2>/dev/null); do
    scontrol show job "$id" 2>/dev/null | tr ' ' '\n' | grep '^StdOut=' | sed 's|^StdOut=||'
  done
}

submit_A2() {
  echo "################  ARM A2 -- 300 s, N_REPS=16, box 0.35  ################"
  BOOT_NOCANCEL=1 BOOT_TAG=_n16 BOOT_NREPS=16 BOOT_ITERS=16 BOOT_HALF=0.35 \
    bash $F/pzap_ci6.sh
}
submit_B() {
  echo "################  ARM B -- 600 s, 10-minute point  ################"
  BOOT_NOCANCEL=1 bash $F/pzap_ci10.sh
}

if [ "$1" = "add" ]; then
  R=$(running_roots)
  echo "jobs currently on the nodes, by output directory:"
  if [ -z "$R" ]; then echo "   (none)"; else
    echo "$R" | sed 's|.*/\(boot_pzap[^/]*\)/.*|   \1|' | sort | uniq -c
  fi
  echo
  a2=0; b=0
  echo "$R" | grep -q '/boot_pzap_n16/'  && a2=1
  echo "$R" | grep -q '/boot_pzap_t600/' && b=1
  [ $a2 = 1 ] && echo "arm A2 is already running -- leaving it."  || submit_A2
  echo
  [ $b  = 1 ] && echo "arm B  is already running -- leaving it."  || submit_B
  echo
  echo "Nothing already running was cancelled.   squeue -u \$USER"
  echo "When they land:  bash \$0 report"
  exit 0
fi

if [ "$1" = "report" ]; then
  for spec in "ARM A2  0-300 s, N_REPS=16, wider box|_n16|pzap_ci6.sh" \
              "ARM B   0-600 s, 10-min point|  |pzap_ci10.sh" \
              "ARM A   0-300 s, last night (reference)||pzap_ci6.sh"; do
    IFS='|' read -r title tag script <<< "$spec"
    echo "################  ${title}  ################"
    BOOT_TAG="$(echo $tag)" BOOT_NOSERVE=1 bash $F/$script report 2>/dev/null \
      || echo "(nothing to report yet)"
    echo
  done
  OUT=~/ci_plots; mkdir -p $OUT
  IP=$(hostname -I | awk '{for(i=1;i<=NF;i++) if($i ~ /^10\.73\./) print $i}')
  echo "Open:  http://${IP:-10.73.170.128}:8000/"
  cd $OUT && $(command -v python3 || command -v python) -m http.server 8000
  exit 0
fi

# clear the queue ONCE here; the arms are told not to (they would kill each other)
if [ "$1" != "add" ]; then
  OLD=$(squeue -u "$USER" -h -o "%i %j" | awk '$2 ~ /^bp/ {print $1}')
  [ -n "$OLD" ] && { echo "$OLD" | xargs -r scancel; echo "cleared queue"; sleep 2; }
else
  echo "add mode: existing pZAP jobs are left running"
fi

submit_A2
echo
submit_B

echo
echo "Both submitted.  arm A2 ~8 h, arm B ~4 h.   squeue -u \$USER"
echo "When they land:  bash \$0 report"
