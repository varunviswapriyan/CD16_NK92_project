#!/usr/bin/env bash
# find_matching_deck.sh -- compare our prod deck against every input deck in
# ~/giuseppe_ref and report how many lines differ, so we find which of
# Giuseppe's files IS our manuscript input (0 diff = identical).
# Run:  cd ~/CD16_NK92_project && git pull && bash find_matching_deck.sh
set -uo pipefail
OURS=~/CD16_NK92_project/SPPARKS_Sex_Bias/rsc_out_2026-07-03/in.full
REF=~/giuseppe_ref
[ -f "$OURS" ] || { echo "our deck not found: $OURS"; exit 1; }

echo "comparing OUR prod deck against each of Giuseppe's input decks"
echo "(difflines = number of differing lines; 0 = identical)"
echo
# ignore comment-only differences too, so we judge the actual directives
for f in "$REF"/in.*; do
  [ -f "$f" ] || continue
  raw=$(diff "$OURS" "$f" 2>/dev/null | grep -c '^[<>]')
  # strip blank + comment lines before comparing, to catch "same rules, different notes"
  code=$(diff <(grep -vE '^\s*(#|$)' "$OURS") <(grep -vE '^\s*(#|$)' "$f") 2>/dev/null | grep -c '^[<>]')
  printf "%-55s  difflines=%-5s  code-only=%-5s\n" "$(basename "$f")" "$raw" "$code"
done | sort -t= -k3 -n

echo
echo ">> the file with code-only=0 is your manuscript deck (same directives)."
echo ">> if none are 0, the smallest code-only is the closest ancestor."
