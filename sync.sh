#!/bin/bash
# =============================================================================
# sync.sh -- move freshly downloaded scripts into the repo, commit, push.
#            Lives in ~/Downloads/indrani_code/ and replaces typing the
#            mv / git add / git commit / git push chain every time.
#
#   bash ~/Downloads/indrani_code/sync.sh              # everything downloaded today
#   bash ~/Downloads/indrani_code/sync.sh pairs.sh     # only these
#   bash ~/Downloads/indrani_code/sync.sh -n           # dry run, show and stop
#   HOURS=72 bash ~/Downloads/indrani_code/sync.sh     # widen the time window
#
# Handles the browser's duplicate naming: "pairs (1).sh" lands as "pairs.sh",
# and if several copies of the same file are in Downloads the NEWEST wins.
# =============================================================================
set -o pipefail

REPO=${REPO:-$HOME/Downloads/indrani_code}
DEST=${DEST:-$REPO/filesCC}
DL=${DL:-$HOME/Downloads}
HOURS=${HOURS:-24}

DRY=""
if [ "$1" = "-n" ] || [ "$1" = "--dry-run" ]; then DRY=1; shift; fi

[ -d "$REPO" ] || { echo "ERROR: repo not found at $REPO"; exit 1; }
[ -d "$DEST" ] || { echo "ERROR: $DEST not found"; exit 1; }

# ---- collect candidates ----------------------------------------------------
tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
if [ $# -gt 0 ]; then
  for f in "$@"; do
    if   [ -f "$f" ];        then echo "$f"     >> "$tmp"
    elif [ -f "$DL/$f" ];    then echo "$DL/$f" >> "$tmp"
    else echo "  skip: $f not found (looked in . and $DL)" >&2
    fi
  done
else
  find "$DL" -maxdepth 1 \( -name '*.sh' -o -name '*.py' \) \
       -mmin -$((HOURS * 60)) 2>/dev/null | sort >> "$tmp"
fi

if [ ! -s "$tmp" ]; then
  echo "Nothing to sync (no .sh/.py in $DL from the last $HOURS h)."
  echo "Widen it with:  HOURS=72 bash \$0"
  exit 0
fi

# ---- normalise "name (1).sh" -> "name.sh", newest copy wins -----------------
# Build "target<TAB>source", keep the newest source per target.
pick=$(mktemp); trap 'rm -f "$tmp" "$pick" "${pick}.go" "${pick}.old"' EXIT
while IFS= read -r f; do
  [ -f "$f" ] || continue
  base=$(basename "$f")
  tgt=$(printf '%s' "$base" | sed -E 's/[ _]?\(([0-9]+)\)(\.[A-Za-z0-9]+)$/\2/')
  printf '%s\t%s\n' "$tgt" "$f" >> "$pick"
done < "$tmp"

targets=$(cut -f1 "$pick" | sort -u)   # basenames only; no spaces survive normalisation

echo "About to sync into $DEST:"
n=0
for t in $targets; do
  # newest source for this target -- bash's -nt, so filenames with spaces
  # like "pairs (1).sh" are safe and no GNU-only tool is needed
  src=""
  while IFS=$'\t' read -r a b; do
    [ "$a" = "$t" ] || continue
    if [ -z "$src" ] || [ "$b" -nt "$src" ]; then src="$b"; fi
  done < "$pick"
  [ -f "$src" ] || continue
  # every OTHER copy of this target is superseded -- park it, or the next run
  # would pick the stale one up and quietly overwrite the newer file
  while IFS=$'\t' read -r a b; do
    [ "$a" = "$t" ] || continue
    [ "$b" = "$src" ] && continue
    echo "$b" >> "${pick}.old"
  done < "$pick"

  if [ -f "$DEST/$t" ] && cmp -s "$src" "$DEST/$t"; then
    printf '  = %-22s unchanged, skipping\n' "$t"; continue
  fi
  printf '  + %-22s <- %s\n' "$t" "$(basename "$src")"
  echo "$t"$'\t'"$src" >> "${pick}.go"
  n=$((n + 1))
done

sweep_old () {                      # park superseded duplicates out of the way
  [ -f "${pick}.old" ] || return 0
  local park="$DL/.synced_old" b k=0
  mkdir -p "$park" 2>/dev/null || return 0
  while IFS= read -r b; do
    [ -f "$b" ] || continue
    mv -f "$b" "$park/" 2>/dev/null && k=$((k + 1))
  done < "${pick}.old"
  rm -f "${pick}.old"
  [ "$k" -gt 0 ] && echo "  moved $k superseded duplicate(s) to $park"
  return 0
}

if [ "$n" -eq 0 ]; then
  echo
  [ -z "$DRY" ] && sweep_old
  echo "Everything already matches what is in the repo.  Nothing to do."
  exit 0
fi

if [ -n "$DRY" ]; then
  echo
  echo "(dry run -- nothing moved.  Rerun without -n to do it.)"
  rm -f "${pick}.go" "${pick}.old"; exit 0
fi

# ---- move ------------------------------------------------------------------
echo
while IFS=$'\t' read -r t src; do
  mv -f "$src" "$DEST/$t" || { echo "ERROR moving $src"; exit 1; }
done < "${pick}.go"
names=$(cut -f1 "${pick}.go" | tr '\n' ' ')
rm -f "${pick}.go"
sweep_old

# ---- commit and push -------------------------------------------------------
cd "$REPO" || exit 1
git add filesCC/ || exit 1
if git diff --cached --quiet; then
  echo "git reports no staged change -- nothing to commit."
  exit 0
fi
MSG=${MSG:-"update: $names"}
git commit -m "$MSG" || exit 1
if ! git push; then
  echo
  echo "Commit succeeded but push FAILED.  Usually someone else pushed first:"
  echo "    cd $REPO && git pull --rebase && git push"
  exit 1
fi

echo
echo "=============================================================="
echo "pushed: $names"
echo
echo "now on the cluster:"
echo "    cd ~/CD16_NK92_project && git pull"
echo "=============================================================="
