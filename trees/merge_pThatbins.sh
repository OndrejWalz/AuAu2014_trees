#!/usr/bin/env bash
# Merge EMBEDDING outputs per pThat bin (one file per job_ptX_Y).
# Reads:  <BASE>/job_ptX_Y/production/*.root
# Writes: <OUT_DIR>/embedding_ptX_Y.root
#
# Usage:
#   ./merge_embedding_per_bin.sh /path/to/trees [OUT_DIR]
#
# Defaults:
#   OUT_DIR = <BASE>/merged_all
#
# Env:
#   MIN_SIZE=5000    # ignore tiny files (bytes)
#   BATCH=200        # files per shard
#   CLEAN_SHARDS=1   # remove shard_*.root after final merge

set -euo pipefail

BASE="${1:-.}"
OUT_DIR="${2:-$BASE/merged_all}"

MIN_SIZE="${MIN_SIZE:-5000}"
BATCH="${BATCH:-200}"
CLEAN_SHARDS="${CLEAN_SHARDS:-1}"

# Fixed list of EMBEDDING bins (no data)
BINS=( job_pt3_5 job_pt5_7 job_pt7_9 job_pt9_11 job_pt11_15 job_pt15_20
       job_pt20_25 job_pt25_30 job_pt30_40 job_pt40_50 job_pt50_-1 )

mkdir -p "$OUT_DIR"

# Detect hadd flags (older ROOT may not support -k)
HADD_FLAGS="-f"
if hadd -h 2>&1 | grep -q -- '-k'; then
  HADD_FLAGS="-fk"
fi

merge_one_bin () {
  local bin="$1"
  local prod="$BASE/$bin/production"
  local tag="${bin#job_}"                 # e.g. pt3_5
  local final="$OUT_DIR/embedding_${tag}.root"

  if [[ ! -d "$prod" ]]; then
    echo "  [skip] $bin: no production/ directory"
    return 0
  fi

  local tmpdir="$OUT_DIR/.${bin}_tmp"
  mkdir -p "$tmpdir"
  local list="$tmpdir/inlist.txt"
  : > "$list"

  # Collect ROOT files above size threshold (sorted)
  local has=0
  while IFS= read -r -d '' f; do
    local sz
    sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
    if [[ "$sz" -gt "$MIN_SIZE" ]]; then
      echo "$f" >> "$list"
      has=1
    fi
  done < <(find "$prod" -type f -name '*.root' -print0 | sort -z)

  if [[ "$has" -eq 0 ]]; then
    echo "  [note] $bin: no usable ROOT files"
    rm -rf "$tmpdir"
    return 0
  fi

  local N
  N=$(wc -l < "$list" | tr -d ' ')
  if [[ "$N" -eq 1 ]]; then
    cp -f "$(head -n1 "$list")" "$final"
    echo "  [ok]   $bin -> $(basename "$final") (1 file)"
    rm -rf "$tmpdir"
    return 0
  fi

  # Shard to avoid argv-too-long; old-ROOT friendly
  rm -f "$tmpdir"/shard_*.root "$tmpdir"/inlist.batch.* || true
  split -l "$BATCH" -d -a 4 "$list" "$tmpdir/inlist.batch."

  local idx=0
  for sub in "$tmpdir"/inlist.batch.*; do
    local shard
    shard=$(printf "%s/shard_%04d.root" "$tmpdir" "$idx")
    hadd $HADD_FLAGS "$shard" $(tr '\n' ' ' < "$sub") >/dev/null
    idx=$((idx+1))
  done

  # Final per-bin merge (trees + histos)
  hadd $HADD_FLAGS "$final" "$tmpdir"/shard_*.root >/dev/null

  # Cleanup
  if [[ "$CLEAN_SHARDS" -eq 1 ]]; then
    rm -rf "$tmpdir"
  else
    rm -f "$tmpdir"/inlist.txt "$tmpdir"/inlist.batch.* || true
  fi

  echo "  [ok]   $bin -> $(basename "$final") ($N files)"
}

echo "Merging per pThat bin into: $OUT_DIR"
for bin in "${BINS[@]}"; do
  merge_one_bin "$bin"
done
echo "Done."
