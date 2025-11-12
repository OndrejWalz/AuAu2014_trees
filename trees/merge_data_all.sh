#!/usr/bin/env bash
# Merge DATA outputs (trees + histograms) from job_pico_low_14 into one file.
#
# Usage:
#   ./merge_data_job.sh /jets/trees [OUT_DIR] [OUT_NAME]
#
# Defaults:
#   OUT_DIR  = <BASE>/merged_all
#   OUT_NAME = data_merged.root
#
# Env:
#   MIN_SIZE=5000   # ignore tiny files (bytes)
#   BATCH=200       # files per shard
#   CLEAN_SHARDS=1  # remove shard_*.root after final merge

set -euo pipefail

BASE="${1:-.}"
OUT_DIR="${2:-$BASE/merged_all}"
OUT_NAME="${3:-data_merged.root}"

MIN_SIZE="${MIN_SIZE:-5000}"
BATCH="${BATCH:-200}"
CLEAN_SHARDS="${CLEAN_SHARDS:-1}"

PROD_DIR="$BASE/job_pico_low_14/production"

mkdir -p "$OUT_DIR"

# Detect hadd flags (old ROOT may not support -k)
HADD_FLAGS="-f"
if hadd -h 2>&1 | grep -q -- '-k'; then
  HADD_FLAGS="-fk"
fi

FINAL="$OUT_DIR/$OUT_NAME"
LIST="$OUT_DIR/inlist_data.txt"
: > "$LIST"

echo "Scanning data production files in: $PROD_DIR"
if [[ ! -d "$PROD_DIR" ]]; then
  echo "No data production directory found (expected: $PROD_DIR)"
  exit 1
fi

# Collect ROOT files above size threshold (sorted)
has_files=0
while IFS= read -r -d '' f; do
  sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
  if [[ "$sz" -gt "$MIN_SIZE" ]]; then
    echo "$f" >> "$LIST"
    has_files=1
  fi
done < <(find "$PROD_DIR" -type f -name '*.root' -print0 | sort -z)

if [[ "$has_files" -eq 0 ]]; then
  echo "No usable ROOT files in data production directory."
  exit 1
fi

N=$(wc -l < "$LIST" | tr -d ' ')

# Fast path for a single file
if [[ "$N" -eq 1 ]]; then
  src=$(head -n1 "$LIST")
  cp -f "$src" "$FINAL"
  rm -f "$LIST"
  echo "Merged 1 file -> $FINAL"
  exit 0
fi

# Shard to avoid argv-too-long; old-ROOT friendly
rm -f "$OUT_DIR"/data_shard_*.root "$OUT_DIR"/inlist_data.batch.* || true
split -l "$BATCH" -d -a 4 "$LIST" "$OUT_DIR/inlist_data.batch."

shard_idx=0
for sub in "$OUT_DIR"/inlist_data.batch.*; do
  shard=$(printf "%s/data_shard_%04d.root" "$OUT_DIR" "$shard_idx")
  # Expand a modest number of args per shard; quiet output
  hadd $HADD_FLAGS "$shard" $(tr '\n' ' ' < "$sub") >/dev/null
  shard_idx=$((shard_idx + 1))
done

# Final merge (trees + histos)
hadd $HADD_FLAGS "$FINAL" "$OUT_DIR"/data_shard_*.root >/dev/null

# Cleanup
rm -f "$OUT_DIR"/inlist_data.txt "$OUT_DIR"/inlist_data.batch.* || true
if [[ "$CLEAN_SHARDS" -eq 1 ]]; then
  rm -f "$OUT_DIR"/data_shard_*.root || true
fi

echo "Merged $N files -> $FINAL"
