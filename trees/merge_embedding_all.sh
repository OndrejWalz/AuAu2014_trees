#!/usr/bin/env bash
# Merge embedding outputs (trees + histograms) from a fixed set of pThat bins.
# Only merges from these folders (if present):
#   job_pt3_5 job_pt5_7 job_pt7_9 job_pt9_11 job_pt11_15 job_pt15_20
#   job_pt20_25 job_pt25_30 job_pt30_40 job_pt40_50 job_pt50_-1
#
# Usage:
#   ./merge_embedding_from_bins.sh /jets/trees [OUT_DIR] [OUT_NAME]
# Defaults:
#   OUT_DIR  = <BASE>/merged_all
#   OUT_NAME = embedding_merged.root
#
# Env:
#   MIN_SIZE=5000   # ignore tiny files (bytes)
#   BATCH=200       # files per shard
#   CLEAN_SHARDS=1  # remove shard_*.root after final merge

set -euo pipefail

BASE="${1:-.}"
OUT_DIR="${2:-$BASE/merged_all}"
OUT_NAME="${3:-embedding_merged.root}"

MIN_SIZE="${MIN_SIZE:-5000}"
BATCH="${BATCH:-200}"
CLEAN_SHARDS="${CLEAN_SHARDS:-1}"

# Fixed list of embedding bins
BINS=( job_pt3_5 job_pt5_7 job_pt7_9 job_pt9_11 job_pt11_15 job_pt15_20
       job_pt20_25 job_pt25_30 job_pt30_40 job_pt40_50 job_pt50_-1 )

mkdir -p "$OUT_DIR"

# Detect hadd flags (old ROOT may lack -k)
HADD_FLAGS="-f"
if hadd -h 2>&1 | grep -q -- '-k'; then
  HADD_FLAGS="-fk"
fi

FINAL="$OUT_DIR/$OUT_NAME"
LIST="$OUT_DIR/inlist.txt"
: > "$LIST"

echo "Scanning embedding production files in fixed pThat bins..."
found_any=0
for bin in "${BINS[@]}"; do
  prod="$BASE/$bin/production"
  if [[ ! -d "$prod" ]]; then
    echo "  [skip] $bin: no production/ directory"
    continue
  fi

  # Collect ROOT files above size threshold
  has_files=0
  while IFS= read -r -d '' f; do
    sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
    if [[ "$sz" -gt "$MIN_SIZE" ]]; then
      echo "$f" >> "$LIST"
      has_files=1
      found_any=1
    fi
  done < <(find "$prod" -type f -name '*.root' -print0 | sort -z)

  if [[ "$has_files" -eq 0 ]]; then
    echo "  [note] $bin: no usable ROOT files"
  else
    echo "  [ok]   $bin"
  fi
done

if [[ "$found_any" -eq 0 ]]; then
  echo "No embedding files found in the specified bins under: $BASE"
  exit 1
fi

N=$(wc -l < "$LIST" | tr -d ' ')
# Single file fast path
if [[ "$N" -eq 1 ]]; then
  src=$(head -n1 "$LIST")
  cp -f "$src" "$FINAL"
  rm -f "$LIST"
  echo "Merged 1 file -> $FINAL"
  exit 0
fi

# Shard to avoid argv-too-long; old-ROOT friendly
rm -f "$OUT_DIR"/shard_*.root "$OUT_DIR"/inlist.batch.* || true
split -l "$BATCH" -d -a 4 "$LIST" "$OUT_DIR/inlist.batch."

shard_idx=0
for sub in "$OUT_DIR"/inlist.batch.*; do
  shard=$(printf "%s/shard_%04d.root" "$OUT_DIR" "$shard_idx")
  # expand modest number of args per shard
  hadd $HADD_FLAGS "$shard" $(tr '\n' ' ' < "$sub") >/dev/null
  shard_idx=$((shard_idx + 1))
done

# Final merge: trees + histos
hadd $HADD_FLAGS "$FINAL" "$OUT_DIR"/shard_*.root >/dev/null

# Cleanup
rm -f "$OUT_DIR"/inlist.txt "$OUT_DIR"/inlist.batch.* || true
if [[ "$CLEAN_SHARDS" -eq 1 ]]; then
  rm -f "$OUT_DIR"/shard_*.root || true
fi

echo "Merged $N files -> $FINAL"
