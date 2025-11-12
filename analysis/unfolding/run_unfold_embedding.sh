#!/usr/bin/env bash
set -euo pipefail

# --- edit these three if your layout changes ---
BASE="/gpfs/mnt/gpfs01/star/pwg/svomich/JetsTrees"
SIF="/gpfs/mnt/gpfs01/star/pwg/prozorov/singularity/roounfold.sif"
MACRO="/gpfs/mnt/gpfs01/star/pwg/svomich/JetsTrees/analysis/unfolding/unfold.cxx"
# ----------------------------------------------

IN="${1:-$BASE/trees/merged_all/embedding_merged.root}"
OUT="${2:-$BASE/analysis/unfolding/out}"

if [[ ! -f "$IN" ]]; then
  echo "Input not found: $IN"
  exit 1
fi

mkdir -p "$OUT"

apptainer exec -e -B /gpfs/mnt/gpfs01 \
  "$SIF" \
  root -l -b -q 'gSystem->Load("libRooUnfold");' "$MACRO++(\"$IN\",\"$OUT\")"


