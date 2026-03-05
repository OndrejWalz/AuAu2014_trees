#!/usr/bin/env bash
set -euo pipefail

########################
# Paths relative to this script
########################

# Directory where *this* script lives (analysis/unfolding)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# JetsTrees base (two levels up)
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Files relative to this structure
SIF="${SCRIPT_DIR}/roounfold.sif"
MACRO="${SCRIPT_DIR}/compare_bayes_svd.C"

########################
# Arguments
########################

# 1st arg: Bayes input file (basename in out_embedding_Bayes/ or absolute path)
if [[ $# -ge 1 ]]; then
  if [[ "$1" = /* ]]; then
    BAYES_INPUT="$1"
  else
    BAYES_INPUT="${SCRIPT_DIR}/$1"
  fi
else
  BAYES_INPUT="${SCRIPT_DIR}/out_embedding_BAYES/responses_embedding.root"
fi

# 2nd arg: SVD input file (basename in out_embedding_SVD/ or absolute path)
if [[ $# -ge 2 ]]; then
  if [[ "$2" = /* ]]; then
    SVD_INPUT="$2"
  else
    SVD_INPUT="${SCRIPT_DIR}/$2"
  fi
else
  SVD_INPUT="${SCRIPT_DIR}/out_embedding_SVD/responses_embedding.root"
fi

# 3rd arg: output directory
OUT_DIR="${3:-${SCRIPT_DIR}/out_comparison}"

########################
# Checks
########################

echo "----------------------------------------"
echo "Running Bayes vs SVD comparison"
echo "SCRIPT_DIR  : $SCRIPT_DIR"
echo "BASE        : $BASE"
echo "SIF         : $SIF"
echo "Macro       : $MACRO"
echo "Bayes input : $BAYES_INPUT"
echo "SVD input   : $SVD_INPUT"
echo "Output dir  : $OUT_DIR"
echo "----------------------------------------"

[[ -f "$SIF"         ]] || { echo "ERROR: SIF not found:         $SIF";         exit 1; }
[[ -f "$MACRO"       ]] || { echo "ERROR: Macro not found:       $MACRO";       exit 1; }
[[ -f "$BAYES_INPUT" ]] || { echo "ERROR: Bayes input not found: $BAYES_INPUT"; exit 1; }
[[ -f "$SVD_INPUT"   ]] || { echo "ERROR: SVD input not found:   $SVD_INPUT";   exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

apptainer exec -e -B /gpfs01 \
  "$SIF" \
  root -l -b <<EOF
gSystem->Load("libRooUnfold");
.x ${MACRO}+("${BAYES_INPUT}","${SVD_INPUT}","${OUT_DIR}");
.q
EOF

echo "----------------------------------------"
echo "Done. Plots written to: $OUT_DIR"
echo "----------------------------------------"