#!/usr/bin/env bash
set -euo pipefail

########################
# Configuration
########################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# reuse the same container you already have in unfolding/
SIF="${BASE}/analysis/unfolding/roounfold.sif"

# macro in this folder
MACRO="${SCRIPT_DIR}/rcp_check.cxx"

########################
# Arguments
########################
# 1st arg: unfolded_data.root (absolute or basename under unfolding/out_data)
if [[ $# -ge 1 ]]; then
  if [[ "$1" = /* ]]; then
    INPUT="$1"
  else
    # default assumption: user passes e.g. out_data/unfolded_data.root or similar
    INPUT="${BASE}/analysis/unfolding/$1"
  fi
else
  INPUT="${BASE}/analysis/unfolding/out_data_SVD_ME/unfolded_data.root"
fi

# 2nd arg: output dir
OUT_DIR="${2:-${SCRIPT_DIR}/out_rcp_SVD_ME}"

# 3rd arg: histogram name in unfolded file
#HNAME="${3:-hUnfoldedTruthBins_rawSVD_reg4}" # use for SVD
HNAME="${3:-hUnfoldedTruthBins_matchCorr}" # use for Bayes

########################
# Checks
########################
echo "----------------------------------------"
echo "Running R_CP quick check"
echo "SCRIPT_DIR : $SCRIPT_DIR"
echo "BASE       : $BASE"
echo "SIF        : $SIF"
echo "MACRO      : $MACRO"
echo "INPUT      : $INPUT"
echo "OUT_DIR    : $OUT_DIR"
echo "HIST       : $HNAME"
echo "----------------------------------------"

[[ -f "$SIF"   ]] || { echo "ERROR: SIF not found: $SIF"; exit 1; }
[[ -f "$MACRO" ]] || { echo "ERROR: MACRO not found: $MACRO"; exit 1; }
[[ -f "$INPUT" ]] || { echo "ERROR: INPUT not found: $INPUT"; exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

apptainer exec -e -B /gpfs01 \
  "$SIF" \
  root -l -b <<EOF
gSystem->Load("libRooUnfold");
.x ${MACRO}+("${INPUT}","${OUT_DIR}","${HNAME}");
.q
EOF

echo "----------------------------------------"
echo "Done."
echo "Outputs:"
echo "  ROOT: $OUT_DIR/rcp.root"
echo "  PDF : $OUT_DIR/pdf/"
echo "----------------------------------------"