#!/usr/bin/env bash
set -euo pipefail

########################
# Configuration
########################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SIF="${SCRIPT_DIR}/roounfold.sif"
MACRO="${SCRIPT_DIR}/unfold_data.cxx"

########################
# Interactive Menu
########################

# Only show the menu if METHOD was not passed as argument
if [[ $# -ge 1 && ( "$1" == "BAYES" || "$1" == "SVD" ) ]]; then
  METHOD="$1"
  shift
else
  echo ""
  echo "╔══════════════════════════════════════════╗"
  echo "║        Jet Unfolding — Method Select     ║"
  echo "╠══════════════════════════════════════════╣"
  echo "║                                          ║"
  echo "║   1)  Bayesian (iterative)               ║"
  echo "║   2)  SVD                                ║"
  echo "║                                          ║"
  echo "╚══════════════════════════════════════════╝"
  echo ""

  while true; do
    read -rp "   Select method [1/2]: " CHOICE
    case "$CHOICE" in
      1) METHOD="BAYES"; break ;;
      2) METHOD="SVD";   break ;;
      *) echo "   Invalid choice — please enter 1 or 2." ;;
    esac
  done

  echo ""
  echo "   Selected: ${METHOD}"
  echo ""
fi

########################
# Arguments
########################

# Remaining positional args (after optional METHOD arg was consumed/shifted)
# $1: data input, $2: response file, $3: output dir, $4: nIter

if [[ $# -ge 1 ]]; then
  if [[ "$1" = /* ]]; then
    INPUT="$1"
  else
    INPUT="${BASE}/trees/$1"
  fi
else
  INPUT="${BASE}/trees/data_merged.root"
fi

RESP_FILE="${2:-${SCRIPT_DIR}/out_embedding_${METHOD}/responses_embedding.root}"
OUT_DIR="${3:-${SCRIPT_DIR}/out_data_${METHOD}}"
NITER="${4:-4}"

########################
# Checks
########################

echo "----------------------------------------"
echo "Running unfolding on REAL DATA"
echo "SCRIPT_DIR  : $SCRIPT_DIR"
echo "BASE        : $BASE"
echo "SIF         : $SIF"
echo "Macro       : $MACRO"
echo "Input data  : $INPUT"
echo "Resp. file  : $RESP_FILE"
echo "Output dir  : $OUT_DIR"
echo "Iterations  : $NITER"
echo "Method      : $METHOD"
echo "----------------------------------------"

[[ -f "$SIF"       ]] || { echo "ERROR: SIF not found:        $SIF";       exit 1; }
[[ -f "$MACRO"     ]] || { echo "ERROR: MACRO not found:      $MACRO";     exit 1; }
[[ -f "$INPUT"     ]] || { echo "ERROR: Input not found:      $INPUT";     exit 1; }
[[ -f "$RESP_FILE" ]] || { echo "ERROR: Resp. file not found: $RESP_FILE"; exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

apptainer exec -e -B /gpfs01 \
  "$SIF" \
  root -l -b <<EOF
gSystem->Load("libRooUnfold");
.x ${MACRO}+("${INPUT}","${RESP_FILE}","${OUT_DIR}","${METHOD}",${NITER});
.q
EOF

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"