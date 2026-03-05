#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 1st arg: input ROOT file basename in ../../trees/ (or absolute path)
IN_ARG="${1:-data_merged.root}"
if [[ "$IN_ARG" = /* ]]; then
  INPUT="$IN_ARG"
else
  INPUT="${SCRIPT_DIR}/../../trees/${IN_ARG}"
fi

MACRO="${SCRIPT_DIR}/make_qa_pdfs.C"

# Where PDFs will be created (macro uses relative paths pdf/QA, pdf/Jets)
OUTDIR="${2:-$SCRIPT_DIR}"
mkdir -p "$OUTDIR"

echo "----------------------------------------"
echo "Running QA PDF production"
echo "Script dir : $SCRIPT_DIR"
echo "Input      : $INPUT"
echo "Macro      : $MACRO"
echo "Output dir : $OUTDIR"
echo "----------------------------------------"

# Run from OUTDIR so pdf/QA and pdf/Jets land there (not wherever you launched from)
pushd "$OUTDIR" >/dev/null

root -l -b -q "${MACRO}(\"${INPUT}\")"

popd >/dev/null

echo "----------------------------------------"
echo "Done. PDFs are in: ${OUTDIR}/pdf/QA and ${OUTDIR}/pdf/Jets"
echo "----------------------------------------"

