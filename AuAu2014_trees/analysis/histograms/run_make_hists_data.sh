#!/bin/bash

# Stop on errors
set -e

# Directory where *this* script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default file names (can be overridden via arguments)
IN_BASENAME="${1:-data_merged.root}"
OUT_BASENAME="${2:-hists_data.root}"

# Paths relative to the script location
INPUT="${SCRIPT_DIR}/../../trees/${IN_BASENAME}"
OUTPUT="${SCRIPT_DIR}/${OUT_BASENAME}"
MACRO="${SCRIPT_DIR}/make_hists_data.C"

echo "----------------------------------------"
echo "Running DATA histogram production"
echo "Script dir : $SCRIPT_DIR"
echo "Input      : $INPUT"
echo "Output     : $OUTPUT"
echo "Macro      : $MACRO"
echo "----------------------------------------"

root -l -b -q "$MACRO(\"$INPUT\",\"$OUTPUT\")"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"
