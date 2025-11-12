#!/usr/bin/env bash

filename=/gpfs01/star/pwg/svomich/JetsTrees/trees/merged_all/embedding_merged.root
outdir=/gpfs01/star/pwg/svomich/JetsTrees/analysis/unfolding/out

root -l -b -q \
  -e 'gSystem->Load("libRooUnfold");' \
  "unfold.cxx+(\"$filename\",\"$outdir\")"

