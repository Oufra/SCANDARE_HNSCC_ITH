#!/bin/bash
# =============================================================================
# Script  : 02_3_run_pyclone-vi.sh
# Project : SCANDARE-ORL (HNSCC intra-tumour heterogeneity study)
# Purpose : Fit PyClone-VI (Gillis & Roth, BMC Bioinformatics 2020) on one
#           pre-processed per-sample variant table and write the per-mutation
#           cluster assignments.
# Inputs  : -i  Pre-processed PyClone-VI input TSV produced by
#               01_preprocess_ITH_4tools.R. Columns:
#                 mutation_id, sample_id, ref_counts, alt_counts,
#                 major_cn, minor_cn, normal_cn, tumour_content
#           -o  Output directory (will be created if missing)
#           -d  Statistical model passed to `pyclone-vi fit -d`
#               (default: beta-binomial; alternative: binomial)
# Outputs : - <sample>_output.h5 : raw PyClone-VI fit
#           - <sample>_output.tsv: per-mutation cluster table
# Run     : pyclone-vi fit  -c 10 (max clusters) -r 100 (restarts)
# Usage   : ./02_3_run_pyclone-vi.sh -i <input.tsv> -o <outdir/> [-d <model>]
# Depends : pyclone-vi (https://github.com/Roth-Lab/pyclone-vi)
# =============================================================================

set -euo pipefail

# Default model
D="beta-binomial"

while getopts i:o:hd: flag
do
    case "${flag}" in
        h)
            echo "Usage:

  ./02_3_run_pyclone-vi.sh -i <pyclone-vi_input.tsv> -o <outdir/> -d [beta-binomial|binomial]

  -i  Pre-processed PyClone-VI input table
  -o  Output directory
  -d  Statistical model (default: beta-binomial)
"
            exit 0 ;;
        i) INPUT=${OPTARG} ;;
        o) OUTPUT_DIR=${OPTARG} ;;
        d) D=${OPTARG} ;;
        *) echo "Invalid option: -$flag" >&2
           exit 1 ;;
    esac
done

mkdir -p "${OUTPUT_DIR}"

OUTPUT_H5=${OUTPUT_DIR}/$(basename "${INPUT}" _input.tsv)_output.h5
OUTPUT_TSV=${OUTPUT_DIR}/$(basename "${INPUT}" _input.tsv)_output.tsv

pyclone-vi fit \
    -i "${INPUT}" \
    -o "${OUTPUT_H5}" \
    -c 10 \
    -r 100 \
    -d "${D}"

pyclone-vi write-results-file \
    -i "${OUTPUT_H5}" \
    -o "${OUTPUT_TSV}"
