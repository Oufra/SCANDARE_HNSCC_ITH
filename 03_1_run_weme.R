# =============================================================================
# Script  : 03_1_run_weme.R
# Project : SCANDARE-ORL (HNSCC intra-tumour heterogeneity study)
# Purpose : Run the WeMe consensus algorithm (Salcedo et al., Nat Biotechnol
#           2020) on the per-tool subclonal-structure files prepared by
#           03_0_preprocess_results_for_WeMe.R.
# Inputs  : Per-sample, per-method `_subclonal_structure.txt` files arranged
#           in the directory layout that WeMe expects:
#             <weme_working_dir>/<sample>/<method>_results/<sample>_subclonal_structure.txt
# Outputs : One consensus subclonal-structure file per sample (`genconsensus`
#           writes to the working directory).
# Usage   : 1. setwd() to the WeMe input directory (the script MUST be run
#              from there, see comment below).
#           2. Edit `weme_source_path` to point at weme.R.
#           3. Rscript 03_1_run_weme.R
# Notes   : `find_sids()` discovers samples by listing sub-directories of the
#           CURRENT WORKING DIRECTORY, so this script must be executed from
#           inside the WeMe working directory.
# Depends : WeMe R sources (https://github.com/morrislab/weme)
# =============================================================================

# CONFIG ----------------------------------------------------------------------

weme_source_path <- "<path/to/bin/weme/weme.R>"

# ----------------------------------------------------------------------------

source(weme_source_path)

sids <- find_sids()

genconsensus(sids, rounddown = FALSE)
