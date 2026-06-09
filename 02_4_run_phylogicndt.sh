#!/bin/bash
# =============================================================================
# Script  : 02_4_run_phylogicndt.sh
# Project : SCANDARE-ORL (HNSCC intra-tumour heterogeneity study)
# Purpose : Submit one PhylogicNDT Cluster job per region to a PBS/Torque
#           cluster via qsub. Uses a Singularity container that ships the
#           PhylogicNDT codebase.
# Inputs  : - One <sample>_phylogicndt.sif file per sample (produced by
#             01_preprocess_ITH_4tools.R) in ${ITH_INPUTS_DIR}.
# Outputs : - One sub-directory per sample under ${PHYLOGIC_OUT} containing
#             PhylogicNDT Cluster outputs (.mut_ccfs.txt, .cluster_ccfs.txt,
#             etc.).
# Usage   : Edit the CONFIG block then run `./02_4_run_phylogicndt.sh`.
#           The script is idempotent at the qsub level: rerunning will resubmit
#           jobs (rename / delete the output sub-directories first if needed).
# Depends : singularity, qsub (Torque / PBS), PhylogicNDT image
# =============================================================================

set -euo pipefail

# --- CONFIG ------------------------------------------------------------------
PHYLOGIC_OUT="<path/to/results/ith_dna/phylogicndt>"
PHYLOGIC_SINGULARITY="<path/to/singularity_images/phylogicndt.sif>"
PHYLOGIC_SCRIPT="<path/to/phylogicndt/PhylogicNDT.py>"
PHYLOGIC_GENE_LIST="<path/to/phylogicndt/data/supplement_data/census_HGNC.txt>"
ITH_INPUTS_DIR="<path/to/data/ith_inputs>"
LOGS_DIR="<path/to/logs>"
# -----------------------------------------------------------------------------

mkdir -p "${PHYLOGIC_OUT}"
mkdir -p "${LOGS_DIR}"

for SIF in $(find "${ITH_INPUTS_DIR}" -name "*_phylogicndt.sif")
do
    SAMPLE_NAME=$(basename "${SIF}" _phylogicndt.sif)
    mkdir -p "${PHYLOGIC_OUT}/${SAMPLE_NAME}/"
    PARAM=""

    echo "singularity exec ${PHYLOGIC_SINGULARITY} ${PHYLOGIC_SCRIPT} Cluster \
        -i ${SAMPLE_NAME} \
        -sif ${SIF} \
        --maf_input_type calc_ccf \
        --impute \
        --driver_genes_file ${PHYLOGIC_GENE_LIST} \
        -rb ${PARAM}" | \
    qsub -N phylogicNDT_${SAMPLE_NAME} \
         -e ${LOGS_DIR} \
         -d "${PHYLOGIC_OUT}/${SAMPLE_NAME}/" \
         -o ${LOGS_DIR} \
         -l nodes=1:ppn=8,mem=8g,walltime=72:00:00
done
