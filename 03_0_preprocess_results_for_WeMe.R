# =============================================================================
# Script  : 03_0_preprocess_results_for_WeMe.R
# Project : SCANDARE-ORL (HNSCC intra-tumour heterogeneity study)
# Purpose : Harmonise the outputs of the four clonal-deconvolution tools
#           (DPClust, PhylogicNDT, PyClone-VI, CloneSig) into the per-tool
#           subclonal-structure tables expected by WeMe (Salcedo et al.,
#           Nat Biotechnol 2020), and write them into the directory layout
#           that WeMe requires.
# Inputs  : - DPClust master TSV (gives the canonical sample list).
#           - Per-tool result directories under `ith_results_path`:
#               dpclust/<sample>/<sample>_DPoutput_*_bestClusterInfo.txt
#               phylogicndt/<sample>/<sample>.mut_ccfs.txt
#               pyclone_vi/<sample>_pyclonevi.tsv_output.tsv
#               clonesig/<sample>_clonesig_out.tsv
# Outputs : `weme_input_dir/<sample>/<method>_results/<sample>_subclonal_structure.txt`
#           with three columns: cluster | n_ssms | proportion
#           Proportion is the cancer cell fraction (CCF), NOT cellular
#           prevalence (CP).
# Usage   : Edit the CONFIG block then `Rscript 03_0_preprocess_results_for_WeMe.R`
# Depends : tidyverse, rjson, furrr, the WeMe R sources
#           (https://github.com/morrislab/weme)
# =============================================================================

# Libraries -------------------------------------------------------------------

library(tidyverse)
library(rjson)
library(furrr)

# CONFIG ----------------------------------------------------------------------

# DPClust master file (also serves as the canonical list of samples).
sample_list_path <- "<path/to/ith_inputs/DPClust_sample_master_file.txt>"

# Root directory holding per-tool output sub-directories.
ith_results_path <- "<path/to/results/ith_dna/>"

# Output root: per-sample, per-method input files for WeMe will be created
# under this directory.
weme_input_dir   <- "<path/to/data/ith_consensus_input/>"

# WeMe source file (provides find_sids(), genconsensus(), ...).
weme_source_path <- "<path/to/bin/weme/weme.R>"

source(weme_source_path)

# Read sample list ------------------------------------------------------------

sample_list <- read_tsv(sample_list_path)

# Per-tool result tables ------------------------------------------------------
#
# Each block returns a list-column of data.frames with the WeMe-compatible
# schema:  cluster | n_ssms | proportion (where proportion = CCF in [0, 1]).

## DPClust --------------------------------------------------------------------

DPClust_results <-
  sample_list %>%
  select(sample) %>%
  mutate(
    dpclust_results = map(sample, \(x) {
      read_tsv(
        str_c(
          ith_results_path,
          "dpclust/", "/", x, "/",
          x, "_DPoutput_2000iters_1000burnin_seed123/",
          x, "_2000iters_1000burnin_bestClusterInfo.txt"
        )
      ) %>%
        rename(
          "cluster"    = cluster.no,
          "proportion" = location,
          "n_ssms"     = no.of.mutations
        ) %>%
        # CCFs above 1 (truncation artefact) are clipped.
        mutate(proportion = ifelse(proportion > 1, 1, proportion)) %>%
        group_by(proportion) %>%
        summarise(n_ssms = sum(n_ssms, na.rm = TRUE)) %>%
        arrange(desc(proportion)) %>%
        mutate(cluster = row_number()) %>%
        select(cluster, n_ssms, proportion)
    })
  )

## PhylogicNDT ----------------------------------------------------------------

PhylogicNTD_results <-
  sample_list %>%
  select(sample) %>%
  mutate(
    phylogicndt_results = map(sample, \(x) {
      tryCatch(
        error = function(e) e,
        read_tsv(
          str_c(ith_results_path, "phylogicndt/", "/", x, "/",
                x, ".mut_ccfs.txt"),
          col_select = c("Cluster_Assignment", "clust_ccf_mean")
        ) %>%
          group_by(Cluster_Assignment) %>%
          mutate(n_ssms = n()) %>%
          distinct() %>%
          rename("cluster"    = Cluster_Assignment,
                 "proportion" = clust_ccf_mean)
      )
    })
  )

## PyClone-VI -----------------------------------------------------------------

Pyclonebeta_results <-
  sample_list %>%
  select(sample) %>%
  mutate(
    pyclone_results = map(sample, \(x) {
      tryCatch(
        error = function(e) e,
        read_tsv(
          str_c(ith_results_path, "pyclone_vi/", x,
                "_pyclonevi.tsv_output.tsv"),
          col_select = c("cluster_id", "cellular_prevalence")
        ) %>%
          group_by(cluster_id) %>%
          mutate(n_ssms = n()) %>%
          distinct() %>%
          rename("cluster"    = cluster_id,
                 "proportion" = cellular_prevalence) %>%
          # PyClone-VI cluster ids are 0-indexed; shift to 1-indexed for WeMe.
          mutate(cluster = cluster + 1)
      )
    })
  )

## CloneSig -------------------------------------------------------------------

Clonesig_results <-
  sample_list %>%
  select(sample) %>%
  mutate(
    clonesig_results = map(sample, \(x) {
      tryCatch(
        error = function(e) e,
        read_tsv(
          str_c(ith_results_path, "clonesig/", x, "_clonesig_out.tsv"),
          col_select = c("clone", "clone_ccf")
        ) %>%
          group_by(clone) %>%
          mutate(n_ssms = n(), clone = clone + 1) %>%
          distinct() %>%
          rename("cluster"    = clone,
                 "proportion" = clone_ccf)
      )
    })
  )

# Merge and write -------------------------------------------------------------

ith_results <-
  DPClust_results %>%
  left_join(PhylogicNTD_results) %>%
  left_join(Pyclonebeta_results) %>%
  left_join(Clonesig_results)

methods <- names(ith_results)[-1]

# Create the {sample}/{method} sub-directory layout expected by WeMe.
map(
  unlist(map(ith_results$sample, \(x) {
    str_c(weme_input_dir, "/", x, "/", methods)
  })),
  \(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)
)

# Helper: write only if the per-sample result is a data.frame (some tools may
# have failed on some samples; we skip those quietly but record a warning).
safe_write_subclonal <- function(result, sample_name, subdir) {
  if (is.data.frame(result)) {
    write_tsv(
      result,
      file.path(weme_input_dir, sample_name, subdir,
                str_c(sample_name, "_subclonal_structure.txt"))
    )
  } else {
    warning(sprintf("No result for sample '%s' / %s; skipping.",
                    sample_name, subdir))
  }
}

walk2(ith_results$dpclust_results,     ith_results$sample,
      ~ safe_write_subclonal(.x, .y, "dpclust_results"))

walk2(ith_results$phylogicndt_results, ith_results$sample,
      ~ safe_write_subclonal(.x, .y, "phylogicndt_results"))

walk2(ith_results$pyclone_results,     ith_results$sample,
      ~ safe_write_subclonal(.x, .y, "pyclone_results"))

walk2(ith_results$clonesig_results,    ith_results$sample,
      ~ safe_write_subclonal(.x, .y, "clonesig_results"))
