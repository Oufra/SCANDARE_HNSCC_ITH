# =============================================================================
# Script  : 01_preprocess_ITH_4tools.R
# Project : SCANDARE-ORL (HNSCC intra-tumour heterogeneity study)
# Purpose : Build tool-specific input tables for four clonal-deconvolution
#           algorithms from DRAGON SNV reports and FACETS CNV calls:
#             - CloneSig       (Abecassis et al.)
#             - PyClone-VI     (Roth et al.)
#             - PhylogicNDT    (Leshchiner et al.)
#             - DPClust        (Nik-Zainal / Sanger)
# Inputs  : - Sample list TSV
#           - Per-run DRAGON SNV reports (`variant_report.table1.tsv`)
#           - Per-run FACETS CNV reports (`cnv_report.table2.tsv`)
#           - hg19 reference FASTA (for trinucleotide context, CloneSig)
#           - OncoKB cancer gene list (for the recurrence-rescue rule)
# Outputs : One sub-directory per region under `local_output_ithinput_dir`,
#           each containing:
#             {Barcode}_clonesig.tsv
#             {Barcode}_pyclonevi.tsv
#             {Barcode}_phylogicndt.tsv
#             {Barcode}_phylogicndt.sif
#             DPClust input files: _loci.txt, _subclones.txt,
#               _alleleFrequencies.txt, _rho_and_psi.txt,
#               _allDirichletProcessInfo.txt
#           Plus a master file `DPClust_sample_master_file.txt`.
# Usage   : Edit the CONFIG block, then `Rscript 01_preprocess_ITH_4tools.R`
#           (or run interactively). Parallel computing uses `furrr` with 9
#           workers — adjust to your machine.
# =============================================================================

# Libraries -------------------------------------------------------------------

library(GenomicRanges)
library(Rsamtools)
library(plyranges)
library(VariantAnnotation)
library(tidyverse)
library(vroom)
library(janitor)
library(ggpubr)
library(ggfortify)   # for autoplot() on prcomp objects
library(dpclust3p)
library(furrr)

# Parallel back-end
plan("multisession", workers = 9)

# CONFIG ----------------------------------------------------------------------

# hg19 reference FASTA (must be indexed: .fai expected alongside).
hg19_fa_path             <- "<path/to/hg19.fa>"

# Path on the HPC cluster where ITH input files will be read from by downstream
# jobs (used to write absolute paths into PhylogicNDT .sif files).
ith_input_dir_cluster    <- "<cluster/path/to/ith_inputs>"

# Local directory where this script writes the per-sample input files.
local_output_ithinput_dir <- "<path/to/ith_inputs>"

# OncoKB cancer gene list.
oncokb_gene_list_path    <- "<path/to/cancerGeneList.tsv>"

# Sample list TSV (must contain ID_SCANDARE and Barcode_DRAGON columns).
SCANDARE_sample_list_path <- "<path/to/SCANDARE_ORL_sample_list.tsv>"

# Root directory containing per-run DRAGON output folders.
SCANDARE_data_path        <- "<path/to/dragon_runs>"

# Resources -------------------------------------------------------------------

hg19.fa <- FaFile(hg19_fa_path)

# Data import -----------------------------------------------------------------

## OncoKB gene list ----------------------------------------------------------

oncokb_gene_list.tsv <-
  read_tsv(oncokb_gene_list_path) %>%
  clean_names() %>%
  separate_rows(gene_aliases)

oncokb_all_genes <-
  bind_rows(
    oncokb_gene_list.tsv %>%
      select("Gene" = hugo_symbol, is_oncogene, is_tumor_suppressor_gene),
    oncokb_gene_list.tsv %>%
      select("Gene" = gene_aliases, is_oncogene, is_tumor_suppressor_gene)
  ) %>%
  distinct()

oncokb_oncogene_list <-
  oncokb_all_genes %>%
  filter(is_oncogene == "Yes") %>%
  pull(Gene)

## Sample metadata ------------------------------------------------------------

SCANDARE_sample_list.tsv <- read_tsv(SCANDARE_sample_list_path)

# CONIPHER requires patient IDs to share a common string prefix followed by a
# number. We build `CASE_ID = "SCDR" + numeric part of ID_SCANDARE`.
SCANDARE_sample_list <-
  SCANDARE_sample_list.tsv %>%
  separate_rows(Barcode_DRAGON, sep = "/") %>%
  mutate(
    id_run_DRAGON = str_extract(Barcode_DRAGON, "D[:digit:]{3}"),
    CASE_ID       = str_c("SCDR", str_remove(ID_SCANDARE, "-"))
  ) %>%
  group_by(CASE_ID) %>%
  mutate(SAMPLE = str_c(CASE_ID, "_", row_number())) %>%
  ungroup()

## File index ----------------------------------------------------------------

# One row per DRAGON run, with paths to the SNV and CNV reports.
snv_files_table <-
  tibble(
    snv_reports_paths = list.files(
      SCANDARE_data_path,
      pattern    = "variant_report.table1.tsv",
      recursive  = TRUE,
      full.names = TRUE
    ),
    id_run_DRAGON = str_extract(snv_reports_paths, "D[:digit:]{3}")
  )

cnv_files_table <-
  tibble(
    cnv_reports_paths = list.files(
      SCANDARE_data_path,
      pattern    = "cnv_report.table2.tsv",
      recursive  = TRUE,
      full.names = TRUE
    ),
    id_run_DRAGON = str_extract(cnv_reports_paths, "D[:digit:]{3}")
  )

files_table <-
  full_join(snv_files_table, cnv_files_table) %>%
  na.omit()

# SNV data import and filtering -----------------------------------------------

snv_columns_to_keep <- c(
  "Barcode", "Gene", "C_point", "P_point", "Type",
  "Chr", "Start", "End", "Ref", "Alt",
  "Cov_Recalc_RF_S1", "Cov_Recalc_RR_S1",
  "Cov_Recalc_AF_S1", "Cov_Recalc_AR_S1",
  "TAG_BIOINFO"
)

# DRAGON in-house tags flagging likely polymorphisms / recurrent artefacts.
TAG_BIOINFO_to_exclude <- c(
  "frameshift recurrent",
  "polym rare decrit faible desequilibre",
  "polym rare decrit faible equilibre",
  "polym rare decrit fort",
  "polym rare non decrit faible",
  "polym rare non decrit fort desequilibre",
  "polym rare non decrit fort equilibre",
  "polymorphisme",
  "redondant faible"
)

# Read all SNV tables. This is the longest single step in the script.
snv_data <-
  files_table %>%
  select(id_run_DRAGON, snv_reports_paths) %>%
  mutate(
    snv = map(
      snv_reports_paths,
      ~ {
        vroom(
          .x,
          col_select = c(all_of(snv_columns_to_keep), starts_with("Freq"))
        ) %>%
          mutate(across(starts_with("Freq"), \(x) {
            as.numeric(str_replace(x, ",", "."))
          })) %>%
          filter(
            # Keep only SCANDARE samples.
            Barcode %in% SCANDARE_sample_list$Barcode_DRAGON &
            # TODO: the filter `!TAG_BIOINFO %in% TAG_BIOINFO_to_exclude` was
            # disabled in the analysis used for the manuscript. Re-enable
            # by uncommenting if you want TAG_BIOINFO-based exclusion.
            # !TAG_BIOINFO %in% TAG_BIOINFO_to_exclude &
            # Remove variants without coverage information.
            if_all(starts_with("Cov"), \(x) !is.na(x)) &
            # Remove variants with 0 reads in forward or reverse.
            !if_any(c(Cov_Recalc_AF_S1, Cov_Recalc_AR_S1), \(x) x == 0) &
            # Remove common polymorphisms (Freq* < 0.001 across populations).
            if_all(starts_with("Freq"), \(x) is.na(x) | x < 0.001)
          )
      }
    )
  ) %>%
  unnest(snv) %>%
  rowwise() %>%
  mutate(
    DEPTH     = sum(Cov_Recalc_RF_S1, Cov_Recalc_RR_S1,
                    Cov_Recalc_AF_S1, Cov_Recalc_AR_S1),
    ALT_DEPTH = sum(Cov_Recalc_AF_S1, Cov_Recalc_AR_S1),
    REF_DEPTH = sum(Cov_Recalc_RF_S1, Cov_Recalc_RR_S1),
    alt_VAF   = ALT_DEPTH / (ALT_DEPTH + REF_DEPTH)
  ) %>%
  mutate(
    pass_filter_depth = DEPTH >= 100,
    pass_filter_VAF   = alt_VAF >= 0.1
  ) %>%
  left_join(SCANDARE_sample_list, by = c("Barcode" = "Barcode_DRAGON")) %>%
  group_by(Chr, Start, End) %>%
  # Recurrence filter: remove variants seen in >=5 patients (likely artefacts),
  # but rescue oncogene variants seen in <30 patients with a meaningful coding
  # change, and always rescue TERT promoter variants.
  mutate(
    pass_filter_recurrence =
      n_distinct(ID_SCANDARE) < 5 |
      Gene %in% oncokb_oncogene_list &
        n_distinct(ID_SCANDARE) < 30 &
        Type %in% c(
          "in-frame_deletion", "in-frame_insertion",
          "nonframeshift_deletion", "nonframeshift_insertion",
          "nonsynonymous_SNV", "SNV", "synonymous_SNV"
        ) |
      Gene == "TERT"
  ) %>%
  # IMPORTANT: ungroup before furrr / future_pmap, otherwise extreme slowdowns.
  # See https://furrr.futureverse.org/articles/gotchas.html
  ungroup() %>%
  mutate(
    # Strand-bias test (Fisher exact on the 2x2 strand x allele table).
    # Phred-scaled p-value; thresholds borrowed from GATK FilterMutectCalls.
    phred_bias = -10 * log10(future_pmap_dbl(
      .progress = TRUE,
      list(Cov_Recalc_RF_S1, Cov_Recalc_RR_S1,
           Cov_Recalc_AF_S1, Cov_Recalc_AR_S1),
      ~ broom::tidy(fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2)))$p.value
    )),
    MUT_TYPE = ifelse(nchar(Ref) > 1 | nchar(Alt) > 1, "INDEL", "SNV")
  ) %>%
  mutate(
    pass_filter_bias =
      MUT_TYPE == "INDEL" & phred_bias <= 200 |
      MUT_TYPE == "SNV"   & phred_bias <  50
  )

# SNV exploratory plots -------------------------------------------------------

snv_filter_plot_data <-
  snv_data %>%
  mutate(
    filter = case_when(
      !pass_filter_recurrence ~ "Recurrence",
      !pass_filter_depth      ~ "Depth < 100X",
      !pass_filter_VAF        ~ "VAF < 0.1",
      !pass_filter_bias       ~ "Strand bias",
      TRUE                    ~ "Pass",
      .ptype = factor(
        levels = c("Recurrence", "Depth < 100X", "VAF < 0.1",
                   "Strand bias", "Pass"),
        ordered = TRUE
      )
    )
  )

snv_filter_plot <-
  snv_filter_plot_data %>%
  group_by(Barcode) %>%
  mutate(n_variants = dplyr::n()) %>%
  ggplot(aes(x = reorder(Barcode, -n_variants), fill = filter)) +
  geom_bar() +
  theme_bw() +
  theme(
    axis.title.x        = element_blank(),
    axis.text.x         = element_blank(),
    axis.ticks.x        = element_blank(),
    panel.grid.major.x  = element_blank()
  ) +
  labs(y = "Number of variants", fill = "Filters",
       title = "Filtered variants") +
  scale_y_continuous(expand = c(0, 0))

## Shared vs. private variants per patient

snv_in_common <-
  snv_filter_plot_data %>%
  group_by(CASE_ID) %>%
  filter(n_distinct(Barcode) > 1) %>%
  group_by(CASE_ID, Chr, Start, End, Ref, Alt) %>%
  mutate(unique = ifelse(n_distinct(Barcode) == 1,
                         "Private variant", "Shared variant"))

snv_in_common_plot <-
  snv_in_common %>%
  ggplot(aes(x = unique, fill = filter)) +
  geom_bar() +
  theme_bw() +
  theme(
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_blank()
  ) +
  facet_wrap(~unique, scales = "free") +
  labs(x = "", y = "")

## Functional impact of filtered vs. retained variants in driver genes

snv_filtered_functional_impact <-
  snv_filter_plot_data %>%
  mutate(
    functional_impact = case_when(
      Type %in% c("frameshif_deletion", "frameshift_insertion",
                  "stopgain") ~ "LoF",
      Type %in% c("SNV", "nonsynonymous_SNV",
                  "in-frame_deletion", "in-frame_insertion")
                              ~ "Missense/in-frame",
      TRUE                    ~ "Other"
    ),
    filtered = ifelse(filter == "Pass", "Pass", "Removed")
  ) %>%
  left_join(oncokb_all_genes) %>%
  mutate(
    gene_function = ifelse(
      is_oncogene == "Yes" & is_tumor_suppressor_gene == "Yes",
      "Oncogene & TSG",
      ifelse(
        is_oncogene == "Yes", "Oncogene",
        ifelse(is_tumor_suppressor_gene == "Yes", "TSG", "Other")
      )
    )
  ) %>%
  filter(gene_function %in% c("Oncogene", "TSG"))

snv_filtered_functional_impact_plot <-
  snv_filtered_functional_impact %>%
  ggplot(aes(x = gene_function, fill = functional_impact)) +
  geom_bar() +
  theme_bw() +
  theme(
    axis.title.x       = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(y = "Number of variants", fill = "Functional impact",
       title = "Functional impact of filtered variants") +
  facet_wrap(~filtered, scales = "free")

snv_filtered_functional_impact_plot
snv_in_common_plot

snv_pass_per_sample_plot <-
  snv_filter_plot_data %>%
  filter(filter == "Pass") %>%
  group_by(Barcode) %>%
  tally() %>%
  ggplot(aes(x = reorder(Barcode, -n), y = n)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.title.x       = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(y = "Number of variants", fill = "Filters",
       title = "Filtered variants") +
  scale_y_continuous(expand = c(0, 0))
snv_pass_per_sample_plot

snv_exploratory_graphs <-
  ggarrange(snv_filter_plot, snv_in_common_plot, common.legend = TRUE)
snv_exploratory_graphs

# CNV data import -------------------------------------------------------------

cnv_columns_to_keep <- c(
  "Barcode", "Chr", "Start", "Stop",
  "BAF",          # encoded as e.g. "AABB" -> derives COPY_NUMBER_A/_B below
  "Ploidy", "Cellularity"
)

cnv_data <-
  files_table %>%
  select(id_run_DRAGON, cnv_reports_paths) %>%
  mutate(
    cnv = map(cnv_reports_paths, ~ {
      vroom(.x, col_select = all_of(cnv_columns_to_keep)) %>%
        filter(Barcode %in% SCANDARE_sample_list$Barcode_DRAGON)
    })
  ) %>%
  unnest(cnv) %>%
  filter(!is.na(Cellularity)) %>%
  mutate(
    COPY_NUMBER_A = str_count(BAF, "A"),
    COPY_NUMBER_B = str_count(BAF, "B")
  )

# Join SNV and CNV data as GRanges -------------------------------------------

snv_grange <-
  snv_data %>%
  filter(Chr != "chrX") %>%
  filter(pass_filter_recurrence, pass_filter_depth,
         pass_filter_VAF, pass_filter_bias) %>%
  rowwise() %>%
  select(-c(id_run_DRAGON.x, id_run_DRAGON.y, snv_reports_paths,
            phred_bias, TAG_BIOINFO, matches("Cov"))) %>%
  group_by(Barcode) %>%
  nest(.key = "snv_table") %>%
  mutate(snv_grange = map(snv_table,
                          ~ makeGRangesFromDataFrame(.x,
                                                     keep.extra.columns = TRUE)))

cnv_grange <-
  cnv_data %>%
  filter(Chr != "chr23") %>%
  filter(!is.na(COPY_NUMBER_A) & !is.na(COPY_NUMBER_B)) %>%
  select(-id_run_DRAGON, -cnv_reports_paths) %>%
  group_by(Barcode) %>%
  nest(.key = "cnv_table") %>%
  mutate(cnv_grange = map(cnv_table,
                          ~ makeGRangesFromDataFrame(.x,
                                                     keep.extra.columns = TRUE)))

snv_with_cnv_data_grange <-
  inner_join(snv_grange, cnv_grange, by = "Barcode") %>%
  select(Barcode, snv_grange, cnv_grange) %>%
  filter(Barcode %in% SCANDARE_sample_list$Barcode_DRAGON &
         Barcode != "D700R19") %>%
  mutate(
    snv_cnv_grange = map2(snv_grange, cnv_grange, \(x, y) {
      join_overlap_inner(x, y)
    })
  ) %>%
  select(Barcode, snv_cnv_grange) %>%
  mutate(snv_cnv_grange = map2(snv_cnv_grange, Barcode, \(x, y) {
    x$Barcode <- y
    x
  }))

# Concatenate all per-sample GRanges into one GRanges.
snv_with_cnv_data_grange_unlisted <-
  unlist(as(snv_with_cnv_data_grange$snv_cnv_grange, "GRangesList"))

# Output directories ---------------------------------------------------------

map(str_c(local_output_ithinput_dir, "/", unique(SCANDARE_sample_list$Barcode_DRAGON)),
    \(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

# CloneSig input --------------------------------------------------------------

# CloneSig only handles SNVs.
snv_with_cnv_data_grange_unlisted_noindel <-
  snv_with_cnv_data_grange_unlisted %>%
  filter(MUT_TYPE == "SNV")

# Generate the trinucleotide context required by CloneSig to detect
# mutational signatures. Reference and alternative bases are pyrimidine-
# normalised (C and T as reference).
generate_pattern <- function(context, alt) {
  alt        <- unlist(alt)
  context    <- as.character(context)
  ref        <- substr(context, 2, 2)
  match_dict <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
  if (ref %in% c("C", "T")) {
    pattern <- str_c(substr(context, 1, 1), "[", ref, ">", alt, "]",
                     substr(context, 3, 3))
  } else {
    pattern <- str_c(
      match_dict[substr(context, 3, 3)], "[",
      match_dict[ref], ">", match_dict[alt], "]",
      match_dict[substr(context, 1, 1)]
    )
  }
  pattern
}

context_snv_noindel <-
  scanFa(hg19.fa, stretch(snv_with_cnv_data_grange_unlisted_noindel, 2))
alt_snv_noindel <- snv_with_cnv_data_grange_unlisted_noindel$Alt

snv_with_cnv_data_grange_unlisted_noindel$trinucleotide_context <-
  map2_chr(context_snv_noindel, alt_snv_noindel,
           \(x, y) generate_pattern(x, y))

clonesig_table <-
  snv_with_cnv_data_grange_unlisted_noindel %>%
  as.data.frame() %>%
  mutate(alteration = str_c(seqnames, ":", start, Ref, ">", Alt)) %>%
  mutate(
    total_read_count = ALT_DEPTH + REF_DEPTH,
    C_tumor_tot      = COPY_NUMBER_A + COPY_NUMBER_B
  ) %>%
  select(
    Barcode,
    alteration,
    "pattern"              = trinucleotide_context,
    "variant_allele_count" = ALT_DEPTH,
    total_read_count,
    "freq"                 = alt_VAF,
    C_tumor_tot,
    "C_tumor_minor"        = COPY_NUMBER_B,
    "purity"               = Cellularity
  ) %>%
  group_by(Barcode) %>%
  mutate(purity = mean(purity)) %>%
  nest(.key = "variant_tables")

# PhylogicNDT input -----------------------------------------------------------

phylogicndt_table <-
  snv_with_cnv_data_grange_unlisted %>%
  as.data.frame() %>%
  rowwise() %>%
  mutate(local_cn_a2 = COPY_NUMBER_B, local_cn_a1 = COPY_NUMBER_A) %>%
  mutate(Protein_change = str_c(C_point, "_", P_point)) %>%
  select(
    Barcode,
    Cellularity,
    "Chromosome"        = seqnames,
    "Start_position"    = start,
    "Reference_Allele"  = Ref,
    "Tumor_Seq_Allele2" = Alt,
    "t_ref_count"       = REF_DEPTH,
    "t_alt_count"       = ALT_DEPTH,
    local_cn_a2,
    local_cn_a1
  ) %>%
  group_by(Barcode) %>%
  mutate(Cellularity = mean(Cellularity)) %>%
  group_by(Barcode, Cellularity) %>%
  nest(.key = "variant_tables")

# Companion .sif file (sample information file) expected by PhylogicNDT.
phylogicndt_sif <-
  phylogicndt_table %>%
  ungroup() %>%
  transmute(
    Barcode,
    sample_id = Barcode,
    maf_fn    = str_c(ith_input_dir_cluster, "/", sample_id, "/",
                      sample_id, "_phylogicndt.tsv"),
    seg_fn    = "",
    purity    = Cellularity,
    timepoint = 1
  ) %>%
  group_by(Barcode) %>%
  nest(.key = "sif_tables")

# PyClone-VI input ------------------------------------------------------------

pyclonevi_table <-
  snv_with_cnv_data_grange_unlisted %>%
  as.data.frame() %>%
  mutate(
    mutation_id = str_c(seqnames, ":", start, Ref, ">", Alt),
    sample_id   = Barcode,
    normal_cn   = 2
  ) %>%
  group_by(sample_id) %>%
  mutate(tumour_content = mean(Cellularity)) %>%
  select(
    Barcode,
    mutation_id,
    sample_id,
    "ref_counts" = REF_DEPTH,
    "alt_counts" = ALT_DEPTH,
    "major_cn"   = COPY_NUMBER_A,
    "minor_cn"   = COPY_NUMBER_B,
    normal_cn,
    tumour_content
  ) %>%
  group_by(Barcode) %>%
  nest(.key = "variant_tables")

# DPClust input ---------------------------------------------------------------
#
# DPClust expects a fixed set of per-sample files. We assemble each component
# in a list-column then iterate row-wise to write them out.

snv_data_dpclust <-
  snv_data %>%
  filter(Chr != "chrX") %>%
  filter(pass_filter_recurrence, pass_filter_depth,
         pass_filter_VAF, pass_filter_bias) %>%
  mutate(Chr = str_remove(Chr, "chr")) %>%
  select(Barcode, Chr, Start, End, Ref, Alt, matches("Cov_Recalc"))

loci_dpclust <-
  snv_data_dpclust %>%
  select(Barcode, Chr, Start, Ref, Alt) %>%
  group_by(Barcode) %>%
  nest(.key = "loci")

subclones_dpclust <-
  cnv_data %>%
  filter(Chr != "chr23") %>%
  mutate(Chr = str_remove(Chr, "chr")) %>%
  select(
    Barcode,
    "chr"     = Chr,
    "startpos" = Start,
    "endpos"   = Stop,
    "nMaj1_A" = COPY_NUMBER_A,
    "nMin1_A" = COPY_NUMBER_B
  ) %>%
  mutate(frac1_A = 1, nMaj2_A = NA, nMin2_A = NA, frac2_A = NA) %>%
  group_by(Barcode) %>%
  nest(.key = "subclones")

# Per-base allele counts (Count_A/C/G/T columns) required by DPClust's
# downstream `GetWTandMutCount()`.
alleleFrequencies_dpclust <-
  snv_data_dpclust %>%
  mutate(
    Count_A = case_when(
      Ref == "A" ~ Cov_Recalc_RF_S1 + Cov_Recalc_RR_S1,
      Alt == "A" ~ Cov_Recalc_AF_S1 + Cov_Recalc_AR_S1,
      TRUE ~ 0
    ),
    Count_C = case_when(
      Ref == "C" ~ Cov_Recalc_RF_S1 + Cov_Recalc_RR_S1,
      Alt == "C" ~ Cov_Recalc_AF_S1 + Cov_Recalc_AR_S1,
      TRUE ~ 0
    ),
    Count_G = case_when(
      Ref == "G" ~ Cov_Recalc_RF_S1 + Cov_Recalc_RR_S1,
      Alt == "G" ~ Cov_Recalc_AF_S1 + Cov_Recalc_AR_S1,
      TRUE ~ 0
    ),
    Count_T = case_when(
      Ref == "T" ~ Cov_Recalc_RF_S1 + Cov_Recalc_RR_S1,
      Alt == "T" ~ Cov_Recalc_AF_S1 + Cov_Recalc_AR_S1,
      TRUE ~ 0
    )
  ) %>%
  rowwise() %>%
  mutate(
    Good_depth = sum(Cov_Recalc_RF_S1, Cov_Recalc_RR_S1,
                     Cov_Recalc_AF_S1, Cov_Recalc_AR_S1)
  ) %>%
  select(Barcode, "#CHR" = Chr, "POS" = Start,
         starts_with("Count"), Good_depth) %>%
  group_by(Barcode) %>%
  nest(.key = "alleleFrequencies")

# Purity/ploidy file expected by DPClust (3x5 matrix per sample).
rho_and_psi_dpclust <-
  cnv_data %>%
  select(Barcode, Cellularity, Ploidy) %>%
  group_by(Barcode) %>%
  mutate(Cellularity = mean(Cellularity)) %>%
  distinct() %>%
  rowwise() %>%
  mutate(
    purity_ploidy = list(
      array(NA, c(3, 5)) %>% {
        colnames(.) <- c("rho", "psi", "ploidy", "distance", "is.best")
        rownames(.) <- c("ASCAT", "FRAC_GENOME", "REF_SEG")
        .["FRAC_GENOME", "rho"] <- Cellularity
        .["FRAC_GENOME", "psi"] <- Ploidy
        .
      }
    )
  ) %>%
  select(Barcode, purity_ploidy) %>%
  group_by(Barcode) %>%
  nest(.key = "rho_and_psi")

outfile.prefix <- local_output_ithinput_dir

data_dpclust <-
  loci_dpclust %>%
  inner_join(subclones_dpclust,         by = "Barcode") %>%
  inner_join(alleleFrequencies_dpclust, by = "Barcode") %>%
  inner_join(rho_and_psi_dpclust,       by = "Barcode") %>%
  mutate(
    subclones_path         = str_c(outfile.prefix, "/", Barcode, "/",
                                   Barcode, "_subclones.txt"),
    loci_path              = str_c(outfile.prefix, "/", Barcode, "/",
                                   Barcode, "_loci.txt"),
    alleleFrequencies_path = str_c(outfile.prefix, "/", Barcode, "/",
                                   Barcode, "_alleleFrequencies.txt"),
    rho_and_psi_path       = str_c(outfile.prefix, "/", Barcode, "/",
                                   Barcode, "_rho_and_psi.txt")
  )

# Write per-sample DPClust intermediate files ---------------------------------

for (i in 1:nrow(data_dpclust)) {
  write.table(data_dpclust$loci[[i]],
              file = data_dpclust$loci_path[i],
              sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  write.table(data_dpclust$subclones[[i]],
              file = data_dpclust$subclones_path[i],
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data_dpclust$alleleFrequencies[[i]],
              file = data_dpclust$alleleFrequencies_path[i],
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data_dpclust$rho_and_psi[[i]]$purity_ploidy,
              file = data_dpclust$rho_and_psi_path[i],
              sep = "\t", quote = FALSE)
}

# Assemble DPClust dpIn file --------------------------------------------------
#
# We override `dpclust3p::GetWTandMutCount()` with a copy that fixes the
# handling of per-base allele-count columns (the upstream implementation does
# not match the column layout we produce above).

poss_runGetDirichletProcessInfo <-
  possibly(.f = runGetDirichletProcessInfo, otherwise = "Error")

GetWTandMutCount <- function(loci_file, allele_frequencies_file) {
  subs.data <- tryCatch(
    read.table(loci_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE),
    error = function(e) NA
  )
  subs.data <- subs.data[order(subs.data[, 1], subs.data[, 2]), ]
  subs.data[, 3] <- apply(as.data.frame(subs.data[, 3]), 1,
                          function(x) substring(x, 1, 1))
  subs.data[, 4] <- apply(as.data.frame(subs.data[, 4]), 1,
                          function(x) substring(x, 1, 1))

  subs.data.gr <- GenomicRanges::GRanges(
    subs.data[, 1],
    IRanges::IRanges(subs.data[, 2], subs.data[, 2]),
    rep("*", nrow(subs.data))
  )
  elementMetadata(subs.data.gr) <- subs.data[, c(3, 4)]

  alleleFrequencies <- read.delim(
    allele_frequencies_file,
    sep = "\t", header = TRUE, quote = NULL, stringsAsFactors = FALSE
  )
  alleleFrequencies <- alleleFrequencies[
    order(alleleFrequencies[, 1], alleleFrequencies[, 2]),
  ]
  print(head(alleleFrequencies))

  alleleFrequencies.gr <- GenomicRanges::GRanges(
    alleleFrequencies[, 1],
    IRanges::IRanges(alleleFrequencies[, 2], alleleFrequencies[, 2]),
    rep("*", nrow(alleleFrequencies))
  )
  elementMetadata(alleleFrequencies.gr) <- alleleFrequencies[, 3:7]

  overlap <- findOverlaps(subs.data.gr, alleleFrequencies.gr)
  alleleFrequencies <- alleleFrequencies[subjectHits(overlap), ]

  nucleotides <- c("A", "C", "G", "T")
  ref.indices <- match(subs.data[, 3], nucleotides)
  alt.indices <- match(subs.data[, 4], nucleotides)

  WT.count <- as.numeric(unlist(sapply(
    1:nrow(alleleFrequencies),
    function(v, a, i) v[i, a[i] + 2],
    v = alleleFrequencies, a = ref.indices
  )))
  mut.count <- as.numeric(unlist(sapply(
    1:nrow(alleleFrequencies),
    function(v, a, i) v[i, a[i] + 2],
    v = alleleFrequencies, a = alt.indices
  )))

  combined <- data.frame(
    chr     = subs.data[, 1],
    pos     = subs.data[, 2],
    WTCount = WT.count,
    mutCount = mut.count
  )
  colnames(combined) <- c("chr", "pos", "WT.count", "mut.count")

  combined.gr <- GenomicRanges::GRanges(
    seqnames(subs.data.gr),
    ranges(subs.data.gr),
    rep("*", nrow(subs.data))
  )
  elementMetadata(combined.gr) <- data.frame(
    WT.count = WT.count, mut.count = mut.count
  )
  combined.gr <- sortSeqlevels(combined.gr)
  combined.gr
}

assignInNamespace("GetWTandMutCount", GetWTandMutCount, ns = "dpclust3p")

pwalk(
  list(
    data_dpclust$loci_path,
    data_dpclust$alleleFrequencies_path,
    data_dpclust$rho_and_psi_path,
    data_dpclust$subclones_path,
    str_c(outfile.prefix, "/", data_dpclust$Barcode, "/",
          data_dpclust$Barcode, "_allDirichletProcessInfo.txt")
  ),
  ~ runGetDirichletProcessInfo(
    loci_file              = ..1,
    allele_frequencies_file = ..2,
    cellularity_file       = ..3,
    subclone_file          = ..4,
    gender                 = "male", # chrX variants are removed beforehand
    SNP.phase.file         = "NA",
    mut.phase.file         = "NA",
    output_file            = ..5
  )
)

dpclust_single_sample_master_file <-
  data_dpclust %>%
  mutate(
    datafile    = str_c(Barcode, "_allDirichletProcessInfo.txt"),
    cellularity = map_dbl(rho_and_psi, \(x) x$purity_ploidy[[1]][2, 1])
  ) %>%
  select("sample" = Barcode, "subsample" = Barcode, datafile, cellularity)

write_tsv(dpclust_single_sample_master_file,
          str_c(outfile.prefix, "/DPClust_sample_master_file.txt"))

# Write per-sample tool input tables -----------------------------------------

## CloneSig
clonesig_table %>%
  rowwise() %>%
  mutate(path = str_c(local_output_ithinput_dir, "/", Barcode,
                      "/", Barcode, "_clonesig.tsv")) %>%
  ungroup() %>%
  { walk2(.$variant_tables, .$path, ~ write_tsv(.x, .y)) }

## PhylogicNDT (variants + sif)
phylogicndt_table %>%
  rowwise() %>%
  mutate(path = str_c(local_output_ithinput_dir, "/", Barcode,
                      "/", Barcode, "_phylogicndt.tsv")) %>%
  ungroup() %>%
  { walk2(.$variant_tables, .$path, ~ write_tsv(.x, .y)) }

phylogicndt_sif %>%
  rowwise() %>%
  mutate(path = str_c(local_output_ithinput_dir, "/", Barcode,
                      "/", Barcode, "_phylogicndt.sif")) %>%
  ungroup() %>%
  { walk2(.$sif_tables, .$path, ~ write_tsv(.x, .y)) }

## PyClone-VI
pyclonevi_table %>%
  rowwise() %>%
  mutate(path = str_c(local_output_ithinput_dir, "/", Barcode,
                      "/", Barcode, "_pyclonevi.tsv")) %>%
  ungroup() %>%
  { walk2(.$variant_tables, .$path, ~ write_tsv(.x, .y)) }

# Exploratory PCA on shared SNV calls (sanity check) --------------------------

snv_data_acp <-
  snv_data %>%
  mutate(mutation = str_c(Gene, "_", C_point, "_", P_point)) %>%
  select(Barcode, mutation) %>%
  distinct() %>%
  mutate(var = 1) %>%
  pivot_wider(names_from = mutation, values_from = var, values_fill = 0) %>%
  column_to_rownames("Barcode") %>%
  select(-`NA`)

snv_data_acp_pca <- prcomp(snv_data_acp, scale = TRUE)

autoplot(snv_data_acp_pca,
         data   = snv_data_acp,
         colour = "Barcode",
         shape  = "Barcode")
