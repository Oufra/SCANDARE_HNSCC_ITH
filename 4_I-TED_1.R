# =============================================================================
# Script  : 4_I-TED_1.R
# Project : SCANDARE-ORL (HNSCC intra-tumour heterogeneity study)
# Purpose : Compute the intratumour expression distance (I-TED) metric
#           (Martínez-Ruiz et al., Nature 2023) on the multi-region bulk
#           RNA-seq cohort. Standard ITH measures are sensitive to the number
#           of regions sampled per tumour; the pairwise correlation distance
#           used here yields a per-tumour estimate (median of pairwise
#           distances) that is independent of the number of regions sampled.
#
#           For each pair of regions (R_i, R_j) within a patient, the
#           correlation distance is computed with `energy::dcor()` on the
#           variance-stabilised expression of the 500 most variable genes
#           (rows of `assay(vsd)`). The metric is converted to a distance
#           with d = 1 - dcor. Tumour-level I-TED is the median pairwise
#           distance across all available pairs for that patient.
#
# Inputs  : - `vsd`     : DESeq2 vst-transformed SummarizedExperiment
#                        (must be available in the workspace; created
#                        upstream by the bulk RNA-seq differential
#                        expression pipeline).
#           - coldata.xlsx (sheet "Listing 98 cas" range "A4:BK295"):
#               clinical metadata
#           - coldata.xlsx (sheet "deg"): one row per sequenced region,
#               with columns patient, run, region, sequenced_region,
#               hpv_status, site, ITH (pathology class).
# Outputs : - "I.TED metric per patient.csv": one row per patient with
#               I.TEDmedian, I.TEDmin, I.TEDmax.
#           - "I.TED.RData": serialised `vsd` and `coldata`.
# Usage   : Source this script AFTER the upstream DE pipeline has populated
#           `vsd` in the environment.
# Depends : tidyverse, readxl, energy (dcor), matrixStats (rowVars),
#           DESeq2 / SummarizedExperiment (assay), crosstable, rstatix
# =============================================================================

# Libraries -------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(energy)              # dcor()
library(matrixStats)         # rowVars()
library(SummarizedExperiment)# assay()
library(crosstable)
library(rstatix)

# CONFIG ----------------------------------------------------------------------

work_dir         <- "<path/to/tablecounts_raw>"
coldata_xlsx     <- "coldata.xlsx"
ited_csv_out     <- "I.TED metric per patient.csv"
ited_rdata_out   <- "I.TED.RData"

# `vsd` must already exist in the environment (vst-transformed DESeq2
# SummarizedExperiment produced upstream by the bulk RNA-seq DE pipeline).
stopifnot(exists("vsd"))

setwd(work_dir)

# Metadata --------------------------------------------------------------------

# Clinical sheet ("Listing 98 cas"): keep patient id, sequenced run, and the
# clinical/pathological annotations between Age and "Niveau de différenciation".
test <- read_excel(coldata_xlsx,
                   sheet = "Listing 98 cas",
                   range = "A4:BK295")
test <- test %>%
  dplyr::select(1, 34, "Age (an)":"Niveau de différenciation") %>%
  drop_na()

coldata <- read_excel(coldata_xlsx, sheet = "deg") %>%
  drop_na(run) %>%
  left_join(test)

coldata$patient <- str_replace_all(coldata$patient, "-", "_")
coldata$ITH     <- str_replace_all(coldata$ITH, "-", "_")
coldata <- coldata %>% mutate(across(ITH, as_factor))
coldata$ITH <- fct_relevel(coldata$ITH, "I", "II", "III_L", "III_H")
coldata$sequenced_region <- fct_recode(
  coldata$sequenced_region,
  "Poorly_differentiated" = "R",
  "Well_differentiated"   = "V"
)
coldata$sequenced_region <- relevel(coldata$sequenced_region,
                                    ref = "Well_differentiated")

# Baseline characteristics ----------------------------------------------------

crosstable(
  coldata %>% group_by(patient) %>% dplyr::slice(1),
  c(`Age (an)`, Sexe, "Consommation alcool", "Consommation tabac",
    "Localisation tumorale", "Stade tumoral", "Stade ganglionnaire"),
  total = "column", percent_digits = 0
) %>% as_flextable()

crosstable(
  coldata,
  c(sequenced_region, hpv_status, site, region),
  total = "column", percent_digits = 0
) %>% as_flextable()

coldata %>%
  group_by(patient) %>%
  summarise(nb_regions = n()) %>%
  get_summary_stats(nb_regions, type = "common")

# Variance-stabilised expression of the top-500 most variable genes ----------
# (preselection: edgeR::filterByExpr; normalisation: DESeq2 median-of-ratios;
#  transformation: DESeq2 vst — all done upstream).

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
vsd.dcor <- assay(vsd)[topVarGenes, ] %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  gather(run, measurement, D935R17:D1024R48)

# Pairwise correlation distances between regions ------------------------------
#
# `energy::dcor()` does not accept pairwise comparisons, so we build a
# dedicated long->wide data frame for each (R_i, R_j) pair.

## R1-R2
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(1:2) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R1.R2 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R1, R2))

## R1-R3
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(1, 3) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R1.R3 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R1, R3))

## R1-R4
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(1, 4) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R1.R4 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R1, R4))

## R1-R5
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(1, 5) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R1.R5 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R1, R5))

## R2-R3
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(2:3) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R2.R3 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R2, R3))

## R2-R4
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(2, 4) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R2.R4 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R2, R4))

## R2-R5
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(2, 5) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R2.R5 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R2, R5))

## R3-R4
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(3, 4) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R3.R4 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R3, R4))

## R3-R5
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(3, 5) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R3.R5 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R3, R5))

## R4-R5
coldata.dcor <- coldata %>% group_by(patient) %>% dplyr::slice(4, 5) %>% ungroup()
list.1       <- coldata.dcor$patient[duplicated(coldata.dcor$patient)]
coldata.dcor <- coldata.dcor[coldata.dcor$patient %in% list.1, ]
merge.dcor   <- coldata.dcor %>% dplyr::select(1, 2, 6) %>% left_join(vsd.dcor)
R4.R5 <- merge.dcor %>% dplyr::select(-1) %>% spread(region, measurement) %>%
  group_by(patient) %>% summarise(dcor = dcor(R4, R5))

# Combine all pairwise distances --------------------------------------------

df.list <- list(R1.R2, R1.R3, R1.R4, R1.R5,
                R2.R3, R2.R4, R2.R5,
                R3.R4, R3.R5, R4.R5) %>%
  purrr::reduce(full_join, by = "patient")

colnames(df.list)[2:11] <- c("R1.R2", "R1.R3", "R1.R4", "R1.R5",
                             "R2.R3", "R2.R4", "R2.R5",
                             "R3.R4", "R3.R5", "R4.R5")

# Convert dcor to distance: d = 1 - dcor
fct      <- function(x) 1 - x
distance <- df.list %>% mutate(across(where(is.numeric), ~ fct(.)))

distance_long <- distance %>% gather(region, distance, `R1.R2`:`R4.R5`)

# Tumour-level summary: median + min + max of pairwise distances per patient.
I.TED <- distance_long %>%
  group_by(patient) %>%
  drop_na(distance) %>%
  dplyr::summarise(
    I.TEDmax = max(distance),
    I.TEDmin = min(distance)
  )

distance_long <- distance_long %>% left_join(I.TED)

# Per-patient I-TED table that gets written out (median + min + max).
I.TED.patients <- distance_long %>%
  group_by(patient) %>%
  drop_na(distance) %>%
  dplyr::summarise(
    I.TEDmedian = median(distance),
    I.TEDmin    = min(distance),
    I.TEDmax    = max(distance)
  ) %>%
  dplyr::arrange(desc(I.TEDmedian))

# Cohort-level position and dispersion of the median I-TED.
I.TED.patients %>% get_summary_stats(I.TEDmedian)

# Fig. 1b style plot (cf. Martínez-Ruiz et al., Nature 2023).
ggplot(distance_long) +
  geom_point(aes(x = fct_reorder(patient, I.TEDmax, .desc = TRUE),
                 y = distance),
             alpha = 0.5, shape = 21, colour = "darkred", size = 2) +
  geom_linerange(aes(x = patient, ymax = I.TEDmax, ymin = I.TEDmin),
                 color = "darkred", alpha = 0.5) +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  ) +
  ylab("I-TED")

# Export ----------------------------------------------------------------------

write_csv(I.TED.patients, file = ited_csv_out)
save(vsd, coldata, file = ited_rdata_out)
