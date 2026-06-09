# =============================================================================
# Script  : 02_1_dpclust_pipeline.R
# Project : SCANDARE-ORL (HNSCC intra-tumour heterogeneity study)
# Purpose : Run the DPClust clustering pipeline on a single sample (one region).
#           This is the DPClust v2.2.6 wrapper distributed by the Wellcome
#           Sanger Institute (sd11 [at] sanger.ac.uk, 2018-11-01),
#           reproduced verbatim except for surrounding documentation and
#           cosmetic formatting. The clustering logic is unchanged.
# Inputs  : - Master TSV (`-i`) with columns: sample, subsample, datafile,
#             cellularity (produced by 01_preprocess_ITH_4tools.R).
#           - `-d` directory containing the `_allDirichletProcessInfo.txt`
#             datafiles referenced by the master TSV.
# Outputs : One sub-directory per sample under `-o`, containing the
#           DPClust clustering results (bestClusterInfo, bestConsensusAssignments,
#           MCMC diagnostics, etc.).
# Usage   : Rscript 02_1_dpclust_pipeline.R \
#             -r <SAMPLE_INDEX> \
#             -d <path/to/dpclust_input_dir> \
#             -o <path/to/dpclust_output_dir> \
#             -i <path/to/DPClust_sample_master_file.txt>
#           (one process per sample; can be parallelised with a job array).
# =============================================================================

# Command line options --------------------------------------------------------

library(optparse)

option_list <- list(
  make_option(c("-r", "--run_sample"),    type = "integer",   default = NULL,
              help = "Sample to run",                metavar = "character"),
  make_option(c("-d", "--data_path"),     type = "character", default = NULL,
              help = "Path to where dpinput data files are stored",
              metavar = "character"),
  make_option(c("-o", "--outputdir"),     type = "character", default = NULL,
              help = "Directory where the output is saved",
              metavar = "character"),
  make_option(c("-i", "--input"),         type = "character", default = NULL,
              help = "Datafile with sample information",
              metavar = "character"),
  make_option(c("-k", "--keep_temp_files"), type = "logical", action = "store_true",
              default = FALSE,
              help = "Keep intermediate output files",
              metavar = "character"),
  make_option(c("-a", "--analysis_type"), type = "character", default = "nd_dp",
              help = "Analysis type to run", metavar = "character"),
  make_option(c("--iterations"),          type = "integer",   default = 2000,
              help = "Number of iterations to run the MCMC chain",
              metavar = "character"),
  make_option(c("--burnin"),              type = "integer",   default = 1000,
              help = "Number of iterations to discard as burnin",
              metavar = "character"),
  make_option(c("--mut_assignment_type"), type = "integer",   default = 1,
              help = "Mutation assignment method", metavar = "character"),
  make_option(c("--min_muts_cluster"),    type = "integer",   default = -1,
              help = paste("Minimum number of mutations per cluster required",
                           "for it to be kept in the final output, set to -1",
                           "to disable (default), see also --min_frac_muts_cluster"),
              metavar = "character"),
  make_option(c("--min_frac_muts_cluster"), type = "numeric", default = 0.01,
              help = paste("Minimum fraction of mutations per cluster required",
                           "for it to be kept in the final output, set to -1",
                           "to disable, see also --min_muts_cluster"),
              metavar = "character"),
  make_option(c("--num_muts_sample"),     type = "integer",   default = 50000,
              help = "Number of mutations from which downsampling starts",
              metavar = "character"),
  make_option(c("--bin_size"),            type = "double",    default = NULL,
              help = paste("Binsize to use when constructing multi-dimensional",
                           "density - only used when number of samples > 1"),
              metavar = "character"),
  make_option(c("--seed"),                type = "integer",   default = 123,
              help = "Provide a seed", metavar = "character"),
  make_option(c("--assign_sampled_muts"), type = "integer",   default = TRUE,
              help = "Whether to assign mutations that have been removed during sampling",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

run                   <- opt$run_sample
datpath               <- opt$data_path
outdir                <- opt$outputdir
purity_file           <- opt$input
analysis_type         <- opt$analysis_type
no.iters              <- opt$iterations
no.iters.burn.in      <- opt$burnin
mut.assignment.type   <- opt$mut_assignment_type
num_muts_sample       <- opt$num_muts_sample
bin.size              <- opt$bin_size
seed                  <- opt$seed
assign_sampled_muts   <- opt$assign_sampled_muts
keep_temp_files       <- opt$keep_temp_files
min_muts_cluster      <- opt$min_muts_cluster
min_frac_muts_cluster <- opt$min_frac_muts_cluster

if (is.null(outdir)) { outdir <- getwd() }

# Fixed parameters ------------------------------------------------------------

is.male                  <- TRUE
sample.snvs.only         <- TRUE   # downsample on SNVs only, not on CNAs
remove.snvs              <- FALSE  # cluster on CNAs only
generate_cluster_ordering <- FALSE
species                  <- "human"

# Co-cluster CNA parameters
co_cluster_cna             <- FALSE
add.conflicts              <- FALSE  # conflicts matrix (CN + mut2mut phasing)
cna.conflicting.events.only <- FALSE
num.clonal.events.to.add   <- 1
min.cna.size               <- 100    # 10 kb units

# Avoid RGL X11 error on headless nodes
options(rgl.useNULL = TRUE)
suppressMessages(library(DPClust))

# Process input ---------------------------------------------------------------

sample2purity <- read.table(purity_file, header = TRUE, stringsAsFactors = FALSE)
samplename    <- unique(sample2purity$sample)[run]
datafiles     <- sample2purity[sample2purity$sample == samplename, ]$datafile
subsamples    <- sample2purity[sample2purity$sample == samplename, ]$subsample
cellularity   <- sample2purity[sample2purity$sample == samplename, ]$cellularity

if ("sex" %in% colnames(sample2purity)) {
  is.male     <- (sample2purity[sample2purity$sample == samplename, ]$sex == "male")[1]
  cndatafiles <- sample2purity[sample2purity$sample == samplename, ]$cndatafile
} else {
  is.male     <- TRUE
  cndatafiles <- NA
}

if ("mutphasing" %in% colnames(sample2purity)) {
  mutphasingfiles <- sample2purity[sample2purity$sample == samplename, ]$mutphasing
} else {
  mutphasingfiles <- NA
}

# Status message --------------------------------------------------------------

print("")
print(paste("Running:",        samplename,     sep = " "))
print(paste("Working dir:",    outdir,         sep = " "))
print(paste("Analysis type:",  analysis_type,  sep = " "))
print("Datafiles:")
print(datafiles)
print("")

# Output path -----------------------------------------------------------------

outdir <- file.path(outdir, paste(samplename, "_DPoutput_", no.iters,
                                  "iters_", no.iters.burn.in,
                                  "burnin_seed", seed, "/", sep = ""))

# Setup parameters ------------------------------------------------------------

run_params <- make_run_params(
  no.iters, no.iters.burn.in, mut.assignment.type, num_muts_sample,
  is.male                   = is.male,
  min_muts_cluster          = min_muts_cluster,
  min_frac_muts_cluster     = min_frac_muts_cluster,
  species                   = species,
  assign_sampled_muts       = assign_sampled_muts,
  keep_temp_files           = keep_temp_files,
  generate_cluster_ordering = generate_cluster_ordering
)
sample_params   <- make_sample_params(datafiles, cellularity, is.male,
                                      samplename, subsamples, mutphasingfiles)
advanced_params <- make_advanced_params(seed)

# Run clustering --------------------------------------------------------------

RunDP(
  analysis_type   = analysis_type,
  run_params      = run_params,
  sample_params   = sample_params,
  advanced_params = advanced_params,
  outdir          = outdir,
  cna_params      = NULL,
  mutphasingfiles = NULL
)
