# =============================================================================
# Script  : 02_2_run_clonesig.py
# Project : SCANDARE-ORL (HNSCC intra-tumour heterogeneity study)
# Purpose : Run CloneSig (Abecassis et al., Nat Commun 2021) on one
#           pre-processed per-sample variant table.
# Inputs  : - Pre-processed CloneSig TSV produced by
#             01_preprocess_ITH_4tools.R, with columns:
#               alteration, pattern, variant_allele_count, total_read_count,
#               freq, C_tumor_tot, C_tumor_minor, purity
# Outputs : - {sample}_clonesig_out.tsv   : per-SNV cluster assignment,
#             estimated CCF, signature, multiplicity
#           - {sample}_clonesig_out.pickle: pickled CloneSig estimator object
# Method  : CloneSig is run with the breast-cancer-specific reference
#           signature set (cancer_type = 15) which restricts COSMIC v2
#           signatures to those reported by Alexandrov et al. (2013) in
#           breast tumours. NOTE: for an HNSCC study cancer_type would
#           normally be 6 ("Head and Neck"); the value 15 (breast) was used
#           to match the configuration applied in the manuscript.
# Usage   : python 02_2_run_clonesig.py <input.tsv> <output_dir/>
# Depends : clonesig, numpy, pandas, scipy
# =============================================================================

import os
import re
import sys
import pickle

import numpy as np
import pandas as pd
from scipy import special, spatial, stats

from clonesig.run_clonesig import get_MU, run_clonesig
from clonesig.data_loader import get_context, PAT_LIST


# CLI -------------------------------------------------------------------------

inputfile = sys.argv[1]
outputdir = sys.argv[2]

print("Input file is " + inputfile)

# Load and clean input --------------------------------------------------------

clonesig_table = pd.read_csv(inputfile, sep="\t")
clonesig_table["normal_cn"] = 2

# Extract sample name from file name (expects "<sample>_clonesig.tsv").
filename = os.path.basename(inputfile)
sample_name = re.findall("(.+)_clonesig", filename)[0]

# Drop homozygous deletions (cause errors). Redundant with the preprocessing
# script but kept here as a safeguard.
clonesig_table = clonesig_table[clonesig_table["C_tumor_tot"] > 0]

# Encode the trinucleotide context as the index in CloneSig's PAT_LIST.
clonesig_table = clonesig_table.assign(
    trinucleotide=clonesig_table.apply(
        lambda x: PAT_LIST.index(x["pattern"]), axis=1
    )
)

# Tumour purity (a single value per sample).
purity = float(np.unique(clonesig_table.purity))

# Reference signature set (see Method note in the header).
default_MU = get_MU()

# Run CloneSig ----------------------------------------------------------------

est, lr, pval, new_inputMU, cst_est, future_sigs = run_clonesig(
    np.array(clonesig_table.trinucleotide),
    np.array(clonesig_table.variant_allele_count),
    np.array(clonesig_table.total_read_count),
    np.array(clonesig_table.normal_cn),
    np.array(clonesig_table.C_tumor_tot),
    np.array(clonesig_table.C_tumor_minor),
    float(purity),
    default_MU,
)

# Assemble per-SNV result table -----------------------------------------------

est_table = pd.DataFrame({
    "alteration":     clonesig_table.alteration,
    "trinucleotide":  est.T,
    "var_counts":     est.B,
    "minor_cn":       est.C_tumor_minor,
    "major_cn":       est.C_tumor_major,
    "total_cn":       est.C_tumor_tot,
    "depth":          est.D,
    "clone":          est.qun.argmax(axis=1),
    "snv_proportion": est.xi[est.qun.argmax(axis=1)],
    "clone_ccf":      est.phi[est.qun.argmax(axis=1)],
    "signature":      np.arange(default_MU.shape[0])[
        est.rnus[np.arange(est.N), est.qun.argmax(axis=1), :].argmax(axis=1)
    ],
    "mult":           est.vmnu[
        np.arange(est.N), est.qun.argmax(axis=1), :
    ].argmax(axis=1) + 1,
})

# Write outputs ---------------------------------------------------------------

est_table.to_csv(outputdir + sample_name + "_clonesig_out.tsv", sep="\t")

with open(outputdir + sample_name + "_clonesig_out.pickle", "wb") as f:
    pickle.dump(est, f)
