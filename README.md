# SCANDARE_HNSCC_ITH

Code accompanying the study of intra-tumour heterogeneity (ITH) in the
**SCANDARE** prospective cohort of patients with resectable head-and-neck
squamous-cell carcinoma (HNSCC) treated by upfront surgery. The cohort comprises
87 patients with 2–5 spatially distinct tumour regions per patient, profiled by
targeted DNA sequencing and bulk RNA sequencing.

The scripts cover:

- **Genomic ITH (gITH)** — clonal deconvolution with four independent tools
  (CloneSig [1], PyClone-VI [2], PhylogicNDT [3], DPClust [4]), consensus
  reconciliation with WeMe [5], and per-sample Shannon-index scoring.
- **Transcriptional ITH (tITH)** — pairwise distance-correlation between
  regions of the same tumour (I-TED metric [6]).

## Repository layout

```
scripts/
├── 00_CNV_QC.Rmd                     # FACETS CNV QC and NRPCC computation
├── 01_preprocess_ITH_4tools.R        # Build tool-specific input tables from
│                                     #   DRAGON SNV + FACETS CNV outputs
├── 02_1_dpclust_pipeline.R           # DPClust clustering (Sanger wrapper)
├── 02_2_run_clonesig.py              # CloneSig
├── 02_3_run_pyclone-vi.sh            # PyClone-VI
├── 02_4_run_phylogicndt.sh           # PhylogicNDT (qsub array)
├── 03_0_preprocess_results_for_WeMe.R# Harmonize the 4 outputs for WeMe
├── 03_1_run_weme.R                   # WeMe consensus
├── 4_I-TED_1.R                       # Transcriptional ITH (I-TED)
└── ITH_DNA_analysis.qmd              # Per-sample Shannon ITH scoring +
                                      #   inter-tool comparison
```

## Pipeline overview

```
                       ┌──────────────────────────────────────────────┐
                       │ 00_CNV_QC.Rmd                                │
                       │  FACETS CNV QC + NRPCC                       │
                       └───────────────────┬──────────────────────────┘
                                           │
                       ┌───────────────────▼──────────────────────────┐
                       │ 01_preprocess_ITH_4tools.R                   │
                       │  Build per-tool input tables                 │
                       └───────────────────┬──────────────────────────┘
                                           │
        ┌──────────────────┬───────────────┴────────────────┬─────────────────────┐
        ▼                  ▼                                ▼                     ▼
┌───────────────┐  ┌───────────────┐               ┌────────────────┐   ┌──────────────────┐
│ 02_1 DPClust  │  │ 02_2 CloneSig │               │ 02_3 PyClone-VI│   │ 02_4 PhylogicNDT │
└───────┬───────┘  └───────┬───────┘               └────────┬───────┘   └────────┬─────────┘
        └──────────────────┴────────────────┬───────────────┴────────────────────┘
                                            │
                       ┌────────────────────▼─────────────────────────┐
                       │ 03_0_preprocess_results_for_WeMe.R           │
                       │  Reshape per-tool outputs for WeMe           │
                       └────────────────────┬─────────────────────────┘
                                            │
                       ┌────────────────────▼─────────────────────────┐
                       │ 03_1_run_weme.R    Consensus subclonal       │
                       │                    structure                 │
                       └────────────────────┬─────────────────────────┘
                                            │
                       ┌────────────────────▼─────────────────────────┐
                       │ ITH_DNA_analysis.qmd                         │
                       │  Per-sample Shannon ITH score + cross-tool   │
                       │  comparison                                  │
                       └──────────────────────────────────────────────┘

   ┌─────────────────────────────────┐
   │ 4_I-TED_1.R   (independent)     │
   │  Transcriptional ITH (I-TED)    │
   │  from bulk RNA-seq vsd matrix   │
   └─────────────────────────────────┘
```

## Inputs expected

The pipeline starts from outputs of the in-house **DRAGON** somatic-variant
calling and copy-number pipeline (one run per region):

- `variant_report.table1.tsv` — SNV/indel calls with strand-specific read
  counts.
- `cnv_report.table2.tsv`     — FACETS segments with `BAF` encoding, ploidy,
  cellularity.
- `mapping.stats.xlsx`        — mapping QC (used to derive NRPCC).

A sample list (TSV) with one row per sequenced region is also required, with at
least the columns `ID_SCANDARE` and `Barcode_DRAGON`. The transcriptional ITH
script (`4_I-TED_1.R`) expects a DESeq2 `vst`-transformed `SummarizedExperiment`
(`vsd`) and a `coldata.xlsx` clinical-metadata file produced by the bulk
RNA-seq pipeline.

## Configuring local paths

Every script exposes a `CONFIG` block at the top with placeholders of the form
`<path/to/...>`. **Edit these before running.** No path needs to be modified
elsewhere in the scripts.

## Dependencies

### R (≥ 4.2)

- CRAN: `tidyverse`, `vroom`, `janitor`, `readxl`, `ggpubr`, `ggfortify`,
  `corrplot`, `optparse`, `rjson`, `furrr`, `energy`, `matrixStats`,
  `crosstable`, `rstatix`, `kableExtra`, `DT`.
- Bioconductor: `GenomicRanges`, `Rsamtools`, `plyranges`,
  `VariantAnnotation`, `SummarizedExperiment`, `DESeq2`.
- GitHub: [`dpclust3p`](https://github.com/Wedge-lab/dpclust3p),
  [`DPClust`](https://github.com/Wedge-lab/dpclust),
  [`weme`](https://github.com/morrislab/weme).

### Python (≥ 3.8)

- [`clonesig`](https://github.com/judithabk6/clonesig)
- [`pyclone-vi`](https://github.com/Roth-Lab/pyclone-vi)

### Other

- [PhylogicNDT](https://github.com/broadinstitute/PhylogicNDT) (run from a
  Singularity image; tested with PhylogicNDT v1.0).
- A PBS/Torque scheduler (only `02_4_run_phylogicndt.sh`).
- hg19 reference FASTA, indexed with `samtools faidx`.

External resources used by `01_preprocess_ITH_4tools.R`:

- [OncoKB cancer gene list](https://www.oncokb.org/cancer-genes) (TSV download).

## Running

The scripts are numbered in the order they should be executed:

```bash
# 0. QC
Rscript -e 'rmarkdown::render("scripts/00_CNV_QC.Rmd")'

# 1. Build per-tool input tables
Rscript scripts/01_preprocess_ITH_4tools.R

# 2. Run the four deconvolution tools (parallelisable across samples)
for i in $(seq 1 N_SAMPLES); do
  Rscript scripts/02_1_dpclust_pipeline.R \
    -r $i \
    -d <path/to/ith_inputs> \
    -o <path/to/results/dpclust> \
    -i <path/to/ith_inputs/DPClust_sample_master_file.txt>
done

for f in <path/to/ith_inputs>/*/*_clonesig.tsv; do
  python scripts/02_2_run_clonesig.py "$f" <path/to/results/clonesig/>
done

for f in <path/to/ith_inputs>/*/*_pyclonevi.tsv; do
  ./scripts/02_3_run_pyclone-vi.sh -i "$f" -o <path/to/results/pyclone_vi/>
done

./scripts/02_4_run_phylogicndt.sh    # submits one qsub job per sample

# 3. Consensus
Rscript scripts/03_0_preprocess_results_for_WeMe.R
(cd <path/to/weme_input> && Rscript /path/to/scripts/03_1_run_weme.R)

# 4. Final ITH scoring and analysis
quarto render scripts/ITH_DNA_analysis.qmd

# Transcriptional ITH (independent of the DNA pipeline; requires upstream
# bulk RNA-seq pipeline to have produced `vsd`)
Rscript scripts/4_I-TED_1.R
```

## Data availability

Raw sequencing data are subject to controlled access; processed per-sample
subclonal structures, ITH scores, and I-TED values are deposited as
supplementary material to the manuscript.

## Citation

> *Manuscript reference to be added upon acceptance.*

## References

The tools applied by this pipeline are described in the following publications. Please cite them in addition to the SCANDARE manuscript when reusing this code.

1. **CloneSig** — Abécassis, J., Reyal, F. & Vert, J.-P. CloneSig can jointly infer intra-tumor heterogeneity and mutational signature activity in bulk tumor sequencing data. *Nature Communications* **12**, 5352 (2021). <https://doi.org/10.1038/s41467-021-24992-y>

2. **PyClone-VI** — Gillis, S. & Roth, A. PyClone-VI: scalable inference of clonal population structures using whole genome data. *BMC Bioinformatics* **21**, 571 (2020). <https://doi.org/10.1186/s12859-020-03919-2>

3. **PhylogicNDT** — Leshchiner, I. *et al.* Comprehensive analysis of tumour initiation, spatial and temporal progression under multiple lines of treatment. *bioRxiv* 508127 (2019). <https://doi.org/10.1101/508127>

4. **DPClust** — Nik-Zainal, S. *et al.* The life history of 21 breast cancers. *Cell* **149**, 994–1007 (2012). <https://doi.org/10.1016/j.cell.2012.04.023>

5. **WeMe consensus** — Dentro, S. C. *et al.* Characterizing genetic intra-tumor heterogeneity across 2,658 human cancer genomes. *Cell* **184**, 2239–2254.e39 (2021). <https://doi.org/10.1016/j.cell.2021.03.009>

6. **I-TED metric** — Martínez-Ruiz, C. *et al.* Genomic–transcriptomic evolution in lung cancer and metastasis. *Nature* **616**, 543–552 (2023). <https://doi.org/10.1038/s41586-023-05706-4>

## Contact

Abderaouf Hamza — Institut Curie, Paris.
