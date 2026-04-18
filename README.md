# analysis

S4 container class `Analysis` and all downstream analysis methods for bulk multi-omic cancer cohort analysis. This is the top-level package in the HaifengPackages system.

## Role in the Pipeline

```
utiltools
    └── assay
         └── bulkexperiment
              └── analysis  <-- YOU ARE HERE
```

`analysis` wraps a fully-processed `BulkExperiment` in an `Analysis` S4 object and provides all `run*Analysis()` methods for downstream statistical analysis and visualization: survival, differential analysis, consensus clustering, mutational signatures, CNV, fusion, clonal evolution, prognostic signature building, and more.

## The Analysis Class

```
Analysis
├── project / analysis_name   — output path components
├── experiment                — BulkExperiment (the processed cohort)
└── output                    — named list of results (matrices, data.frames, plot objects)
```

All `run*Analysis()` methods:
1. Store results in `obj@output[["result_name"]]` as R objects (data.frames, matrices, heatmap objects)
2. Write CSV + SVG files to `{project}/{analysis_name}/` when `save_analysis = TRUE`

Three internal helper functions back multiple public methods: `deltaAnalysis()` (longitudinal delta), `difAnalysis()` (differential), `signatureAnalysis()` (group signatures).

## Analysis Methods

### Survival and Endpoint Analyses

| Method | Description | Key params |
|---|---|---|
| `runSurvivalAnalysis()` | Univariate + multivariate Cox, KM plots | `survival_column`, `status_column`, `covariate_columns`, `extra_covariates` |
| `runCINTMBSurvivalAnalysis()` | KM survival by CIN-High/Low and TMB-High/Low | CIN/TMB column names, cutoffs (median used if omitted) |
| `runAssaySurvivalSignatureAnalysis()` | LASSO-Cox prognostic signature, risk score, time-ROC AUC | `assay_name`, `survival_column`, `uni_cox_pvalue_cutoff`, `beta_filter` |

### Variant and Genomic Analyses

| Method | Description | Key params |
|---|---|---|
| `runVariantStatisticsAnalysis()` | Per-sample CIN, TMB, COSMIC-TMB, AMP/DEL gene counts | `cosmic_tier12_cancer_genes`, `genetic_cancer_genes`, `oncokb_cancer_genes` |
| `runCNVAnalysis()` | GISTIC2 lesion heatmap, stable CNV clusters, KM; genome-wide and arm CN plots | `top_anno`, `clusters`, `survival_column` |
| `runCNVDifAnalysis()` | Differential CNV by group at cytoband/arm/gene level | `group_column`, `level`, `contrasts` |
| `runMutationalSignatureAnalysis()` | NMF signature extraction, COSMIC v3 matching, per-sample contributions | `profile_group`, `rank` (review rank plot first) |
| `runMutationalMafAnalysis()` | maftools oncoplot with cancer hallmarks, lollipop plots, somatic interactions | `cancer_genes`, `selected_genes`, `pathways` |
| `runMDLAnalysis()` | MDL panel + WES harmonized oncoplot (timepoint-aware) | `mdl_info`, `mutation_panel_info`, `oncokb_cancer_genes_info` |
| `runFusionAnalysis()` | Multi-panel fusion heatmap (type, caller, gene function) | `top_anno`, `top_fusions` |

### Transcriptomic Analyses

| Method | Description | Key params |
|---|---|---|
| `runAssayDifAnalysis()` | Differential analysis for any assay (DGE, immune dif, activity dif, etc.) | `assay_name`, `group_column`, `test_method`, `block_column` (paired), `contrasts` |
| `runAssayDeltaAnalysis()` | Longitudinal delta (TP2 minus Baseline) for any assay + group comparison on delta | `assay_name`, `timepoint_column`, `timepoint_terms`, `group_column` |
| `runAssayunSupervisedAnalysis()` | Stable consensus clustering (13 algorithms) for any assay | `assay_name`, `scale`, `clusters`, `method` |
| `runAssaySignatureAnalysis()` | Per-group gene/feature signatures with optional GO enrichment | `assay_name`, `group_column`, `sign_filter`, `ora` |
| `runICGDifAnalysis()` | Immune checkpoint gene inhibitory/stimulatory scoring + differential + survival | `icg_info`, `survival_column`, `group_column` |

### Clonal Evolution Analyses

| Method | Description | Notes |
|---|---|---|
| `runCloneAnalysis()` | Load PyClone-VI CCF results into `clone_assay` | Requires prior PyClone-VI run (external Python) |
| `runCCFDifAnalysis()` | Differential CCF between groups | Uses `clone_assay` |
| `runCloneunSupervisedAnalysis()` | Unsupervised clustering of clone fractions + KM | Uses `clone_assay` |
| `runpyNBSAnalysis()` | pyNBS cluster survival visualization | Requires prior pyNBS run (external Python) |

## Installation

```r
# Install all dependencies in order
devtools::install("utiltools")
devtools::install("assay")
devtools::install("bulkexperiment")
devtools::install("analysis")

# For interactive development:
devtools::load_all("analysis")
```

## Basic Usage

```r
library(Analysis)

# Wrap a processed BulkExperiment in an Analysis object
obj <- new("Analysis",
    project       = "MyProject",
    analysis_name = "cohort_analysis_v1",
    experiment    = exp_sub              # a processed, subsetted BulkExperiment
)

# Survival analysis
obj <- runSurvivalAnalysis(obj,
    survival_column   = "OS_days",
    status_column     = "VitalStatus",
    covariate_columns = c("Age", "Sex", "TreatmentArm"),
    save_analysis     = TRUE
)

# Differential mRNA by treatment arm
obj <- runAssayDifAnalysis(obj,
    assay_name    = "mRNA_assay",
    group_column  = "TreatmentArm",
    test_method   = "limma",
    contrasts     = list(c("ArmA", "ArmB")),
    save_analysis = TRUE
)

# Consensus clustering of immune cell fractions
obj <- runAssayunSupervisedAnalysis(obj,
    assay_name    = "icf_assay",
    clusters      = 3,
    method        = "combined",
    save_analysis = TRUE
)

# Longitudinal delta analysis (paired baseline/TP2)
obj <- runAssayDeltaAnalysis(obj,
    assay_name       = "mRNA_assay",
    timepoint_column = "Timepoint",
    timepoint_terms  = list(baseline = "Baseline", tp2 = "TP2"),
    group_column     = "TreatmentArm",
    save_analysis    = TRUE
)

# Mutational signatures (NMF + COSMIC v3)
# First run without specifying rank to review the rank plot
obj <- runMutationalSignatureAnalysis(obj, rank = NULL)
# REVIEW: NMF rank estimation plot (cophenetic, RSS)
obj <- runMutationalSignatureAnalysis(obj, rank = 3, save_analysis = TRUE)

# LASSO-Cox prognostic signature from mRNA
obj <- runAssaySurvivalSignatureAnalysis(obj,
    assay_name          = "mRNA_assay",
    survival_column     = "OS_days",
    status_column       = "VitalStatus",
    uni_cox_pvalue_cutoff = 0.01,
    save_analysis       = TRUE
)

# Access results
obj@output[["survival"]]
obj@output[["mRNA_assay_differential"]]
```

## Adding a New Analysis Method

1. Add `setGeneric` + `setMethod` to `analysis/R/Analysis.R` following the existing pattern.
2. Put heavy computation in a standalone helper function (not an S4 method) — see `deltaAnalysis()`, `difAnalysis()`, `signatureAnalysis()` as examples.
3. Store results in `obj@output[["result_name"]]` and write to disk when `save_analysis = TRUE`.
4. Run `devtools::document("analysis")` to regenerate NAMESPACE.

## Dependencies

Internal: `bulkexperiment`, `assay`, `utiltools`

External: `ComplexHeatmap`, `survivalAnalysis`, `glmnet`, `maftools`, `MutationalPatterns`, `NMF`, `cosmicsig`, `ClassDiscovery`, `BSgenome.Hsapiens.UCSC.hg19`, `GenomeInfoDb`, `clonevol`, `deconstructSigs`, `dplyr`, `tidyr`, `purrr`, `RColorBrewer`, `grid`, `assertthat`

## Author

Haifeng Zhu, MD Anderson Cancer Center (`zhuhf77@mdanderson.org`)
