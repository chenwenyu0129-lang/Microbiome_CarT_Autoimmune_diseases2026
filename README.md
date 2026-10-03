# Microbiome_CarT_Autoimmune_diseases2026
# Analysis Workflows for Longitudinal Microbiome and Single-cell Data

## Overview

This repository contains selected R and Python workflows for longitudinal microbiome and single-cell analyses in a clinical study. The scripts cover microbial diversity, community composition, outcome-associated bacterial features, host–microbiome relationships, metagenomic functions, and B-cell transcriptional profiles.

The workflows are organized around analysis questions and figure outputs rather than as a single end-to-end pipeline. Each R Markdown script documents its settings, reads its required inputs, and exports analysis-specific results. The single-cell notebook contains downstream analyses of a processed, annotated AnnData object, including cell-type visualization, marker expression, trajectory analysis, cell-composition PCA, and functional enrichment.

## Repository organization

```text
Code/
  01.alpha_diversity_pairwise_heatmap.Rmd
  02.16sPCoA.Rmd
  03.16sPCoA.PERMANOVA.Rmd
  04.16s_BrayCurtis_and_metadata_distributions.Rmd
  05.lefseCompareWithBaseline.Rmd
  06.lefse_Bsub_clone.Rmd
  07.outcomeM1_AUC.Rmd
  008.bacteriaScoreClinicalCorrelation.Rmd
  009.pathway.Rmd
  010.scRNA-seq.ipynb
  corresponding microbiome HTML reports
Demo/
  scRNAseq_demo.py
  Demodata/
README.md
```

The `.Rmd` files contain the microbiome analysis code; the corresponding `.html` files provide rendered reports. Input files and generated output directories are separate from the analysis source files.

## Analysis workflows

| Script | Analysis focus |
|---|---|
| [01.alpha_diversity_pairwise_heatmap.Rmd](Code/01.alpha_diversity_pairwise_heatmap.Rmd) | Longitudinal alpha diversity and pairwise comparisons. |
| [02.16sPCoA.Rmd](Code/02.16sPCoA.Rmd) | Microbial community ordination using CLR-PCA. |
| [03.16sPCoA.PERMANOVA.Rmd](Code/03.16sPCoA.PERMANOVA.Rmd) | Bray–Curtis PERMANOVA and sensitivity analyses. |
| [04.16s_BrayCurtis_and_metadata_distributions.Rmd](Code/04.16s_BrayCurtis_and_metadata_distributions.Rmd) | Bray–Curtis distances and sample-characteristic distributions. |
| [05.lefseCompareWithBaseline.Rmd](Code/05.lefseCompareWithBaseline.Rmd) | LEfSe comparisons with baseline. |
| [06.lefse_Bsub_clone.Rmd](Code/06.lefse_Bsub_clone.Rmd) | B-cell reconstitution-associated bacterial features. |
| [07.outcomeM1_AUC.Rmd](Code/07.outcomeM1_AUC.Rmd) | Bacterial score models for M1 outcome. |
| [008.bacteriaScoreClinicalCorrelation.Rmd](Code/008.bacteriaScoreClinicalCorrelation.Rmd) | Bacterial scores and clinical-indicator correlations. |
| [009.pathway.Rmd](Code/009.pathway.Rmd) | Metagenomic pathway and KO functional analyses. |
| [010.scRNA-seq.ipynb](Code/010.scRNA-seq.ipynb) | Single-cell visualization, trajectory, and functional analyses. |

## Software requirements

### R workflows

Use R with R Markdown and Pandoc support, for example through RStudio. Required packages include:

- CRAN: `rmarkdown`, `knitr`, `dplyr`, `tidyr`, `tibble`, `purrr`, `ggplot2`, `ggpubr`, `scales`, `stringr`, `forcats`, `openxlsx`, `vegan`, `lme4`, `lmerTest`, `emmeans`, `factoextra`, `rstatix`, and `circlize`.
- Bioconductor: `lefser`, `SummarizedExperiment`, `ComplexHeatmap`, `limma`, `clusterProfiler`, and `enrichplot`.

Install the packages in R:

```r
install.packages(c(
  "rmarkdown", "knitr", "dplyr", "tidyr", "tibble", "purrr",
  "ggplot2", "ggpubr", "scales", "stringr", "forcats", "openxlsx",
  "vegan", "lme4", "lmerTest", "emmeans", "factoextra", "rstatix", "circlize"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "lefser", "SummarizedExperiment", "ComplexHeatmap",
  "limma", "clusterProfiler", "enrichplot"
))
```

### Python workflows

The notebook uses `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`, `gseapy`, and `adjustText`, together with Jupyter. The trajectory section additionally requires the `py_monocle` module exposing `learn_graph` and `order_cells`; provide that module in the Python environment before running the section.

```bash
pip install scanpy anndata pandas numpy matplotlib seaborn scikit-learn gseapy adjustText jupyterlab
```

To execute the additional example steps present in the currently included demo, install Harmony, Leiden, and Louvain dependencies:

```bash
pip install harmonypy igraph leidenalg louvain
```

Use mutually compatible package versions and record the environment used for each analysis. Runtime and memory requirements depend on input size and hardware.

## Running the analyses

1. Download the repository and install the required R or Python packages.
2. Place the required input files at the paths specified in the corresponding scripts.
3. Run the selected workflow and review its generated figures, tables, and model outputs.

For microbiome analyses, open `Code/` in RStudio and knit the desired `.Rmd` file to HTML, or render it from R. For example, from the repository root:

```r
setwd("Code")
rmarkdown::render("01.alpha_diversity_pairwise_heatmap.Rmd")
```

Run `07.outcomeM1_AUC.Rmd` before `008.bacteriaScoreClinicalCorrelation.Rmd`, which uses its bacterial score model.

For single-cell analyses, open `Code/010.scRNA-seq.ipynb` in Jupyter, load the processed and annotated AnnData object as `adata_B`, and run the relevant cells. The script in `Demo/` illustrates the data-merging procedure. To run the demo, check its configured input paths and run from the repository root:

```bash
python Demo/scRNAseq_demo.py
```

The demo is for method illustration. To reproduce Harmony integration and the full single-cell analyses, request the original data and necessary analysis materials from the study team.

## Reproducibility and data availability

- Analysis settings, thresholds, and random seeds are specified in the corresponding scripts where applicable.
- Retain consistent identifiers across abundance matrices, metadata, clinical measurements, and cell metrics.
- Required analysis inputs must be supplied separately when they are not included in the repository. Rendered HTML reports are not substitutes for input data.
- Original data and materials needed for the full single-cell analysis should be requested from the study team. The demonstration inputs and outputs do not replace those materials.
- Statistical models and sensitivity analyses should be interpreted in the context of the study design and available observations.
- Record software versions, input preprocessing, and annotation procedures when reproducing results.

The repository provides selected analysis workflows and reports; it does not package every upstream processing step or every required input.
