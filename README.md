# Integrative transcriptomic and network analysis across inherited muscle disorders

This repository contains the analysis code supporting *Integrative transcriptomic and network analysis identifies AP-1 as a shared regulatory axis across inherited muscle disorders*.

The study integrates publicly available transcriptomic data from Duchenne muscular dystrophy (DMD), Pompe disease, Emery-Dreifuss muscular dystrophy (EDMD), and DNM2-related centronuclear myopathy (CNM). It identifies convergent expression and pathway-level signatures, infers transcription-factor activity, and reconstructs an AP-1/JUN-centred regulatory network.

## Repository contents

- `integrative_transcriptomic_network_analysis.Rmd` - main analysis and figure generation.
- `function.R` - helper functions used by the report.
- `data/` - processed inputs used directly by the main analysis.
- `DEGs and preprocess/` - disease-specific preprocessing scripts for the four source datasets.

## Running the main analysis

Open `integrative_transcriptomic_network_analysis.Rproj` in RStudio, install the packages loaded at the start of the R Markdown file, then render:

```r
rmarkdown::render("integrative_transcriptomic_network_analysis.Rmd")
```

The report reads the processed `.RDS` files in `data/` and writes generated figures to `figures/`.

## Reproducing the processed inputs

The preprocessing scripts correspond to DMD (GSE38417), EDMD (GSE204804), DNM2-related CNM (GSE160078), and Pompe disease (GSE38680). To prevent accidental replacement of the included processed inputs, each script uses `write_processed_data <- FALSE` by default. Set it to `TRUE` only when you intend to regenerate the corresponding `.RDS` file.

The DNM2 preprocessing script additionally requires the raw count matrix `data/GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx`, which is not included in this repository.

## Notes

The main analysis requires internet access on first use for STRING interactions and KEGG pathway resources used by `pathview`. These resources may subsequently be cached by their respective R packages.
