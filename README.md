# 🧬 CanisLupus2.0: Microbiome Analysis Dashboard

![App logo](https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot%202025-09-04%20223841.png)

**A Shiny-based interactive dashboard for analyzing and visualizing microbiome data.**
**Built with R, Shiny, and phyloseq. Runs locally, so your data never leaves your machine.**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21287079.svg)](https://doi.org/10.5281/zenodo.21287079)
[![R tests](https://github.com/barah123/canis_lupus2.0/actions/workflows/r-tests.yml/badge.svg)](https://github.com/barah123/canis_lupus2.0/actions/workflows/r-tests.yml)
![License: MIT](https://img.shields.io/badge/license-MIT-green)

---

## 📌 Table of Contents
- [Project Overview](#project-overview)
- [Features](#features)
- [How to Use](#how-to-use)
- [Key Analyses](#key-analyses)
- [Testing](#testing)
- [Recognition](#recognition)
- [How to Cite](#how-to-cite)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## 🧪 Project Overview <a name="project-overview"></a>
CanisLupus2.0 is a user-friendly, web-based tool designed for researchers and bioinformaticians to explore, analyze, and visualize microbiome datasets. Built with R, Shiny, and phyloseq, this app allows users to upload their own amplicon sequencing data (ASV/OTU tables, taxonomy, metadata, and phylogenetic trees) **or** shotgun metagenomics taxonomic profiles (e.g., MetaPhlAn, Kraken2/Bracken), and perform comprehensive analyses without writing any code.

- **Author**: Philip Yamoah Appiah
- **Affiliated Institution**: Student at George Washington University – MS Health Data Science
- **Tools**: R, Shiny, phyloseq, vegan

---

## 🚀 Features <a name="features"></a>
- ✅ **Open the app in R**: Run the `app.R` file in R and the Shiny app will open in a browser
- ✅ **Two data types**: Amplicon (ASV/OTU table + taxonomy) or shotgun taxonomic profiles, selectable from the Data tab
- ✅ **Preprocessing & QC**: Singleton removal and minimum-read-depth filtering (automatically skipped for relative-abundance profiles), with a sequencing-depth plot and downloadable dataset summary
- ✅ **Visual Exploration**: Interactive bar plots, pie charts, rarefaction curves, phylogenetic trees, and metacoder heat trees
- ✅ **Community Analysis**: Alpha/beta diversity, PERMANOVA, differential abundance, core microbiome, abundance heatmaps, and co-occurrence networks
- ✅ **Repeated-measures aware**: Automatically switches to a Friedman test instead of Kruskal-Wallis when a complete blocked/repeated-measures design is detected
- ✅ **User-Friendly UI**: Intuitive interface with themed styling
- ✅ **Runs Locally**: No server upload required, so unpublished data stays on your machine
- ✅ **Tested**: Automated `testthat` suite covering core helpers, run on every push via GitHub Actions

---

## 📖 How to Use <a name="how-to-use"></a>
1. Install R and the required libraries (see below), then run `app.R` — the Shiny app opens in your browser
2. On the **Data** tab, choose a data type:
   - **Amplicon (ASV/OTU table)** — upload an abundance table and taxonomy table (plus optional metadata and phylogenetic tree)
   - **Shotgun (taxonomic profile)** — upload one merged profile from MetaPhlAn, Kraken2/Bracken, or similar
   - Leave the upload fields empty to explore a bundled example dataset instead
3. Click **Load & process data** and review the sequencing-depth plot and dataset summary
4. Explore the **Diversity**, **Taxonomy**, **Phylogenetics**, **Networks**, and **Statistics** menus with different taxonomic levels, grouping variables, and statistical parameters

### Sample datasets
Three example datasets ship with the app so you can try it before uploading your own data:

| Dataset | Folder | Data type | Notes |
|---|---|---|---|
| Human skin microbiome | `skin_data/` | Amplicon (ASV table) | Body-site comparison, used for validation |
| Inflammatory bowel disease (IBD) | `IBD_data/` | Amplicon (ASV table) | Includes a phylogenetic tree |
| Colorectal cancer vs. control / HMP body sites | `sample_data_shotgun/` | Shotgun (MetaPhlAn profiles) | See [`sample_data_shotgun/README.md`](sample_data_shotgun/README.md) for details |

<div align="center">
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot%202025-11-04%20134710.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot2.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot3.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot4.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot (13).png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot (11).png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot 2026-03-12 210309.jpg" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot 2026-03-12 210529.png" width="45%" />
</div>

---

## 🔬 Key Analyses <a name="key-analyses"></a>
- ✅ **Taxonomic profiling**: interactive stacked bar plots, pie charts, and metacoder heat trees at user-selected taxonomic levels
- ✅ **Alpha diversity**: Observed richness, Chao1, ACE, Shannon, Simpson, InvSimpson, and Fisher indices, with per-sample, by-group, and group-comparison (Kruskal-Wallis or Friedman) views
- ✅ **Beta diversity**: PCoA, NMDS, RDA, and CCA ordination over Bray-Curtis, Jaccard, (weighted) UniFrac, or Euclidean distance, with PERMANOVA group testing
- ✅ **Phylogenetics**: interactive phylogenetic tree viewer (rectangular, circular, fan, or radial layout)
- ✅ **Core microbiome**: detection of taxa consistently present across samples at user-defined abundance and prevalence thresholds
- ✅ **Networks & clustering**: CLR-normalised abundance heatmaps, hierarchical clustering/dendrograms, and taxon co-occurrence networks
- ✅ **Differential abundance**: Kruskal-Wallis (or Friedman for repeated measures) per taxon with Benjamini-Hochberg correction, shown as volcano plots and results tables
- ✅ **Data transformations**: CLR and compositional transformation views (ordination and distribution)

---

### Required R Libraries
To run CanisLupus2.0, ensure the following libraries are installed in your R environment:
```r
  library(shiny)
  library(phyloseq)
  library(microbiome)
  library(plotly)
  library(DT)
  library(vegan)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(tidyverse)
  library(metacoder)
  library(networkD3)
  library(heatmaply)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(igraph)
  library(cowplot)
  library(ggdendro)
  library(dendextend)
  library(reshape2)
  library(scales)
```
---

## ✅ Testing <a name="testing"></a>
CanisLupus2.0 has an automated `testthat` suite covering key helpers (e.g., the amplicon/shotgun data coercion, the phylogenetic/UniFrac guard, and the repeated-measures test-selection logic). Run it locally with:
```r
testthat::test_dir("tests/testthat")
```
Tests run automatically on every push and pull request via GitHub Actions (see [`.github/workflows/r-tests.yml`](.github/workflows/r-tests.yml)).

---

## Recognition <a name="recognition"></a>

- 🏆 **1st Place, George Washington University Open Source Software Awards**
- 🎤 Presented at the **GW Open Source Program Office (OSPO) Conference**, 2026
- 📋 Accepted as a **poster presentation** at **ISMB 2026** (Intelligent Systems
  for Molecular Biology), Microbiome COSI track, hosted by the International
  Society for Computational Biology — Washington, D.C.

---

## 📑 How to Cite <a name="how-to-cite"></a>
If you use CanisLupus2.0 in your work, please cite the software release:

> Appiah, P., Fadaq, A., Banda, E., Perez-Losada, M., Ortega-Santos, C., & Crandall, K. (2026). *CanisLupus 2.0: An Interactive R Shiny Application for Microbiome Data Analysis* (Version v2.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.21287079

Machine-readable citation metadata is also available in [`CITATION.cff`](CITATION.cff).

---

## 🤝 Contributing <a name="contributing"></a>
Contributions are welcome! See [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to report bugs, request features, and submit changes, including the code style guide and how to run the test suite before opening a pull request.

---

## License <a name="license"></a>

Copyright (c) 2026 Philip Yamoah Appiah

This project is licensed under the terms of the [MIT license](LICENSE).
You are free to use, modify, and distribute this work, subject to the
conditions specified in the LICENSE file.

---

### ✉️ Contact <a name="Contact"></a>
For questions or collaboration, please contact:
📧 [pyappiah561@gmail.com](mailto:pyappiah561@gmail.com)

![App logo_2](https://github.com/barah123/canis_lupus2.0/blob/main/www/canis_logo.png)
