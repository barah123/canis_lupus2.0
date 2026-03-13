# 🧬 CanisLupus2.0: Microbiome Analysis Dashboard

![App logo](https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot%202025-09-04%20223841.png)

**A Shiny-based interactive dashboard for analyzing and visualizing microbiome data.**
**Built with R, Shiny, phyloseq, and Dockerised for reproducibility and ease of use.**

---

## 📌 Table of Contents
- [Project Overview](#project-overview)
- [Features](#features)
- [How to Use](#how-to-use)
- [Key Analyses](#key-analyses)
- [License](#license)
- [Contact](#contact)

---

## 🧪 Project Overview <a name="project-overview"></a>
CanisLupus2.0 is a user-friendly, web-based tool designed for researchers and bioinformaticians to explore, analyze, and visualize microbiome datasets. Built with R, Shiny, and phyloseq, this app allows users to upload their own data (ASV tables, taxonomy, metadata, and phylogenetic trees) and perform comprehensive analyses.

- **Author**: Philip Yamoah Appiah
- **Affiliated Institution**: Student at George Washington University – MS Health Data Science
- **Tools**: R Shiny, phyloseq, Docker

---

## 🚀 Features <a name="features"></a>
- ✅ **Open the app in R**: Run the `app.R` file in R and the Shiny app will open in a browser
- ✅ **Data Upload**: Upload ASV, taxonomy, metadata, and phylogenetic tree files
- ✅ **Visual Exploration**: Interactive bar plots, pie charts, rarefaction curves, and phylogenetic trees
- ✅ **Community Analysis**: Alpha/beta diversity, core microbiome, and correlation networks
- ✅ **User-Friendly UI**: Intuitive interface with themed styling
- ✅ **Docker Support**: Pre-configured Docker image for easy deployment

---

## 📖 How to Use <a name="how-to-use"></a>
1. Run `app.R` in R
2. Upload ASV, taxonomy, metadata, and phylogenetic tree files (sample dataset available in the `skin_data` folder)
3. Click on the **Load** button and wait for the summary statistics
4. Explore the **Visual Exploration**, **Community Profile**, and **Network Analysis** with different taxonomic levels and features


<div align="center">
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot%202025-11-04%20134710.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot2.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot3.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot4.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot (13).png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/newplot (11).png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot 2026-03-12 210309.png" width="45%" />
  <img src="https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot 2026-03-12 210529.png" width="45%" />
</div>
---

## 🔬 Key Analyses <a name="key-analyses"></a>
- ✅ Taxonomic profiling at phylum/genus levels
- ✅ Alpha diversity (Shannon, Simpson indices)
- ✅ Beta diversity (Bray-Curtis distance, PCoA plots)
- ✅ Visualization of taxa abundance (bar plots, heatmaps)
- ✅ Sample clustering and dendrograms

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
![License: MIT](https://img.shields.io/badge/license-MIT-green)

Copyright (c) [2025] [Philip Y. Appiah]

This project is licensed under the terms of the MIT license.
You are free to use, modify, and distribute this work, subject to the
conditions specified in the LICENSE file.

---

          
### ✉️ Contact <a name="Contact"></a>
For questions or collaboration, please contact:
📧 [pyappiah561@gmail.com


![App logo_2](https://github.com/barah123/canis_lupus2.0/blob/main/www/canis_logo.png)

                      
