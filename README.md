# 🧬 CanisLupus2.0: Microbiome Analysis Dashboard

**A Shiny-based interactive dashboard for visualizing and analyzing microbiome data.**  
**Built with R, Shiny, phyloseq, and Docker for reproducibility and ease of use.**

![CanisLupus Logo](canis1.png)

---

## 📌 Table of Contents
- [Project Overview](#project-overview)
- [Features](#features)
- [How to Use](#how-to-use)
- [Key Analyses](#key-analyses)
- [Technologies Used](#technologies-used)
- [Installation](#installation)
- [License](#license)
- [Contact](#contact)

---

## 🧪 Project Overview

CanisLupus2.0 is a user-friendly, web-based tool designed for researchers and bioinformaticians to explore, analyze, and visualize microbiome datasets. Built with R, Shiny, and phyloseq, this app allows users to upload their own data (ASV tables, taxonomy, metadata, and phylogenetic trees) and perform comprehensive analyses.

- **Author**: Philip Yamoah Appiah
- **Affiliated Institution**: Student at George Washington University – MS Health Data Science
- **Tools**: R Shiny, phyloseq, Docker

---

## 🚀 Features

✅ **Open the app in R**: Run the app.R file in R and the Shiny app will open in a browser  
✅ **Data Upload**: Upload ASV, taxonomy, metadata, and phylogenetic tree files  
✅ **Visual Exploration**: Interactive bar plots, pie charts, rarefaction curves, and phylogenetic trees  
✅ **Community Analysis**: Alpha/beta diversity, core microbiome, and correlation networks  
✅ **User-Friendly UI**: Intuitive interface with themed styling  
✅ **Docker Support**: Pre-configured Docker image for easy deployment  

---

## 📖 How to Use

1. Run `app.R` in R
![App interface](https://github.com/barah123/canis_lupus2.0/blob/main/www/Screenshot%202025-11-04%20134710.pn)
3. Upload ASV, taxonomy, metadata, and phylogenetic tree files (sample dataset available in the `skin_data` folder)
4. Click on the Load button and wait for the summary statistics
5. Explore the Visual Exploration, Community Profile and Network Analysis with different taxonomic levels and features

---

## 🔬 Key Analyses

✅ Taxonomic profiling at phylum/genus levels  
✅ Alpha diversity (Shannon, Simpson indices)  
✅ Beta diversity (Bray-Curtis distance, PCoA plots)  
✅ Visualization of taxa abundance (bar plots, heatmaps)  
✅ Sample clustering and dendrograms  

---

## 💻 Technologies Used

- `R Shiny`
- `Bioconductor-phyloseq`
- `Docker`

---

## 📥 Installation

### Load the required libraries in your R environment:

```r
library(shiny)
library(phyloseq)
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
