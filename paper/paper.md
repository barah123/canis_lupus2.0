---
title: 'CanisLupus 2.0: An Interactive R Shiny Application for Microbiome Data Analysis'
tags:
  - R
  - Shiny
  - microbiome
  - bioinformatics
  - amplicon sequencing
  - diversity analysis
authors:
  - name: Philip Y. Appiah
    orcid: 0009-0002-0706-2506
    affiliation: 1
affiliations:
  - name: Department of Biostatistics and Bioinformatics, George Washington University
    index: 1
date: 05 July 2026
bibliography: paper.bib
---

# Summary

Microbiome research increasingly relies on amplicon sequencing data (e.g., 16S rRNA, ITS)
to characterize microbial community composition and structure across sample groups.
While this data type is common, downstream analysis typically requires combining
multiple specialized R packages with substantial custom scripting, creating a
barrier for researchers without a programming background. `CanisLupus 2.0` is an
R Shiny web application that provides an integrated, no-code interface for a
standard microbiome exploratory analysis workflow: users upload amplicon sequence
variant (ASV) tables, taxonomy tables, sample metadata, and phylogenetic trees,
and the application handles preprocessing and interactive visualization through
a single dashboard.

# Statement of need

Many microbiome analysis pipelines exist as command-line R or Python packages that
assume fluency in scripting, package-specific data structures (e.g., `phyloseq`
objects), and statistical methodology. This creates a steep barrier for wet-lab
researchers, clinicians, and students who generate or need to interpret microbiome
data but lack a computational background. Existing web-based microbiome tools
often support only a subset of the standard exploratory workflow (e.g., taxonomic
summaries alone) or require submitting data to a third-party server, raising data
privacy concerns for unpublished or sensitive datasets.

`CanisLupus 2.0` addresses this gap by (1) supporting a broad set of standard
exploratory analyses — preprocessing, taxonomic profiling, alpha and beta
diversity, phylogenetic visualization, rarefaction, core microbiome
identification, and correlation network analysis — in a single interface;
(2) running locally or on a self-hosted server so data never needs to leave the
user's institution; and (3) requiring no R programming knowledge to operate, while
remaining fully open source so computationally inclined users can extend or audit
the underlying methods.

# State of the field

Several existing tools address parts of the microbiome analysis workflow.
Command-line R packages such as `phyloseq` [@mcmurdie2013phyloseq] and `vegan`
[@oksanen2013package] provide the underlying statistical and ecological methods
that most downstream tools, including `CanisLupus 2.0`, build on, but require
scripting proficiency to use directly. Packages such as `metacoder`
[@foster2017metacoder] extend this further with hierarchical taxonomic
visualization ("heat trees"), and `ggtree` [@yu2017ggtree] provides
publication-quality phylogenetic tree annotation, both of which require
familiarity with R's plotting and data-manipulation idioms. Web-based platforms
(e.g., hosted microbiome analysis portals) lower this barrier but typically
require uploading data to a third-party server and often expose only a
constrained subset of possible analyses. `CanisLupus 2.0` differentiates itself
by combining a broad set of standard exploratory analyses in a single
self-hostable dashboard that requires no scripting, while remaining fully open
source.

# Functionality

`CanisLupus 2.0` is built in R using the `shiny` framework [@chang2015shiny]
and the `phyloseq` package [@mcmurdie2013phyloseq] as its core data structure,
with phylogenetic handling via `ape` [@paradis2019ape] and `ggtree`
[@yu2017ggtree], taxonomic visualization via `metacoder` [@foster2017metacoder],
diversity statistics via `vegan` [@oksanen2013package], and general data
manipulation via the `tidyverse` [@wickham2019tidyverse]. Key features include:

- **Data input and preprocessing**: upload of ASV tables, taxonomy tables, sample
  metadata, and phylogenetic tree files.
- **Taxonomic profiling**: interactive stacked bar plots, pie charts, and
  heatmaps of community composition at user-selected taxonomic levels.
- **Alpha diversity**: Shannon and Simpson indices with group comparisons.
- **Beta diversity**: PCoA ordination using Bray-Curtis distance.
- **Phylogenetic visualization**: tree plots and rendered phylogenies via
  `ggtree`.
- **Rarefaction curves**: assessment of sequencing depth sufficiency.
- **Core microbiome identification**: detection of taxa consistently present
  across samples at user-defined abundance and prevalence thresholds.
- **Correlation networks**: visualization of taxon co-occurrence patterns.
- **Sample clustering**: hierarchical clustering and dendrogram visualization.

The application is Dockerized for reproducible deployment and ships with
example datasets (skin and IBD microbiome data) so new users can explore its
functionality before uploading their own data.

# Research impact statement

`CanisLupus 2.0` has been used to analyze gut microbiome data for the Fit Gut
Lab at George Washington University and formed the basis of the author's
thesis research. The software received first place at the George Washington
University Open Source Software Awards, was presented at the GW Open Source
Program Office (OSPO) Conference in 2026, and was accepted as a poster
presentation at the Intelligent Systems for Molecular Biology (ISMB) 2026
conference, Microbiome COSI track, hosted by the International Society for
Computational Biology, in Washington, D.C.

# AI usage disclosure

The application code for `CanisLupus 2.0` was written entirely by the author
without the use of AI coding assistants. AI assistance (Claude, Anthropic) was
used to help draft supporting project documentation for this submission,
including this paper, `CITATION.cff`, `CONTRIBUTING.md`, and scaffolding for
the automated test suite; all AI-assisted content was reviewed and edited by
the author before inclusion.

# Acknowledgements

The author thanks the George Washington University Department of Biostatistics
and Bioinformatics, and the Fit Gut Lab, for support during the development and
research application of this software.

# References
