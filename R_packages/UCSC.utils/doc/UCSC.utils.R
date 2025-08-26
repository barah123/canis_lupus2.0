## ----eval=FALSE---------------------------------------------------------------
#  if (!require("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("UCSC.utils")

## ----list_UCSC_genomes--------------------------------------------------------
suppressPackageStartupMessages(library(UCSC.utils))

list_UCSC_genomes("cat")

## ----get_UCSC_chrom_sizes-----------------------------------------------------
felCat9_chrom_sizes <- get_UCSC_chrom_sizes("felCat9")
head(felCat9_chrom_sizes)

## ----list_UCSC_tracks---------------------------------------------------------
list_UCSC_tracks("felCat9", group="varRep")

## ----fetch_UCSC_track_data----------------------------------------------------
mm9_cytoBandIdeo <- fetch_UCSC_track_data("mm9", "cytoBandIdeo")
head(mm9_cytoBandIdeo)

## ----UCSC_dbselect------------------------------------------------------------
felCat9_refGene <- UCSC_dbselect("felCat9", "refGene")
head(felCat9_refGene)

## ----UCSC_dbselect2-----------------------------------------------------------
columns <- c("chrom", "strand", "txStart", "txEnd", "exonCount", "name2")
UCSC_dbselect("felCat9", "refGene", columns=columns, where="chrom='chrA1'")

## -----------------------------------------------------------------------------
sessionInfo()

