## ----objectCreation, message=FALSE--------------------------------------------
library(Biostrings)
origMAlign <-
  readDNAMultipleAlignment(filepath =
                           system.file("extdata",
                                       "msx2_mRNA.aln",
                                       package="Biostrings"),
                           format="clustal")

phylipMAlign <-
  readAAMultipleAlignment(filepath =
                          system.file("extdata",
                                      "Phylip.txt",
                                      package="Biostrings"),
                          format="phylip")

## ----renameRows---------------------------------------------------------------
rownames(origMAlign)
rownames(origMAlign) <- c("Human","Chimp","Cow","Mouse","Rat",
                          "Dog","Chicken","Salmon")
origMAlign

## ----detail, eval=FALSE-------------------------------------------------------
# detail(origMAlign)

## ----usingMasks---------------------------------------------------------------
maskTest <- origMAlign
rowmask(maskTest) <- IRanges(start=1,end=3)
rowmask(maskTest)
maskTest

colmask(maskTest) <- IRanges(start=c(1,1000),end=c(500,2343))
colmask(maskTest)
maskTest

## ----nullOut masks------------------------------------------------------------
rowmask(maskTest) <- NULL
rowmask(maskTest)
colmask(maskTest) <- NULL
colmask(maskTest)
maskTest

## ----invertMask---------------------------------------------------------------
rowmask(maskTest, invert=TRUE) <- IRanges(start=4,end=8)
rowmask(maskTest)
maskTest
colmask(maskTest, invert=TRUE) <- IRanges(start=501,end=999)
colmask(maskTest)
maskTest

## ----setup--------------------------------------------------------------------
## 1st lets null out the masks so we can have a fresh start.
colmask(maskTest) <- NULL
rowmask(maskTest) <- NULL

## ----appendMask---------------------------------------------------------------
## Then we can demonstrate how the append argument works
rowmask(maskTest) <- IRanges(start=1,end=3)
maskTest

rowmask(maskTest,append="intersect") <- IRanges(start=2,end=5)
maskTest

rowmask(maskTest,append="replace") <- IRanges(start=5,end=8)
maskTest

rowmask(maskTest,append="replace",invert=TRUE) <- IRanges(start=5,end=8)
maskTest

rowmask(maskTest,append="union") <- IRanges(start=7,end=8)
maskTest

## ----maskMotif----------------------------------------------------------------
tataMasked <- maskMotif(origMAlign, "TATA")
colmask(tataMasked)

## ----maskGaps-----------------------------------------------------------------
autoMasked <- maskGaps(origMAlign, min.fraction=0.5, min.block.width=4)
autoMasked

## ----asmatrix-----------------------------------------------------------------
full = as.matrix(origMAlign)
dim(full)
partial = as.matrix(autoMasked)
dim(partial)

## ----alphabetFreq-------------------------------------------------------------
alphabetFrequency(autoMasked)

## ----consensus----------------------------------------------------------------
consensusMatrix(autoMasked, baseOnly=TRUE)[, 84:90]
substr(consensusString(autoMasked),80,130)
consensusViews(autoMasked)

## ----cluster------------------------------------------------------------------
sdist <- stringDist(as(origMAlign,"DNAStringSet"), method="hamming")
clust <- hclust(sdist, method = "single")
png(file="badTree.png")
plot(clust)
dev.off()

## ----figure, echo=FALSE, fig=TRUE, eps=FALSE, fig.align = 'center', fig.cap='Funky tree produced by using unmasked strings'----
knitr::include_graphics("badTree.png")

## ----cluster2-----------------------------------------------------------------
sdist <- stringDist(as(autoMasked,"DNAStringSet"), method="hamming")
clust <- hclust(sdist, method = "single")
png(file="goodTree.png")
plot(clust)
dev.off()
fourgroups <- cutree(clust, 4)
fourgroups

## ----figure1, echo=FALSE, fig=TRUE, eps=FALSE, width=0.6,fig.align = 'center', height=5, fig.cap='A tree produced by using strings with masked gaps'----
knitr::include_graphics("goodTree.png")

## ----fastaExample, eval=FALSE-------------------------------------------------
# DNAStr = as(origMAlign, "DNAStringSet")
# writeXStringSet(DNAStr, file="myFile.fa")

## ----write.phylip, eval=FALSE-------------------------------------------------
# write.phylip(phylipMAlign, filepath="myFile.txt")

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

