## ----setup,include=FALSE------------------------------------------------------
library(BiocStyle)

## ----loadPackages,message=FALSE,warning=FALSE---------------------------------
library(Biostrings) 
library(hgu95av2probe) 
library(hgu95av2cdf) 

## ----hgu95av2probe,eval=FALSE-------------------------------------------------
# ?hgu95av2probe

## ----reverseComplement,eval=FALSE---------------------------------------------
# ?reverseComplement

## ----MatchingSetsofProbesAgasintEachOther-------------------------------------
pm <- DNAStringSet(hgu95av2probe) 
dict <- pm[3801:4000]
pdict <- PDict(dict)
m <- vcountPDict(pdict, pm)
dim(m) 
table(rowSums(m))
which(rowSums(m) == 3)
ii <- which(m[77, ] != 0)
pm[ii]

## ----alphabetFrequency--------------------------------------------------------
bcpm <- alphabetFrequency(pm, baseOnly=TRUE)
head(bcpm) 
alphabetFrequency(pm, baseOnly=TRUE, collapse=TRUE)

## ----hgu95av2dimncol----------------------------------------------------------
nc = hgu95av2dim$NCOL
nc

## ----hgu95av2dimnrow----------------------------------------------------------
nr = hgu95av2dim$NROW
nr

## ----abseq--------------------------------------------------------------------
library(affy) 
abseq = rep(as.character(NA), nc*nr) 
ipm = with(hgu95av2probe, xy2indices(x, y, nc=nc)) 
any(duplicated(ipm)) # just a sanity check 
abseq[ipm] = hgu95av2probe$sequence
table(is.na(abseq))

## ----pm2mm--------------------------------------------------------------------
mm <- pm
subseq(mm, start=13, width=1) <- complement(subseq(mm, start=13, width=1))
cat(as.character(pm[[1]]), as.character(mm[[1]]), sep="\n")

## ----imm----------------------------------------------------------------------
imm = with(hgu95av2probe, xy2indices(x, y+1, nc=nc))
intersect(ipm, imm) # just a sanity check
abseq[imm] = as.character(mm)
table(is.na(abseq))

## ----alphabetFrequency2-------------------------------------------------------
freqs <- alphabetFrequency(DNAStringSet(abseq[!is.na(abseq)]), baseOnly=TRUE)
bc <- matrix(nrow=length(abseq), ncol=5)
colnames(bc) <- colnames(freqs)
bc[!is.na(abseq), ] <- freqs
head(na.omit(bc))

## ----gc-----------------------------------------------------------------------
GC = ordered(bc[,"G"] + bc[,"C"])
colores = rainbow(nlevels(GC))

## ----bap, fig.cap="Distribution of probe GC content. The height of each bar corresponds to the number of probes with the corresponding GC content."----
library(affydata)
f <- system.file("extracelfiles", "CL2001032020AA.cel", package="affydata")
pd <- new("AnnotatedDataFrame", data=data.frame(fromFile=I(f), row.names="f"))
abatch <- read.affybatch(filenames=f, compress=TRUE, phenoData=pd)
barplot(table(GC), col=colores, xlab="GC", ylab="number")

## ----bxp, fig.cap="Boxplots of log~2~ intensity stratifed by probe GC content."----
boxplot(log2(exprs(abatch)[,1]) ~ GC, outline=FALSE,
        col=colores, , xlab="GC", ylab=expression(log[2]~intensity))

## ----p2p, fig.cap="Scatterplot of PM vs MM intensities, colored by probe GC content."----
plot(exprs(abatch)[ipm,1], exprs(abatch)[imm,1], asp=1, pch=".", log="xy",
     xlab="PM", ylab="MM", col=colores[GC[ipm]])
abline(a=0, b=1, col="#404040", lty=3)

## ----devoff,include=FALSE-----------------------------------------------------
dev.off()

