## ----install, eval=FALSE------------------------------------------------------
#  if (!require("BiocManager"))
#    install.packages("BiocManager")
#  BiocManager::install("S4Vectors")

## ----message=FALSE------------------------------------------------------------
library(S4Vectors)

## ----Rle-extends-Vector-------------------------------------------------------
showClass("Rle")

## ----initialize---------------------------------------------------------------
set.seed(0)
lambda <- c(rep(0.001, 4500), seq(0.001, 10, length=500),
            seq(10, 0.001, length=500))
xVector <- rpois(1e7, lambda)
yVector <- rpois(1e7, lambda[c(251:length(lambda), 1:250)])
xRle <- Rle(xVector)
yRle <- Rle(yVector)

## ----basic-ops----------------------------------------------------------------
length(xRle)
xRle[1]
zRle <- c(xRle, yRle)

## ----seq-extraction-----------------------------------------------------------
xSnippet <- window(xRle, 4751, 4760)
xSnippet
head(xSnippet)
tail(xSnippet)
rev(xSnippet)
rep(xSnippet, 2)
subset(xSnippet, xSnippet >= 5L)

## ----seq-concatenate----------------------------------------------------------
c(xSnippet, rev(xSnippet))
append(xSnippet, xSnippet, after=3)

## ----aggregate----------------------------------------------------------------
xSnippet
aggregate(xSnippet, start=1:8, width=3, FUN=median)

## ----shiftApply-cor-----------------------------------------------------------
cor(xRle, yRle)
shifts <- seq(235, 265, by=3)
corrs  <- shiftApply(shifts, yRle, xRle, FUN=cor)

## ----figshiftcorrs, eps=FALSE, fig.align='center', fig.cap='Correlation between `xRle` and `yRle` for various shifts'----
plot(shifts, corrs)

## ----Rle-vector-compare-------------------------------------------------------
as.vector(object.size(xRle) / object.size(xVector))
identical(as.vector(xRle), xVector)

## ----Rle-accessors------------------------------------------------------------
head(runValue(xRle))
head(runLength(xRle))

## ----Rle-ops------------------------------------------------------------------
xRle > 0
xRle + yRle
xRle > 0 | yRle > 0

## ----Rle-summary--------------------------------------------------------------
range(xRle)
sum(xRle > 0 | yRle > 0)

## ----Rle-math-----------------------------------------------------------------
log1p(xRle)

## ----Rle-cor------------------------------------------------------------------
cor(xRle, yRle)
shiftApply(249:251, yRle, xRle,
           FUN=function(x, y) {var(x, y) / (sd(x) * sd(y))})

## ----DataFrame-extends-List---------------------------------------------------
showClass("DataFrame")

## ----DataFrame, warning=FALSE-------------------------------------------------
df <- DataFrame(x=xRle, y=yRle)
sapply(df, class)
sapply(df, summary)
sapply(as.data.frame(df), summary)
endoapply(df, `+`, 0.5)

## ----SessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

