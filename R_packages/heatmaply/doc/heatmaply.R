## ----echo = FALSE, message = FALSE--------------------------------------------
library("heatmaply")
library("knitr")
knitr::opts_chunk$set(
  warning=FALSE, # because of a bug in plotly: https://github.com/ropensci/plotly/issues/1670
  # see also: https://github.com/talgalili/heatmaply/issues/226
  # cache = TRUE,
  dpi = 60,
  comment = "#>",
  tidy = FALSE)

# http://stackoverflow.com/questions/24091735/why-pandoc-does-not-retrieve-the-image-file
# < ! -- rmarkdown v1 -->


## ----eval = FALSE-------------------------------------------------------------
#  install.packages('heatmaply')

## ----eval = FALSE-------------------------------------------------------------
#  # You'll need devtools
#  install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
#  install.packages.2('remotes')
#  
#  remotes::install_github("ropensci/plotly")
#  remotes::install_github('talgalili/heatmaply')

## -----------------------------------------------------------------------------
library("heatmaply")

## -----------------------------------------------------------------------------
library("heatmaply")
heatmaply(mtcars)

## -----------------------------------------------------------------------------
heatmaply(mtcars, xlab = "Features", ylab = "Cars", 
  main = "An example of title and xlab/ylab")

## -----------------------------------------------------------------------------
heatmaply_cor(
  cor(mtcars),
  xlab = "Features",
  ylab = "Features",
  k_col = 2,
  k_row = 2
)

## -----------------------------------------------------------------------------

r <- cor(mtcars)
## We use this function to calculate a matrix of p-values from correlation tests
## https://stackoverflow.com/a/13112337/4747043
cor.test.p <- function(x){
    FUN <- function(x, y) cor.test(x, y)[["p.value"]]
    z <- outer(
      colnames(x), 
      colnames(x), 
      Vectorize(function(i,j) FUN(x[,i], x[,j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}
p <- cor.test.p(mtcars)

heatmaply_cor(
  r,
  node_type = "scatter",
  point_size_mat = -log10(p), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
)

## ----eval = F-----------------------------------------------------------------
#  heatmaply(
#    mtcars,
#    xlab = "Features",
#    ylab = "Cars",
#    scale = "column",
#    main = "Data transformation using 'scale'"
#  )

## ----eval = F-----------------------------------------------------------------
#  heatmaply(
#    normalize(mtcars),
#    xlab = "Features",
#    ylab = "Cars",
#    main = "Data transformation using 'normalize'"
#  )

## -----------------------------------------------------------------------------
heatmaply(
  percentize(mtcars),
  xlab = "Features",
  ylab = "Cars", 
  main = "Data transformation using 'percentize'"
)

## -----------------------------------------------------------------------------
heatmaply_na(
  airquality[1:30, ],
  showticklabels = c(TRUE, FALSE),
  k_col = 3,
  k_row = 3
)


# warning - using grid_color cannot handle a large matrix! For example:
# airquality[1:10, ] %>% is.na10 %>% 
#   heatmaply(color = c("white", "black"), grid_color = "grey",
#             k_col =3, k_row = 3,
#             margins = c(40, 50)) 
# airquality %>% is.na10 %>% 
#   heatmaply(color = c("grey80", "grey20"), # grid_color = "grey",
#             k_col =3, k_row = 3,
#             margins = c(40, 50)) 
# 


## ----eval = F-----------------------------------------------------------------
#  heatmaply_cor(
#    cor(mtcars),
#    k_col = 2,
#    k_row = 2
#  )
#  

## ----eval = FALSE-------------------------------------------------------------
#  heatmaply(
#    percentize(mtcars),
#    colors = heat.colors(100)
#  )

## -----------------------------------------------------------------------------
heatmaply(
  mtcars,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    midpoint = 200, 
    limits = c(0, 500)
  )
)


## -----------------------------------------------------------------------------
# The default of heatmaply:
heatmaply(
  percentize(mtcars)[1:10, ],
  seriate = "OLO"
)

## ----eval = F-----------------------------------------------------------------
#  # Similar to OLO but less optimal (since it is a heuristic)
#  heatmaply(
#    percentize(mtcars)[1:10, ],
#    seriate = "GW"
#  )

## ----eval = F-----------------------------------------------------------------
#  # the default by gplots::heatmaply.2
#  heatmaply(
#    percentize(mtcars)[1:10, ],
#    seriate = "mean"
#  )

## ----eval = F-----------------------------------------------------------------
#  # the default output from hclust
#  heatmaply(
#    percentize(mtcars)[1:10, ],
#    seriate = "none"
#  )

## -----------------------------------------------------------------------------

x  <- as.matrix(datasets::mtcars)

# now let's spice up the dendrograms a bit:
library("dendextend")

row_dend  <- x %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram %>%
  set("branches_k_color", k = 3) %>% 
  set("branches_lwd", c(1, 3)) %>%
  ladderize
# rotate_DendSer(ser_weight = dist(x))
col_dend  <- x %>% 
  t %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram %>%
  set("branches_k_color", k = 2) %>% 
  set("branches_lwd", c(1, 2)) %>%
  ladderize
#    rotate_DendSer(ser_weight = dist(t(x)))

heatmaply(
  percentize(x),
  Rowv = row_dend,
  Colv = col_dend
)



## -----------------------------------------------------------------------------
x  <- as.matrix(datasets::mtcars)
gplots::heatmap.2(
  x,
  trace = "none",
  col = viridis(100),
  key = FALSE
)

## ----eval = FALSE-------------------------------------------------------------
#  heatmaply(
#    x,
#    seriate = "mean"
#  )

## -----------------------------------------------------------------------------
heatmaply(x,
  seriate = "mean", 
  row_dend_left = TRUE,
  plot_method = "plotly"
)

## -----------------------------------------------------------------------------
# Example for using RowSideColors

x  <- as.matrix(datasets::mtcars)
rc <- colorspace::rainbow_hcl(nrow(x))

library("gplots")
library("viridis")
heatmap.2(
  x,
  trace = "none",
  col = viridis(100),
  RowSideColors = rc,
  key = FALSE
)


## ----eval = FALSE-------------------------------------------------------------
#  heatmaply(
#    x,
#    seriate = "mean",
#    RowSideColors = rc
#  )

## -----------------------------------------------------------------------------
heatmaply(
  x[, -c(8, 9)],
  seriate = "mean",
  col_side_colors = c(rep(0, 5), rep(1, 4)),
  row_side_colors = x[, 8:9]
)


## -----------------------------------------------------------------------------
heatmaply(
  mtcars,
  cellnote = mtcars
)

## -----------------------------------------------------------------------------
mat <- mtcars
mat[] <- paste("This cell is", rownames(mat))
mat[] <- lapply(colnames(mat), function(colname) {
    paste0(mat[, colname], ", ", colname)
})
heatmaply(
  mtcars,
  custom_hovertext = mat
)

## -----------------------------------------------------------------------------
ggheatmap(
  mtcars,
  scale = "column",
  row_side_colors = mtcars[, c("cyl", "gear")]
)

## ----eval = FALSE-------------------------------------------------------------
#  dir.create("folder")
#  heatmaply(mtcars, file = "folder/heatmaply_plot.html")
#  browseURL("folder/heatmaply_plot.html")

## ----eval = FALSE-------------------------------------------------------------
#  dir.create("folder")
#  # Before the first time using this code you may need to first run:
#  # webshot::install_phantomjs() or to install
#  # [plotly's orca](https://github.com/plotly/orca) program.
#  heatmaply(mtcars, file = "folder/heatmaply_plot.png")
#  browseURL("folder/heatmaply_plot.png")

## ----eval = FALSE-------------------------------------------------------------
#  # This saves the file, but does not plot it in the RStudio viewer
#  tmp <- heatmaply(mtcars, file = "folder/heatmaply_plot.png")
#  rm(tmp)

## ----eval = FALSE-------------------------------------------------------------
#  
#  library("microbenchmark")
#  
#  
#  library("heatmaply")
#  x <- matrix(1:1000, 500, 2)
#  
#  microbenchmark(
#    heatmaply(x, hclustfun = stats::hclust),
#    heatmaply(x, hclustfun = fastcluster::hclust),
#    times = 10
#  )
#  
#  x <- matrix(1:1000, 1000, 2)
#  microbenchmark(
#    stats::hclust(dist(x)),
#    fastcluster::hclust(dist(x)),
#    times = 10
#  )
#  
#  

## ----cache = FALSE------------------------------------------------------------
sessionInfo()

