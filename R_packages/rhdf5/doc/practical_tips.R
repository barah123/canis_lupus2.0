## ----setup, echo = FALSE, include=FALSE---------------------------------------
set.seed(1234)

library(rhdf5)
library(dplyr)
library(ggplot2)
library(BiocParallel)

run_on_platform <- (Sys.info()['sysname'] == "Linux")

## ----create data, echo=TRUE, warning=FALSE------------------------------------
m1 <- matrix(rep(1:20000, each = 100), ncol = 20000, byrow = FALSE)
ex_file <- tempfile(fileext = ".h5")
h5write(m1, file = ex_file, name = "counts", level = 6)

## ----extract1, echo = TRUE----------------------------------------------------
system.time(
  res1 <- h5read(file = ex_file, name = "counts", 
                 index = list(NULL, 1:10000))
)

## ----extract2, echo = TRUE----------------------------------------------------
index <- list(NULL, seq(from = 1, to = 20000, by = 2))
system.time(
  res2 <- h5read(file = ex_file, name = "counts", 
                 index = index)
)

## ----extract3, echo = TRUE----------------------------------------------------
start <- c(1,1)
stride <- c(1,2)
block <- c(100,1)
count <- c(1,10000)
system.time(
  res3 <- h5read(file = ex_file, name = "counts", start = start,
                 stride = stride, block = block, count = count)
)
identical(res2, res3)

## ----chooseColunns, eval = TRUE-----------------------------------------------
columns <- sample(x = seq_len(20000), size = 1000, replace = FALSE) %>%
  sort()
f1 <- function(cols, name) { 
  h5read(file = ex_file, name = name, 
         index = list(NULL, cols))
  }

## ----singleReads, eval = run_on_platform, cache = TRUE------------------------
# system.time(res4 <- vapply(X = columns, FUN = f1,
#                            FUN.VALUE = integer(length = 100),
#                            name = 'counts'))

## ----createChunked, echo = TRUE, eval = TRUE, results='hide'------------------
h5createDataset(file = ex_file, dataset = "counts_chunked", 
                dims = dim(m1), storage.mode = "integer", 
                chunk = c(100,100), level = 6)
h5write(obj = m1, file = ex_file, name = "counts_chunked")

## ----read_chunked, eval = TRUE------------------------------------------------
system.time(res5 <- vapply(X = columns, FUN = f1, 
                           FUN.VALUE = integer(length = 100), 
                           name = 'counts_chunked'))

## -----------------------------------------------------------------------------
f2 <- function(block_size = 100) {
  cols_grouped <- split(columns,  (columns-1) %/% block_size)
  res <-  lapply(cols_grouped, f1, name = 'counts_chunked') %>%
    do.call('cbind', .)
}
system.time(f2())

## ----benchmark, echo = FALSE, cache = TRUE------------------------------------
bm <- bench::mark(
  f2(10), f2(25), f2(50), f2(100), 
  f2(250), f2(500), f2(1000), 
  f2(2000), f2(5000), f2(10000),
  iterations = 3, check = FALSE, time_unit = "s", memory = FALSE, filter_gc = FALSE
)

bm2 <- data.frame(block_size = gsub(".*\\(([0-9]+)\\)", "\\1", bm$expression) |>
                    as.integer() |>
                    rep(each = 3), 
                  time = unlist(bm$time))

## ----echo = FALSE, fig.width=6, fig.height=3, fig.wide = TRUE-----------------
ggplot(bm2, aes(x = block_size, y = time)) + 
  geom_point() + 
  scale_x_log10() +
  theme_bw() + 
  ylab('time (seconds)')

## ----hyperslab-benchmark, eval = run_on_platform, echo = FALSE, fig.width=6, fig.height=3, fig.wide = TRUE, fig.cap='The time taken to join hyperslabs increases expontentially with the number of join operations.  These timings are taken with no reading occuring, just the creation of a dataset selection.'----
# ## this code demonstrates the exponential increase in time as the
# ## number of hyberslab unions increases
# 
# select_index <- function(n = 1) {
# 
#   ## open the dataspace for the count table
#   fid <- H5Fopen(ex_file)
#   did  <- H5Dopen(fid, name = "counts")
#   sid <- H5Dget_space(did)
# 
#   ## column choice based on number of unions required
#   columns <- c(head(1:10001, n = -n), head(seq(10001-n+2, 20000, 2), n = n-1))
#   index <- list(100, columns)
#   H5Sselect_index(sid, index = index)
# 
#   ## tidy up
#   H5Sclose(sid)
#   H5Dclose(did)
#   H5Fclose(fid)
# }
# 
# bm <- bench::mark(
#   select_index(1), select_index(2), select_index(5),
#   select_index(10), select_index(20), select_index(50),
#   select_index(100), select_index(200), select_index(500),
#   select_index(1000), select_index(2000), select_index(5000),
#   select_index(10000),
#   iterations = 3, check = FALSE, time_unit = "s", memory = FALSE, filter_gc = FALSE
# )
# 
# bm2 <- data.frame(n = gsub(".*\\(([0-9]+)\\)", "\\1", bm$expression) |>
#                             as.integer() |>
#                             rep(each = 3),
#                   time = unlist(bm$time))
# 
# ggplot(bm2,aes(x = n, y = time)) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_bw() +
#   xlab('number of hyperslab unions') +
#   ylab('time (seconds)')

## ----example-dsets------------------------------------------------------------
dsets <- lapply(1:10, FUN = \(i) { matrix(runif(10000000), ncol = 100)} )
names(dsets) <- paste0("dset_", 1:10)

## -----------------------------------------------------------------------------
simple_writer <- function(file_name, dsets) {
  
  fid <- H5Fcreate(name = file_name)
  on.exit(H5Fclose(fid))
  
  for(i in seq_along(dsets)) {
    dset_name = paste0("dset_", i)
    h5createDataset(file = fid, dataset = dset_name, 
                    dims = dim(dsets[[i]]), chunk = c(10000, 10))
    h5writeDataset(dsets[[i]], h5loc = fid, name = dset_name)
  }
  
}

## -----------------------------------------------------------------------------
## Write a single dataset to a temporary file
## Arguments: 
## - dset_name: The name of the dataset to be created
## - dset: The dataset to be written
split_tmp_h5 <- function(dset_name, dset) {

  ## create a tempory HDF5 file for this dataset  
  file_name <- tempfile(pattern = "par", fileext = ".h5")
  fid <- H5Fcreate(file_name)
  on.exit(H5Fclose(fid))
  
  ## create and write the dataset
  ## we use some predefined chunk sizes 
  h5createDataset(file = fid, dataset = dset_name, 
                  dims = dim(dset), chunk = c(10000, 10))
  h5writeDataset(dset, h5loc = fid, name = dset_name)
  
  return(c(file_name, dset_name))
}

## Gather scattered datasets into a final single file
## Arguments: 
## - output_file: The path to the final HDF5 to be created
## - input: A data.frame with two columns containing the paths to the temp
## files and the name of the dataset inside that file
gather_tmp_h5 <- function(output_file, input) {
  
  ## create the output file
  fid <- H5Fcreate(name = output_file)
  on.exit(H5Fclose(fid))
  
  ## iterate over the temp files and copy the named dataset into our new file
  for(i in seq_len(nrow(input))) {
    fid2 <- H5Fopen(input$file[i])
    H5Ocopy(fid2, input$dset[i], h5loc_dest = fid, name_dest = input$dset[i])
    H5Fclose(fid2)
  }
  
}

## ----define-split-gather------------------------------------------------------
split_and_gather <- function(output_file, input_dsets, BPPARAM = NULL) {
  
  if(is.null(BPPARAM)) { BPPARAM <- BiocParallel::SerialParam() }
  
  ## write each of the matrices to a separate file
  tmp <- 
    bplapply(seq_along(input_dsets), 
           FUN = function(i) {
             split_tmp_h5(dset_name = names(input_dsets)[i], 
                          dset = input_dsets[[i]])
           }, 
           BPPARAM = BPPARAM) 
  
  ## create a table of file and the dataset names
  input_table <- do.call(rbind, tmp) |> 
    as.data.frame()
  names(input_table) <- c("file", "dset")
  
  ## copy all datasets from temp files in to final output
  gather_tmp_h5(output_file = output_file, input = input_table)
  
  ## remove the temporary files
  file.remove(input_table$file)
}

## ----eval = FALSE-------------------------------------------------------------
# split_and_gather(tempfile(), input_dsets = dsets,
#                  BPPARAM = MulticoreParam(workers = 2))

## ----run-writing-benchmark, eval = run_on_platform, cache = TRUE, message=FALSE, echo=FALSE----
# bench_results <- bench::mark(
#   "simple writer" = simple_writer(file_name = tempfile(), dsets = dsets),
#   "split/gather - 1 core" = split_and_gather(tempfile(), input_dsets = dsets,
#                                              BPPARAM = NULL),
#   "split/gather - 2 cores" = split_and_gather(tempfile(), input_dsets = dsets,
#                                               BPPARAM = MulticoreParam(workers = 2)),
#   "split/gather - 4 cores" = split_and_gather(tempfile(), input_dsets = dsets,
#                                               BPPARAM = MulticoreParam(workers = 4)),
#   iterations = 2, check = FALSE, time_unit = "s", memory = FALSE, filter_gc = FALSE
# )
# bench_results |> select(expression, min, median)

## ----fake-timings, eval = !run_on_platform, message = FALSE, echo = FALSE-----
  bench_results <- data.frame(median = 4:1)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

