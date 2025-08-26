## ----filter-pipeline, echo = FALSE, fig.cap="The default compression pipeline used by rhdf5"----
knitr::include_graphics("filter_pipeline.png")

## ----plugin-path, eval = TRUE-------------------------------------------------
rhdf5filters::hdf5_plugin_path()

## ----warning = FALSE, eval = FALSE--------------------------------------------
# ## blosc compressed file
# blosc_file <- system.file("h5examples/h5ex_d_blosc.h5",
#                           package = "rhdf5filters")

## ----h5dump-1, warning = FALSE, eval = FALSE----------------------------------
# h5dump_out <- system2('h5dump',
#                       args = c('-p', '-d /dset', blosc_file),
#                       stdout = TRUE, stderr = TRUE)
# cat(h5dump_out, sep = "\n")

## ----h5dump-1-out, eval = TRUE, echo = FALSE----------------------------------
cat(
'HDF5 "rhdf5filters/h5examples/h5ex_d_blosc.h5" {
DATASET "/dset" {
   DATATYPE  H5T_IEEE_F32LE
   DATASPACE  SIMPLE { ( 30, 10, 20 ) / ( 30, 10, 20 ) }
   STORAGE_LAYOUT {
      CHUNKED ( 10, 10, 20 )
      SIZE 3347 (7.171:1 COMPRESSION)
   }
   FILTERS {
      USER_DEFINED_FILTER {
         FILTER_ID 32001
         COMMENT blosc
         PARAMS { 2 2 4 8000 4 1 0 }
      }
   }
   FILLVALUE {
      FILL_TIME H5D_FILL_TIME_IFSET
      VALUE  H5D_FILL_VALUE_DEFAULT
   }
   ALLOCATION_TIME {
      H5D_ALLOC_TIME_INCR
   }
   DATA {h5dump error: unable to print data

   }
}
}'
)

## ----h5dump-2, eval = FALSE---------------------------------------------------
# ## set environment variable to hdf5filter location
# Sys.setenv("HDF5_PLUGIN_PATH" = rhdf5filters::hdf5_plugin_path())
# h5dump_out <- system2('h5dump',
#                       args = c('-p', '-d /dset', '-w 50', blosc_file),
#                       stdout = TRUE,  stderr = TRUE)
# 
# ## find the data entry and print the first few lines
# DATA_line <- grep(h5dump_out, pattern = "DATA \\{")
# cat( h5dump_out[ (DATA_line):(DATA_line+2) ], sep = "\n" )

## ----h5dump-2-out, eval = TRUE, echo = FALSE----------------------------------
cat(
'  DATA {
   (0,0,0): 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
   (0,0,11): 11, 12, 13, 14, 15, 16, 17, 18,'
)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

