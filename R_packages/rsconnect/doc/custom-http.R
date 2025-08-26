## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(rsconnect)

## -----------------------------------------------------------------------------
options(
  rsconnect.http.headers = c(
    "CustomHeader1" = "CustomValue", "CustomHeader2" = "CustomValue2"
  )
)

## -----------------------------------------------------------------------------
options(
  rsconnect.http.headers = c("cookie1=value1", "cookie2=value2")
)

## -----------------------------------------------------------------------------
options(
  rsconnect.http.headers = "cookie1=value1; Expires=Thu, 31 Oct 2021 07:28:00 GMT; Secure"
)

