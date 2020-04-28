library(drake)
library(magrittr)
library(fs)
library(ggplot2)
library(here)
library(readr)
library(dplyr)

stopifnot(
  requireNamespace("fst", quietly = TRUE),
  requireNamespace("data.table", quietly = TRUE)
)

pkgload::load_all(".", attach_testthat = FALSE)
