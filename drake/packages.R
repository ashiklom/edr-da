library(drake)
library(magrittr)
library(fs)
library(ggplot2)
library(here)
library(readr)
library(dplyr)

stopifnot(
  requireNamespace("fst", quietly = TRUE)
)

pkgload::load_all(".", attach_testthat = FALSE)
drake::expose_imports("redr")
