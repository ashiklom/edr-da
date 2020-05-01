library(drake)
library(magrittr)
library(fs)
library(ggplot2)
library(here)
library(readr)
library(dplyr)
library(patchwork)

stopifnot(
  requireNamespace("fst", quietly = TRUE),
  requireNamespace("tidyr", quietly = TRUE),
  requireNamespace("forcats", quietly = TRUE)
)

pkgload::load_all(".", attach_testthat = FALSE)
drake::expose_imports("redr")
