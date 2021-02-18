library(conflicted)

conflict_prefer("filter", "dplyr")

library(drake)
library(magrittr)
library(fs)
library(ggplot2)
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)
library(forcats)
library(stringr)

stopifnot(
  requireNamespace("fst", quietly = TRUE),
  requireNamespace("ggmap", quietly = TRUE),
  requireNamespace("scales", quietly = TRUE)
)

pkgload::load_all(".", attach_testthat = FALSE)
## drake::expose_imports("redr")
