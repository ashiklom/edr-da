library(coda)
library(here)

samples <- readRDS(here("ed-outputs/multi_site_pda/progress.rds"))

coda_samps <- BayesianTools::getSample(samples, coda = TRUE)
