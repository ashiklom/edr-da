library(BayesianTools)
library(here)

samples <- readRDS(here("ed-outputs/multi_site_pda_allom/progress.rds"))
gelmanDiagnostics(samples)

coda_samps <- getSample(samples, coda = TRUE)
