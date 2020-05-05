library(magrittr)
source("drake/functions.R")

## f <- last_result_file("multi_site_pda_results")
f <- last_result_file("multi_site_pda_results-exp")
s <- BayesianTools::getSample(readRDS(f), coda = TRUE, thin = "auto")

pdf(fs::path(fs::path_dir(f), "trace.pdf"))
plot(s)
dev.off()
