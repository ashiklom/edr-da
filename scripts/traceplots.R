library(magrittr)
source("drake/functions.R")

for (conf in c("homo-pooled", "hetero-pooled",
               "homo-sitespecific", "hetero-sitespecific")) {
  message("Plotting ", conf)
  f <- last_result_file(file.path("multi_site_pda_results", conf))
  param_names <- readLines(file.path(dirname(f), "param_names.txt"))
  s <- BayesianTools::getSample(readRDS(f), coda = TRUE, thin = "auto")
  coda::varnames(s) <- param_names
  pdf(fs::path(fs::path_dir(f), "trace.pdf"))
  plot(s)
  dev.off()
}
