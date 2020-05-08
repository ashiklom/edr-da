library(magrittr)
source("drake/functions.R")

for (conf in c("homo-pooled", "hetero-pooled",
               "homo-sitespecific", "hetero-sitespecific")) {
  message("Plotting ", conf)
  f <- last_result_file(file.path("multi_site_pda_results", conf))
  s <- BayesianTools::getSample(readRDS(f), coda = TRUE, thin = "auto")
  pdf(fs::path(fs::path_dir(f), "trace.pdf"))
  plot(s)
  dev.off()
}
