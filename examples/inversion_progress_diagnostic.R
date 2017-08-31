#!/usr/bin/env Rscript
library(BayesianTools)

samples <- readRDS('.edr_inversion/inversion_samples_progress.rds')

pdf('progress_samples.pdf')
plot(samples)
dev.off()
