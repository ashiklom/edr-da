library(BayesianTools)
library(here)
library(tidyverse)

source(here("scripts/multi_site_pda/summarize_results.R"))
source(here("scripts/multi_site_pda/priors.R"))

#pda_dir <- here("ed-outputs", "multi_site_pda_allom")
#source("scripts/multi_site_pda/setup_sites.R")

message("Drawing priors")
prior_draws <- map(1:1000, ~prior$sampler()) %>% invoke(rbind, .)
message("Done drawing priors")

prior_summary <- tibble(
  params = colnames(prior_draws),
  `Mean` = colMeans(prior_draws),
  `2.5%` = apply(prior_draws, 2, quantile, 0.025),
  `97.5%` = apply(prior_draws, 2, quantile, 0.975),
  type = "prior"
) %>%
  split_params("params")

message("Processing results")
samples <- readRDS(here("ed-outputs/multi_site_pda/progress.rds"))
prefix <- "multi_site_pda"
burnin <- 10000

summarize_results(samples, prefix, burnin)
