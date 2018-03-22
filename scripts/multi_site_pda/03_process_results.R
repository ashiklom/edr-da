library(BayesianTools)
library(here)
library(tidyverse)

source(here("scripts/multi_site_pda/summarize_results.R"))
source(here("scripts/multi_site_pda/priors.R"))

#pda_dir <- here("ed-outputs", "multi_site_pda_allom")
#source("scripts/multi_site_pda/setup_sites.R")

message("Drawing priors")
prior_draws <- map(1:1000, ~prior$sampler()) %>% invoke(rbind, .)

message("Processing allom")
results_allom <- readRDS(here("pda_results/multi_site_pda_allom_2018-03-14.rds"))
summarize_results(results_allom, "multi_site_pda_allom")

message("Processing noallom")
results_noallom <- readRDS(here("pda_results/multi_site_pda_noallom_2018-03-14.rds"))
summarize_results(results_noallom, "multi_site_pda_noallom")
