#!/usr/bin/env Rscript

import::from("magrittr", "%>%")

load("priors/mvtraits_priors.RData")
pfts <- rownames(means)
pfts <- pfts[!grepl("Southern_Pine", pfts)]

# Fix row and column names of means and covariances
mv_priors <- c(
  "prospect_N", "prospect_Cab", "prospect_Car",
  "prospect_Cw", "prospect_Cm", "SLA"
)
colnames(means) <- rownames(covars) <- colnames(covars) <- mv_priors

prospect_means <- means[pfts, ]
prospect_covars <- covars[, , pfts]
prospect_names <- colnames(prospect_means)

# Load allometry priors
allometry_stats <- readRDS("priors/allometry_stats.rds") %>%
  purrr::map(list(18, "statistics")) %>%
  .[pfts]

allom_names <- c("b1Bl", "b2Bl")

allom_mu <- purrr::map(
  allometry_stats,
  ~.[c("mu0", "mu1"), "Mean"] %>% setNames(allom_names)
)

allom_Sigma <- purrr::map(
  allometry_stats,
  ~matrix(.[c("tau11", "tau12", "tau12", "tau22"), "Mean"], 2, 2)
)

b1Bl_means <- purrr::map_dbl(allom_mu, "b1Bl")
b1Bl_sds <- purrr::map_dbl(allom_Sigma, ~.[1, 1]) %>% sqrt()
npfts <- length(b1Bl_means)
pfts <- names(b1Bl_means)


usethis::use_data(prospect_means, prospect_covars, prospect_names,
                  allometry_stats, allom_mu, allom_Sigma, allom_names,
                  b1Bl_means, b1Bl_sds,
                  npfts, pfts,
                  internal = TRUE, overwrite = TRUE)
