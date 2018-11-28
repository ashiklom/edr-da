#!/usr/bin/env Rscript

## Priors
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

# Wood allometries
wood_allometry_stats <- readRDS("priors/wood_allometry_stats.rds") %>%
  purrr::map(list(16, "statistics")) %>%
  .[pfts]

wallom_names <- c("b1Bw", "b2Bw")

wallom_mu <- purrr::map(
  wood_allometry_stats,
  ~.[c("mu0", "mu1"), "Mean"] %>% setNames(wallom_names)
)

wallom_Sigma <- purrr::map(
  wood_allometry_stats,
  ~matrix(.[c("tau11", "tau12", "tau12", "tau22"), "Mean"], 2, 2)
)

# Wood allometry uncertainties are WAY too large due to uninformative
# prior. Reduce them to something reasonable.
default_wallom_sigma <- diag(3, 2)
wallom_Sigma[["temperate.Northern_Pine"]] <- default_wallom_sigma
wallom_Sigma[["temperate.Late_Conifer"]] <- default_wallom_sigma

b1Bw_means <- purrr::map_dbl(wallom_mu, "b1Bw")
b1Bw_sds <- purrr::map_dbl(wallom_Sigma, ~.[1, 1]) %>% sqrt()

## Wet and dry soil spectra

pecan_soil_file <- "~/Projects/pecan_project/pecan/modules/rtm/src/RTM/modules/dataSpec/dataSpec_soil.f90"
soil_raw <- readLines(pecan_soil_file)
start_dry_soil <- grep("^[[:space:]]+DATA \\(Rsoil1", soil_raw)[1]
start_wet_soil <- grep("^[[:space:]]+DATA \\(Rsoil2", soil_raw)[1]

raw_dry_soil <- soil_raw[start_dry_soil:(start_wet_soil - 1)]
raw_dry_soil_num <- grep("^[[:space:]]+[[:digit:]]", raw_dry_soil, value = TRUE)
raw_dry_soil_num <- gsub("[&/]", "", raw_dry_soil_num)
dry_soil <- as.numeric(unlist(strsplit(raw_dry_soil_num, ",")))

raw_wet_soil <- soil_raw[-(1:start_wet_soil)]
raw_wet_soil_num <- grep("^[[:space:]]+[[:digit:]]", raw_wet_soil, value = TRUE)
raw_wet_soil_num <- gsub("[&/]", "", raw_wet_soil_num)
wet_soil <- as.numeric(unlist(strsplit(raw_wet_soil_num, ",")))

### Wood spectra
wood_spec_file <- system.file("extdata/wood_reflect_par.dat", package = "PEcAnRTM")
wood_spec <- scan(wood_spec_file)
wood_spec <- c(wood_spec, tail(wood_spec, 1))

stopifnot(
  length(wet_soil) == 2101,
  length(dry_soil) == 2101,
  length(wood_spec) == 2101
)

usethis::use_data(prospect_means, prospect_covars, prospect_names,
                  allometry_stats, allom_mu, allom_Sigma, allom_names,
                  b1Bl_means, b1Bl_sds,
                  wood_allometry_stats, wallom_mu, wallom_Sigma, wallom_names,
                  b1Bw_means, b1Bw_sds,
                  npfts, pfts,
                  dry_soil, wet_soil, wood_spec,
                  internal = TRUE, overwrite = TRUE)
