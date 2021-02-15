library(conflicted)
conflict_prefer("filter", "dplyr")

library(here)
library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(fs)
library(future)
library(foreach)
library(doFuture)

pkgload::load_all()

# Translate from ED PFTs to my PFTs
# Use as `match(pft, pft_dict)`
pft_dict <- c(9, 10, 11, 6, 8)

b2Bl <- purrr::map_dbl(allom_mu, "b2Bl")
b2Bw <- purrr::map_dbl(wallom_mu, "b2Bw")

aviris_data <- read_csv(here("other_site_data", "aviris-processed.csv"))

# Read site data
sites <- readLines(here::here("other_site_data", "site_list"))
nsite <- length(sites)
site_dirs <- file.path("sites", sites)
site_files <- map(site_dirs, list.files, pattern = "css$", full.names = TRUE) %>%
  map_chr(tail, n = 1)
site_data_list <- map(site_files, read.table, header = TRUE)
names(site_data_list) <- sites
site_data <- site_data_list %>%
  map(~.x[, c("pft", "dbh", "n")]) %>%
  map(as_tibble)

# Subset AVIRIS wavelengths
aviris_waves_all <- read_csv(here("aviris", "NASA_FFT", "aviris_c_wavelength.csv")) %>%
  pull(wavelength)
drop_wl <- aviris_waves_all < 400 | aviris_waves_all > 1300
waves <- aviris_waves_all[!drop_wl]
drop_bands <- paste0("band_", seq_along(aviris_waves_all))[drop_wl]

aviris_data_sub <- aviris_data %>%
  filter(iPLOT %in% names(site_data)) %>%
  select(-!!drop_bands) %>%
  rowwise() %>%
  mutate(refl = list(c_across(starts_with("band_")) / 10000)) %>%
  ungroup() %>%
  select(-starts_with("band_"))

# Constant parameters (allometry slopes)
b2Bl <- purrr::map_dbl(allom_mu, "b2Bl")
b2Bw <- purrr::map_dbl(wallom_mu, "b2Bw")

# AR1 autocorrelation coefficient calculated from earlier simulations
rho <- 0.6995862
# Convert to spectral autocorrelation matrix
times <- seq_along(waves)
rho_H <- rho ^ abs(outer(times, times, "-"))

site_likelihood <- function(params, inputs) {
  npft <- 5

  site <- inputs[["site"]]
  refl <- inputs[["refl"]][[1]]
  czen <- inputs[["czen"]][[1]]
  direct_sky_frac <- inputs[["direct_sky_frac"]]
  pft <- inputs[["pft"]][[1]]
  dbh <- inputs[["dbh"]][[1]]
  nplant <- inputs[["nplant"]][[1]]
  hite <- inputs[["hite"]][[1]]

  rho_H <- inputs[["rho_H"]][[1]]
  b2Bl <- inputs[["b2Bl"]][[1]]
  b2Bw <- inputs[["b2Bw"]][[1]]
  waves <- inputs[["waves"]][[1]]

  # Arrange cohorts by height (shortest first)
  ihite <- order(hite, decreasing = FALSE)
  pft <- pft[ihite]
  dbh <- dbh[ihite]
  nplant <- nplant[ihite]

  ll <- 0
  param_names <- names(params)
  soil_moisture <- params[grepl(paste0("sitesoil_", site, "$"), names(params))]
  # Remaining parameters are PFT-specific
  pft_params_v <- params[!grepl("residual|sitesoil", param_names)]
  # Create a matrix nparam (rows) x npft (cols)
  pft_params <- matrix(pft_params_v, ncol = npft)
  N <- pft_params[1, ]
  Cab <- pft_params[2, ]
  Car <- pft_params[3, ]
  Cw <- pft_params[4, ]
  Cm <- pft_params[5, ]
  SLA <- pft_params[6, ]
  b1Bl <- pft_params[7, ]
  b1Bw <- pft_params[8, ]
  clumping_factor <- pft_params[9, ]
  orient_factor <- pft_params[10, ]

  # Calculate allometries
  bleaf <- size2bl(dbh, b1Bl[pft], b2Bl[pft])
  lai <- nplant * bleaf * SLA[pft]
  wai <- wai_allometry(dbh, nplant, b1Bw[pft], b2Bw[pft])

  # "Uniform" LAI --- 0-10 is equal probability; gt 10 is illegal
  if (sum(lai) > 10) {
    warning("Excessive LAI value: ", sum(lai))
    return(-Inf)
  }
  # Similar uniform WAI
  if (sum(wai) > 5) {
    warning("Excessive WAI value: ", sum(wai))
   return(-Inf)
  }

  # Incorporate LAI values in the likelihood.
  ## ll <- ll +
  ##   dlnorm(sum(lai), 1, 0.5, log = TRUE)
  ## if (!is.finite(ll)) return(-Inf)

  # Cohort area index is constant (no crown radius model)
  cai <- rep(1, length(lai))

  model <- tryCatch(
    edr_r(pft, lai, wai, cai,
          N, Cab, Car, Cw, Cm,
          orient_factor, clumping_factor,
          soil_moisture,
          direct_sky_frac,
          czen,
          wavelengths = waves),
    error = function(e) message(conditionMessage(e)))
  if (is.null(model)) {
    stop("Error in model execution")
  }
  albedo <- model[["albedo"]]
  if (!is.finite(albedo) || any(albedo > 1) || any(albedo < 0)) {
    stop("Model execution produced invalid albedo")
  }

  # Heteroskedastic residual model
  rs <- params[grep(paste0("residual_slope$"), param_names)]
  ri <- params[grep(paste0("residual_intercept$"), param_names)]
  rss <- diag(ri + rs * albedo)

  # AR1 covariance matrix
  covar <- rss %*% rho_H %*% rss
  site_ll <- mvtnorm::dmvnorm(albedo, refl, covar, log = TRUE, checkSymmetry = FALSE)
  ll <- ll + site_ll
  return(ll)
}

input_data <- tibble(site = names(site_data), site_data = site_data) %>%
  left_join(aviris_data_sub, c("site" = "iPLOT")) %>%
  mutate(
    dbh = map(site_data, "dbh"),
    pft = map(map(site_data, "pft"), match, c(9, 10, 11, 6, 8)),
    nplant = map(site_data, "n"),
    hite = map2(dbh, pft, dbh2h),
    rho_H = list(rho_H),
    b2Bl = list(b2Bl),
    b2Bw = list(b2Bw),
    waves = list(waves)
  )

overall_likelihood <- function(params) {
  site_lls <- foreach(s = sample(nrow(input_data)), .combine = c) %dopar% {
    inputs <- input_data[s,]
    ll <- site_likelihood(params, inputs)
    ll
  }
  sum(site_lls)
}

site_names <- unique(input_data$site)
nsite <- length(site_names)

# Set up prior
message("Creating prior...")
prior <- create_prior(nsite = nsite, heteroskedastic = TRUE,
                      limits = TRUE, site_specific_var = FALSE)
psamps <- prior$sampler()
param_names_raw <- names(psamps)

# Fix site-soil names so they use the actual site
param_names <- param_names_raw
for (i in seq_along(site_names)) {
  param_names <- gsub(sprintf("(?<=sitesoil_)%d$", i), site_names[i], param_names,
                      perl = TRUE)
}

# Re-create prior, but with parameter names
prior <- create_prior(nsite = nsite, heteroskedastic = TRUE,
                      limits = TRUE, param_names = param_names,
                      site_specific_var = FALSE)
message("Testing prior...")
psamps <- check_prior(prior, error = TRUE)

message("Running on ", availableCores(), " cores.")
registerDoFuture()
plan("multicore")
options(future.rng.onMisuse = "ignore")
for (i in 1:5) {
  message("Testing ", i)
  params <- psamps[i,]
  params[grep("b1Bl", names(params))] <- 0.001
  params[grep("b1Bw", names(params))] <- 0.001
  print(overall_likelihood(params))
}

inversion_data_dir <- dir_create(here("data-raw", "inversion-data"))

save(
  input_data, overall_likelihood, site_likelihood,
  prior, psamps,
  file = path(inversion_data_dir, "inversion-data.RData")
)
