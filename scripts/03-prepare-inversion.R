library(conflicted)
conflict_prefer("filter", "dplyr")

library(here)
library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(fs)
library(future)
library(furrr)

pkgload::load_all()

# Translate from ED PFTs to my PFTs
# Use as `match(pft, pft_dict)`
pft_dict <- c(9, 10, 11, 6, 8)
npft <- length(pft_dict)

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

make_likelihood <- function(site_data, refl, site, czen, direct_sky_frac) {
  dbh <- site_data[["dbh"]]
  nplant <- site_data[["n"]]
  pft_orig <- site_data[["pft"]]
  pft_dict <- c(9, 10, 11, 6, 8)
  pft <- match(pft_orig, pft_dict)

  # Calculate heights and height order (shortest first)
  hite <- dbh2h(dbh, pft)
  ihite <- order(hite, decreasing = FALSE)

  # Order cohorts by decreasing height (shortest first)
  dbh <- dbh[ihite]
  pft <- pft[ihite]
  nplant <- nplant[ihite]
  hite <- hite[ihite]

  ncohort <- length(dbh)

  # Return a likelihood function
  function(params) {
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

    # Incorporate LAI values in the likelihood.
    ll <- ll +
      dlnorm(sum(lai), 1, 0.5, log = TRUE)
    if (!is.finite(ll)) return(-Inf)

    # Cohort area index is constant (no crown radius model)
    cai <- rep(1, ncohort)

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
}

closure_data <- tibble(site = names(site_data), site_data = site_data) %>%
  left_join(aviris_data_sub, c("site" = "iPLOT")) %>%
  mutate(likelihood = pmap(
    list(site_data, refl, site, czen, direct_sky_frac),
    make_likelihood
  ))

closures <- closure_data %>%
  select(site, likelihood) %>%
  deframe()

overall_likelihood <- function(params) {
  site_lls <- furrr::future_map_dbl(
    sample(closures),
    ~purrr::possibly(.x, -Inf, quiet = FALSE)(params)
  )
  sum(site_lls)
}

site_names <- unique(names(closures))
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

# Test likelihood evaluation
psamps <- check_prior(prior, error = TRUE)
for (i in 1:5) {
  message("Testing ", i)
  print(overall_likelihood(psamps[i,]))
}

inversion_data_dir <- dir_create(here("data-raw", "inversion-data"))

save.image(file = path(inversion_data_dir, "inversion-data.RData"))
