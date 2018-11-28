library(redr)
requireNamespace("PEcAnRTM")
## devtools::load_all(".")

# Constants
npft <- 5       # Early, mid, late hardwood, early & late conifer

# HACK: These are not really constants -- need to be site-specific
direct_sky_frac <- 0.9                  # Relatively clear day
czen <- 1                               # Directly overhead

# Translate from ED PFTs to my PFTs
# Use as `match(pft, pft_dict)`
pft_dict <- c(9, 10, 11, 6, 8)

# Constant parameters (allometry slopes)
b2Bl <- purrr::map_dbl(allom_mu, "b2Bl")
b2Bw <- purrr::map_dbl(wallom_mu, "b2Bw")

# Read AVIRIS data
sites <- readLines(here::here("other_site_data", "site_list"))
nsite <- length(sites)
observed <- load_observations(sites)
waves <- PEcAnRTM::wavelengths(observed)
aviris_inds <- match(colnames(observed), sites)

# Read site data
site_dirs <- file.path("sites", sites)
site_files <- purrr::map(site_dirs, list.files, pattern = "css$", full.names = TRUE) %>%
  purrr::map_chr(tail, n = 1)
site_data_list <- purrr::map(site_files, PEcAn.ED2::read_css)
names(site_data_list) <- sites

# Set up prior
prior <- create_prior(nsite = nsite, heteroskedastic = FALSE, limits = TRUE)
psamps <- check_prior(prior, error = TRUE)

# Define likelihood
likelihood <- function(params) {
  ll <- 0
  for (i in seq_len(nsite)) { # site loop
    ## i <- 1
    site <- sites[i]
    site_data <- site_data_list[[site]]
    site_obs <- observed[, colnames(observed) == site]
    dbh <- site_data[["dbh"]]
    pft_orig <- site_data[["pft"]]
    pft <- match(pft_orig, pft_dict)
    if (any(is.na(pft))) {
      stop("Problem with PFTs in site ", site)
    }
    nplant <- site_data[["n"]]
    ncohort <- length(dbh)

    # Extract site-specific soil moisture
    soil_moisture <- params[grepl(paste0("sitesoil_", i, "$"), names(params))]

    # Extract residuals
    rss <- params["residual"]
    ## rs <- params["residual_slope"]
    ## ri <- params["residual_intercept"]

    # Remaining parameters are PFT-specific
    pft_params_v <- params[!grepl("residual|sitesoil", names(params))]
    # Create a matrix nparam (rows) x npft (cols)
    pft_params <- matrix(pft_params_v, ncol = npft)

    # Extract parameters
    N <- pft_params[1, ]
    Cab <- pft_params[2, ]
    Car <- pft_params[3, ]
    Cw <- pft_params[4, ]
    Cm <- pft_params[5, ]
    SLA <- pft_params[6, ]
    b1Bl <- pft_params[7, ]
    b1Bw <- pft_params[8, ]
    orient_factor <- pft_params[9, ]
    clumping_factor <- pft_params[10, ]

    # Calculate allometries
    bleaf <- size2bl(dbh, b1Bl[pft], b2Bl[pft])
    lai <- nplant * bleaf * SLA[pft]
    wai <- wai_allometry(dbh, nplant, b1Bw[pft], b2Bw[pft])

    # Cohort area index is constant (no crown radius model)
    cai <- rep(1, ncohort)

    # Call RTM
    result <- tryCatch(
      edr_r(pft, lai, wai, cai,
            N, Cab, Car, Cw, Cm,
            orient_factor, clumping_factor,
            soil_moisture,
            direct_sky_frac,
            czen,
            wavelengths = waves),
      error = function(e) NULL)
    if (is.null(result)) {
      message("Failed for site: ", site)
      return(-1e20)
    }
    albedo <- result[["albedo"]]
    if (any(!is.finite(albedo)) || any(albedo < 0) || any(albedo > 1)) {
      message("Bad albedo for site: ", site)
      return(-1e20)
    }
    site_ll <- sum(dnorm(albedo, site_obs, rss, log = TRUE))
    if (!is.finite(site_ll)) {
      message("Likelihood calculation failed for site ", site)
      return(-1e20)
    }
    ll <- ll + site_ll
  } # end site loop
  message("success")
  ll
}

## param <- readRDS("good_param.rds")

## likelihood(param)

## profvis::profvis({
##   for (i in 1:5) {
##     l <- likelihood(param)
##   }
## })

# Run inversion
setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
samples <- BayesianTools::runMCMC(setup)
samples <- BayesianTools::runMCMC(samples)
BayesianTools::gelmanDiagnostics(samples)
summary(samples)
