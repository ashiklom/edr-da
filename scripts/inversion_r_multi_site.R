#!/usr/bin/env Rscript
## utils::install.packages(".", repos = NULL, type = "source")
stopifnot(
  requireNamespace("PEcAnRTM", quietly = TRUE),
  requireNamespace("redr", quietly = TRUE)
)
pkgload::load_all(".", attach_testthat = FALSE)

arg <- commandArgs(trailingOnly = TRUE)
resume <- isTRUE(arg[1] == "resume")
message("Value of resume is: ", resume)

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
site_data_list <- purrr::map(site_files, read.table, header = TRUE)
names(site_data_list) <- sites

# Set up prior
prior <- create_prior(nsite = nsite, heteroskedastic = FALSE, limits = TRUE)
psamps <- prior$sampler()
param_names <- names(psamps)
# Re-create prior, but with parameter names
message("Creating prior...")
prior <- create_prior(nsite = nsite, heteroskedastic = FALSE, limits = TRUE, param_names = param_names)
message("Testing prior...")
psamps <- check_prior(prior, error = TRUE)

# Define likelihood
likelihood <- function(params) {
  ## params <- psamps[1,]
  ll <- 0
  npft_param <- 10  # Number of PFT-specific parameters
  has_names <- !is.null(names(params))
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

    # Calculate heights and height order (tallest first)
    hite <- dbh2h(dbh, pft)
    ihite <- order(hite, decreasing = TRUE)

    # Order cohorts by decreasing height (tallest first)
    dbh <- dbh[ihite]
    pft <- pft[ihite]
    nplant <- nplant[ihite]
    hite <- hite[ihite]

    # Extract site-specific soil moisture
    soil_moisture <- if (has_names) {
      params[grepl(paste0("sitesoil_", i, "$"), names(params))] 
    } else {
      params[npft * npft_param + i]
    }

    # Extract residuals
    rss <- if (has_names) params["residual"] else tail(params, 1)
    ## rs <- params["residual_slope"]
    ## ri <- params["residual_intercept"]

    # Remaining parameters are PFT-specific
    pft_params_v <- if (has_names) {
      params[!grepl("residual|sitesoil", names(params))]
    } else {
      head(params, -(nsite + 1))
    }
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
      ## cat("-")
      ## message("Failed for site: ", site)
      return(-1e20)
    }
    albedo <- result[["albedo"]]
    if (any(!is.finite(albedo)) || any(albedo < 0) || any(albedo > 1)) {
      ## cat("-")
      ## message("Bad albedo for site: ", site)
      return(-1e20)
    }
    site_ll <- sum(dnorm(albedo, site_obs, rss, log = TRUE))
    if (!is.finite(site_ll)) {
      ## cat("-")
      ## message("Likelihood calculation failed for site ", site)
      return(-1e20)
    }
    ll <- ll + site_ll
  } # end site loop
  ## cat("x")
  ## message("success")
  ll
}

if (FALSE) {
  param <- readRDS("good_param.rds")
  likelihood(param)
  profvis::profvis({
    for (i in 1:5) {
      l <- likelihood(param)
    }
  })
}

# Create directory for storage
sampdir <- strftime(Sys.time(), "%Y-%m-%d-%H%M")
base_outdir <- "multi_site_pda_results"
outdir <- file.path(base_outdir, sampdir)
dir.create(outdir, recursive = TRUE)
message("Storing results in: ", outdir)

# Inversion settings
niter <- 500
max_iter <- 5e6
max_attempts <- floor(max_iter / niter)
attempt <- 0
threshold <- 1.2
target_neff <- 500
ncores <- 8

# Create BayesianTools setup
message("Creating setup")
samples <- BayesianTools::createBayesianSetup(
  likelihood,
  prior,
  parallel = ncores
)
settings <- list(
  iterations = niter,
  consoleUpdates = 10,
  startValue = prior$sampler(ncores)
)

if (resume) {
  message("Resuming from previous sampling")
  last_samplefile <- tail(list.files(base_outdir, ".rds",
                          recursive = TRUE, full.names = TRUE), 1)
  message("Resuming from sample file: ", last_samplefile)
  stopifnot(length(last_samplefile) > 0,
            file.exists(last_samplefile))
  samples <- readRDS(last_samplefile)
} else {
  message("Starting fresh inversion")
}

repeat {
  attempt <- attempt + 1
  message("Sampling attempt: ", attempt)
  samples <- BayesianTools::runMCMC(samples, settings = settings)
  saveRDS(samples, file.path(outdir, "current_samples.rds"))
  nsamp <- nrow(samples$chain[[1]])
  coda_samples <- BayesianTools::getSample(
    samples,
    start = if (nsamp > 5000) 5000 else floor(nsamp / 2),
    thin = "auto",
    coda = TRUE
  )
  gd <- coda::gelman.diag(
    coda_samples,
    multivariate = FALSE,
    autoburnin = FALSE
  )
  psrf <- gd[["psrf"]][, 1]
  names(psrf) <- param_names
  exceeds <- psrf > threshold
  if (any(exceeds)) {
    message("The following parameters have not converged:")
    print(psrf[exceeds])
  } else {
    neff <- coda::effectiveSize(coda_samples)
    too_few <- neff < target_neff
    if (!any(too_few)) {
      message("Converged!")
      break
    }
    message("Converged, but the following parameters have too few samples:")
    print(neff[too_few])
  }

  message("Resuming sampling...")
  if (attempt > max_attempts) {
    message("Failed to converge after attempt: ", attempt)
    break
  }
}
