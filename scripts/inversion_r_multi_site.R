#!/usr/bin/env Rscript
## utils::install.packages(".", repos = NULL, type = "source")
## commandArgs <- function(...) c("ss")
stopifnot(
  requireNamespace("PEcAnRTM", quietly = TRUE),
  requireNamespace("redr", quietly = TRUE)
)
pkgload::load_all(".", attach_testthat = FALSE)

argv <- commandArgs(trailingOnly = TRUE)

stopifnot(all(argv %in% c("resume", "hetero", "ss")))
resume <- "resume" %in% argv
HETEROSKEDASTIC <- "hetero" %in% argv
SITE_SPECIFIC <- "ss" %in% argv
NCORES <- 8

message("Value of resume is: ", resume)
message("Value of heteroskedastic is: ", HETEROSKEDASTIC)
message("Value of site_specific is: ", SITE_SPECIFIC)

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
prior <- create_prior(nsite = nsite, heteroskedastic = HETEROSKEDASTIC,
                      limits = TRUE, site_specific_var = SITE_SPECIFIC)
psamps <- prior$sampler()
param_names <- names(psamps)
# Re-create prior, but with parameter names
message("Creating prior...")
prior <- create_prior(nsite = nsite, heteroskedastic = HETEROSKEDASTIC,
                      limits = TRUE, param_names = param_names,
                      site_specific_var = SITE_SPECIFIC)
message("Testing prior...")
psamps <- check_prior(prior, error = TRUE)

# Define likelihood
likelihood <- function(params) {
  ll <- 0
  npft_param <- 10  # Number of PFT-specific parameters
  has_names <- !is.null(names(params))
  # NOTE: Sample here to improve efficiency by quickly rejecting site where
  # params fail
  for (i in sample(nsite)) { # site loop
    site <- sites[i]
    site_data <- site_data_list[[site]]
    site_obs <- observed[, colnames(observed) == site]
    nobs <- NCOL(site_obs)
    nwl <- NROW(site_obs)

    dbh <- site_data[["dbh"]]
    pft_orig <- site_data[["pft"]]
    pft <- match(pft_orig, pft_dict)
    if (any(is.na(pft))) {
      stop("Problem with PFTs in site ", site)
    }
    nplant <- site_data[["n"]]
    ncohort <- length(dbh)

    # Calculate heights and height order (shorest first)
    hite <- dbh2h(dbh, pft)
    ihite <- order(hite, decreasing = FALSE)

    # Order cohorts by decreasing height (shortest first)
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

    # Remaining parameters are PFT-specific
    pft_params_v <- params[!grepl("residual|sitesoil", param_names)]
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
    clumping_factor <- pft_params[9, ]
    orient_factor <- pft_params[10, ]

    # Calculate allometries
    bleaf <- size2bl(dbh, b1Bl[pft], b2Bl[pft])
    lai <- nplant * bleaf * SLA[pft]
    wai <- wai_allometry(dbh, nplant, b1Bw[pft], b2Bw[pft])

    # Cohort area index is constant (no crown radius model)
    cai <- rep(1, ncohort)

    # Call RTM
    result <- edr_r(pft, lai, wai, cai,
                    N, Cab, Car, Cw, Cm,
                    orient_factor, clumping_factor,
                    soil_moisture,
                    direct_sky_frac,
                    czen,
                    wavelengths = waves)
    albedo <- result[["albedo"]]
    stopifnot(
      all(is.finite(albedo)),
      all(albedo >= 0),
      all(albedo <= 1)
    )
    # Extract residuals
    if (HETEROSKEDASTIC) {
      rs <- params[grep(
        paste0("residual_slope", if (SITE_SPECIFIC) i, "$"),
        param_names
      )]
      ri <- params[grep(
        paste0("residual_intercept", if (SITE_SPECIFIC) i, "$"),
        param_names
      )]
      rss <- ri + rs * albedo
    } else {
      rss <- params[grep(
        paste0("residual", if (SITE_SPECIFIC) i, "$"),
        param_names
      )]
    }

    site_ll <- sum(dnorm(albedo, site_obs, rss, log = TRUE))
    # NOTE: Down-weight sites with multiple observations
    site_ll <- site_ll - log(nobs) * log(nwl)
    stopifnot(is.finite(site_ll))
    ll <- ll + site_ll
  } # end site loop
  ll
}

# Test likelihood evaluation
for (i in 1:5) {
  message("Testing ", i)
  print(likelihood(psamps[i,]))
}

# Create directory for storage
sampdir <- strftime(Sys.time(), "%Y-%m-%d-%H%M")
outtag <- paste(
  if (HETEROSKEDASTIC) "hetero" else "homo",
  if (SITE_SPECIFIC) "sitespecific" else "pooled",
  sep = "-"
)
base_outdir <- file.path("multi_site_pda_results", outtag)
outdir <- file.path(base_outdir, sampdir)
dir.create(outdir, recursive = TRUE)
message("Storing results in: ", outdir)
writeLines(param_names, file.path(outdir, "param_names.txt"))

# Inversion settings
niter <- 5000
max_iter <- 1e7
max_attempts <- floor(max_iter / niter)
attempt <- 0
threshold <- 1.2
target_neff <- 500
ncores <- NCORES

# Create BayesianTools setup
message("Creating setup")
newsamples <- BayesianTools::createBayesianSetup(
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
  # Need this to reset parallelism
  samples$setup$likelihood <- newsamples$likelihood
  samples$setup$posterior <- newsamples$posterior
} else {
  message("Starting fresh inversion")
  samples <- newsamples
}

repeat {
  attempt <- attempt + 1
  message("Sampling attempt: ", attempt)
  samples <- BayesianTools::runMCMC(samples, settings = settings, sampler = "DREAMzs")
  saveRDS(samples, file.path(outdir, "current_samples.rds"))
  nsamp <- nrow(samples$chain[[1]])
  coda_samples <- BayesianTools::getSample(
    samples,
    start = if (nsamp > 5000) 5000 else floor(nsamp / 2),
    thin = "auto",
    coda = TRUE
  )
  gd <- tryCatch(
    coda::gelman.diag(
      coda_samples,
      multivariate = TRUE,
      autoburnin = FALSE
    ),
    error = function (e) NULL)
  if (is.null(gd)) {
    message("GD calc failed. No convergence.")
  } else {
    psrf <- c(gd[["psrf"]][, 1], gd[["mpsrf"]])
    names(psrf) <- c(param_names, "mpsrf")
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
  }

  message("Resuming sampling...")
  if (attempt > max_attempts) {
    message("Failed to converge after attempt: ", attempt)
    break
  }
}
