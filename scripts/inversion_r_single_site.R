## library(redr)
devtools::load_all(".")

load_local <- function(file) {
  menv <- new.env(parent = baseenv())
  load(file, envir = menv)
  as.list(menv)
}

site_code <- "BH02"

# Read AVIRIS observation
aviris_data <- readRDS("aviris/aviris.rds")
aviris_all_wl <- round(PEcAnRTM::wavelengths(aviris_data))
aviris_use_wl_p <- aviris_all_wl >= 400 & aviris_all_wl <= 1300
aviris_use_wl <- aviris_all_wl[aviris_use_wl_p]
observation <- aviris_data[aviris_use_wl_p, colnames(aviris_data) == site_code]
which_wl <- aviris_use_wl - 399

# Read site file
site_dir <- list.files("sites", site_code, full.names = TRUE)
site_file <- tail(list.files(site_dir, "css$", full.names = TRUE), 1)
site_dat <- PEcAn.ED2::read_css(site_file)

dbh <- site_dat[["dbh"]]
pft <- site_dat[["pft"]] - 7
nplant <- site_dat[["n"]]
npft <- max(pft)

# Define likelihood
create_likelihood <- function(observed, waves, pft, dbh, nplant) {
  ncohort <- length(dbh)
  npft <- max(pft)
  stopifnot(
    length(pft) == ncohort,
    length(dbh) == ncohort,
    length(nplant) == ncohort
  )
  function(params) {
    # Define constants
    ## direct_sky_frac <- 0.9
    ## czen <- 1
    # Pull out site-specific params
    ssigma <- params[1]
    soil_moisture <- params[2]
    direct_sky_frac <- params[3]
    czen <- params[4]

    # Remaining params are pft-specific
    # Construct a parameter matrix -- nparam x npft
    # Each column contains the parameters in the order:
    # b1leaf, b2leaf, sla,
    # b1wood, b2wood,
    # N, Cab, Car, Cw, Cm,
    # clumping_factor, orient_factor
    pft_params <- matrix(params[-(1:4)], ncol = npft)

    # Calculate allometries
    b1leaf <- pft_params[1, pft]
    b2leaf <- pft_params[2, pft]
    sla <- pft_params[3, pft]
    bleaf <- size2bl(dbh, b1leaf, b2leaf)
    lai <- nplant * bleaf * sla

    b1wood <- pft_params[4, pft]
    b2wood <- pft_params[5, pft]
    wai <- wai_allometry(dbh, nplant, b1wood, b2wood)
    cai <- rep(1, ncohort)

    # Extract remaining parameters
    N <- pft_params[6, ]
    Cab <- pft_params[7, ]
    Car <- pft_params[8, ]
    Cw <- pft_params[9, ]
    Cm <- pft_params[10, ]
    orient_factor <- pft_params[11, ]
    clumping_factor <- pft_params[12, ]

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
    if (is.null(result)) return(-1e20)
    albedo <- result[["albedo"]]
    if (any(!is.finite(albedo))) return(-1e20)
    if (any(albedo < 0)) return(-1e20)
    if (any(albedo > 1)) return(-1e20)

    # Calculate likelihood
    ll <- sum(dnorm(albedo, observed, ssigma, log = TRUE))
    ll
  }
}

# Define priors
## prior_mvtraits <- load_local("priors/mvtraits_priors.RData")

pft_lowers <- rep(c(
  b1leaf = 0, b2leaf = -5, sla = 0,
  b1wood = 0, b2wood = -5,
  N = 1, Cab = 0, Car = 0, Cw = 0, Cm = 0,
  orient_factor = -0.75, clumping_factor = 0
), npft)
pft_uppers <- rep(c(
  b1leaf = 3, b2leaf = -2, sla = 100,
  b1wood = 3, b2wood = -2,
  N = 3, Cab = 100, Car = 40, Cw = 0.05, Cm = 0.05,
  orient_factor = 0.75, clumping_factor = 1
), npft)
lowers <- c(
  "ssigma" = 0, "soil_moisture" = 0, "direct_sky_frac" = 0, "czen" = 0,
  pft_lowers
)
uppers <- c(
  "ssigma" = 1, "soil_moisture" = 1, "direct_sky_frac" = 1, "czen" = 1,
  pft_uppers
)
prior <- BayesianTools::createUniformPrior(lower = lowers, upper = uppers)

likelihood <- create_likelihood(observation, aviris_use_wl, pft, dbh, nplant)

## debug(edr_r)
## likelihood(prior$sampler())

# Run inversion
setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = TRUE)
samples <- BayesianTools::runMCMC(setup)
samples <- BayesianTools::runMCMC(samples)
BayesianTools::gelmanDiagnostics(samples)
summary(samples)

## Local Variables:
## ess-r-package--project-cache: (redr . "/Users/shik544/Box Sync/Projects/edr_pda/edr-da/")
## End:
