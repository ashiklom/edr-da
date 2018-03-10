library(PEcAnRTM)
library(PEcAn.ED2)
library(here)
library(purrr)
library(foreach)
library(doParallel)

hostname <- system2("hostname", stdout = TRUE)

# Configuration for geo
img_path <- NULL
edr_exe_path <- "/projectnb/dietzelab/ashiklom/ED2/EDR/build/ed_2.1-dbg"

#if (hostname == "ashiklom") {
  #img_path <- "~/Projects/ED2/ed2.simg"
  #edr_exe_path <- NULL
#}

pda_dir <- here("ed-outputs", "multi_site_pda_allom")
stopifnot(file.exists(pda_dir))
all_sites <- list.files(pda_dir, "_site_")
stopifnot(length(all_sites) > 0)

############################################################
# Load observations
############################################################
message("Loading observations")
aviris_specfile <- here("aviris/aviris.h5")
aviris_h5 <- hdf5r::H5File$new(aviris_specfile)
site_tags <- gsub("([[:alnum:]]+)_.*", "\\1", all_sites)
av_inds <- map(site_tags, ~which(aviris_h5[["iPLOT"]][] == .))
# Subset sites to those that have observations
# Problem sites were causing trouble with the inversion
#problem_sites <- c("BH02_site_1-25665", "IDS05_site_1-25682")
problem_sites <- c("IDS05_site_1-25682")
keep_sites <- map_lgl(av_inds, ~length(.) > 0) &
  !(all_sites %in% problem_sites)

sites <- all_sites[keep_sites]

# Grab AVIRIS data
aviris_spec_raw <- map(av_inds[keep_sites], ~aviris_h5[["spectral_data"]][, .])
aviris_wl <- aviris_h5[["wavelengths"]][]
aviris_h5$close_all()
aviris_spec <- map(aviris_spec_raw, spectra, wavelengths = aviris_wl)
names(aviris_spec) <- sites

# Resample AVIRIS data to integer wavelengths, and keep only VIS and NIR
use_wl <- seq(400, 1300, by = 10)
aviris_spec <- map(aviris_spec, resample, to = use_wl)
aviris_nspec <- map_int(aviris_spec, ncol)
aviris_inds <- rep(seq_along(sites), aviris_nspec)
observed <- do.call(cbind, aviris_spec) / 10000
stopifnot(!any(is.na(observed)))

############################################################
# Setup EDR
############################################################
message("Setting up EDR")
ed2in_paths <- map(
  sites,
  ~list.files(file.path(pda_dir, .), "ED2IN", full.names = TRUE)
)

ed2in_sites <- map(ed2in_paths, read_ed2in)
edr_paths <- file.path(pda_dir, sites, "edr")
site_setup <- map2(ed2in_sites, edr_paths, setup_edr)

############################################################
# Define RTM inversion model
############################################################
# setup parallel
if (!exists("ncores")) {
  arg <- commandArgs(trailingOnly = TRUE)
  if (length(arg) == 0) {
    ncores <- 1
  } else {
    ncores <- as.numeric(arg)
  }
}

message("Creating cluster with ", ncores, " cores.")
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Returns a matrix of spectra
model <- function(params) {
  edr_in <- params2edr(params, prospect = TRUE, version = 5)
  #result <- matrix(0.0, 2101, length(site_setup))
  #for (i in seq_along(site_setup)) {
  pkgs <- c("PEcAnRTM", "PEcAn.ED2")
  exports <- c("img_path", "edr_in", "edr_exe_path")
  result <- tryCatch({
    foreach(s = site_setup, .combine = cbind, .packages = pkgs, .export = exports) %dopar% {
      EDR(
        img_path = img_path,
        ed2in_path = s,
        spectra_list = edr_in$spectra_list,
        trait.values = edr_in$trait.values,
        verbose_error = FALSE,
        edr_exe_path = edr_exe_path
      )
    }
  }, error = function(e) {
    message("EDR failed for at least one site. Returning -1e10 and continuing.")
    NULL
  })
  # Make result match the dimensions of the aviris spectra
  if (!is.null(result)) {
    colnames(result) <- sites
    result[use_wl - 399, aviris_inds]
  } else {
    -1e10
  }
}

############################################################
# Set up prior
############################################################
message("Setting up prior")
load(here("priors/mvtraits_priors.RData"))
pfts <- rownames(means)

# Fix row and column names of means and covariances
mv_priors <- c(
  "prospect_N", "prospect_Cab", "prospect_Car",
  "prospect_Cw", "prospect_Cm", "SLA"
)
colnames(means) <- rownames(covars) <- colnames(covars) <- mv_priors

# Draw only positive values from multivariate prior
rmvnorm_positive <- function(mu, Sigma) {
  draw <- -1
  while(any(draw < 0)) {
    draw <- mvtnorm::rmvnorm(1, mu, Sigma)
  }
  draw
}

# Load allometry priors
allom_stats <- readRDS(here("priors/allometry_stats.rds")) %>%
  map(list(18, "statistics"))

allom_names <- c("b1Bl", "b2Bl")
allom_mu <- map(allom_stats, ~.[c("mu0", "mu1"), "Mean"] %>% setNames(allom_names))
allom_Sigma <- map(allom_stats, ~solve(matrix(.[c("tau11", "tau12", "tau12", "tau22"), "Mean"], 2, 2)))

# Define priors for other parameters
prior_clumping <- c(0, 1)
rclumping <- function(n) runif(length(pfts), prior_clumping[1], prior_clumping[2])
dclumping <- function(x, log = TRUE) dunif(x, prior_clumping[1], prior_clumping[2], log = log)

prior_orient <- c(6, 4) * 3
rorient <- function(n) 2 * rbeta(n, prior_orient[1], prior_orient[2]) - 1
dorient <- function(x, log = TRUE) dbeta((x + 1) / 2, prior_orient[1], prior_orient[2], log = log)

prior_residual <- c(0.01, 0.01)

prior_sample <- function() {
  mv_draws <- -1
  mv_draws <- map(
    pfts,
    ~rmvnorm_positive(means[.,], covars[, , .]) %>% setNames(mv_priors)
  )
  allom_draws <- map2(
    allom_mu,
    allom_Sigma,
    ~mvtnorm::rmvnorm(1, .x, .y)[1, ]
  ) %>%
    map(~c(exp(.[1]), .[2]))
  all_draws <- pmap(
    list(
      mv_draws,
      allom_draws,
      clumping_factor = rclumping(length(pfts)),
      orient_factor = rorient(length(pfts))
    ),
    c
  )
  names(all_draws) <- pfts
  c(unlist(all_draws), residual = rgamma(1, prior_residual[1], prior_residual[2]))
}

prior_density <- function(params) {
  residual <- params["residual"]
  res_dens <- dgamma(residual, prior_residual[1], prior_residual[2], log = TRUE)
  if (!is.finite(res_dens)) {
    #message("Residual density is not finite.")
    #print(residual)
    return(-Inf)
  }
  traits <- params2edr(params, prospect = FALSE)$trait.values
  mvdens <- imap_dbl(
    traits,
    ~mvtnorm::dmvnorm(.x[mv_priors], means[.y, ], covars[, , .y], log = TRUE)
  )
  if (any(!is.finite(mvdens))) {
    #message("Multivariate density is not finite.")
    #print(params)
    return(-Inf)
  }
  allom_dens <- imap_dbl(
    traits,
    ~mvtnorm::dmvnorm(c(log(.x[allom_names[1]]), .x[allom_names[2]]),
                      allom_mu[[.y]], allom_Sigma[[.y]], log = TRUE)
  )
  if (any(!is.finite(allom_dens))) {
    return(-Inf)
  }
  cf_dens <- dclumping(map_dbl(traits, "clumping_factor"), log = TRUE)
  if (any(!is.finite(cf_dens))) {
    #message("Clumping factor density not finite.")
    #print(map_dbl(traits, "clumping_factor")[!is.finite(cf_dens)])
    return(-Inf)
  }
  of_dens <- dorient(map_dbl(traits, "orient_factor"), log = TRUE)
  if (any(!is.finite(of_dens))) {
    #message("Orient factor density not finite.")
    #print(map_dbl(traits, "orient_factor")[!is.finite(of_dens)])
    return(-Inf)
  }
  logdens <- sum(mvdens, allom_dens, cf_dens, of_dens, res_dens)
  if (!is.finite(logdens)) {
    #message("Log density is not finite.")
    return(-Inf)
  }
  logdens
}

prior <- BayesianTools::createPrior(
  density = prior_density,
  sampler = prior_sample
)

############################################################
# Test
############################################################
message("Testing prior sampler and density functions.")
for (i in 1:100) {
  test_params <- prior$sampler()
  test_priors <- prior$density(test_params)
  stopifnot(is.numeric(test_priors), is.finite(test_priors))
}
message("Priors seem reliable")

#debugonce(model)
#debugonce(EDR)
tm <- -1e10
i <- 0
imax <- 5
message("Testing EDR")
#debugonce(model)
while(length(tm) == 1){
  i <- i + 1
  if (i > imax) {
    stop("No successful EDR runs after several attempts. Something is probably wrong.")
  }
  message("Testing multi_site EDR. Attempt # ", i)
  tm <- model(prior$sampler())
}
message("At least one run of EDR was successful!")

############################################################
# Run PDA
############################################################
message("Initializing PDA")
PEcAn.logger::logger.setLevel("INFO")

custom_settings <- list(
  init = list(iterations = 50),
  loop = list(iterations = 50),
  other = list(
    save_progress = file.path(pda_dir, "progress.rds")
  )
)
samples <- invert_bt(observed, model, prior, custom_settings = custom_settings)
saveRDS(samples, file.path(pda_dir, "results.rds"))
