library(PEcAnRTM)
library(PEcAn.ED2)
library(here)
library(purrr)

hostname <- system2("hostname", stdout = TRUE)
if (hostname == "ashiklom") {
  img_path <- "~/Projects/ED2/ed2.simg"
  edr_exe_path <- NULL
} else if (hostname == "geo") {
  img_path <- NULL
  edr_exe_path <- "~/dietzelab/ED2/EDR/build/ed_2.1-dbg"
} else {
  stop("Set img_path and/or edr_exe_path manually")
}

pda_dir <- here("ed-outputs", "multi_site_pda")
all_sites <- list.files(pda_dir, "_site_")

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
problem_sites <- c("BH02_site_1-25665", "IDS05_site_1-25682")
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

# Returns a matrix of spectra
model <- function(params) {
  edr_in <- params2edr(params, prospect = TRUE, version = 5)
  result <- matrix(0.0, 2101, length(site_setup))
  colnames(result) <- sites
  for (i in seq_along(site_setup)) {
    spec <- tryCatch({
      EDR(
        img_path = img_path,
        ed2in_path = site_setup[[i]],
        spectra_list = edr_in$spectra_list,
        trait.values = edr_in$trait.values,
        verbose_error = FALSE,
        edr_exe_path = edr_exe_path
      )
    }, error = function(e) {
      message(
        "EDR failed for site ", i, sites[i],
        ". Moving to next iteration."
      )
    }
    )
    if (is.null(spec)) {
      return(-1e10)
    }
    result[, i] <- spec
  }
  # Make result match the dimensions of the aviris spectra
  result[use_wl - 399, aviris_inds]
}

############################################################
# Set up prior
############################################################
message("Setting up prior")
load("priors/mvtraits_priors.RData")
pfts <- rownames(means)

# Fix row and column names of means and covariances
mv_priors <- c(
  "prospect_N", "prospect_Cab", "prospect_Car",
  "prospect_Cw", "prospect_Cm", "SLA"
)
colnames(means) <- rownames(covars) <- colnames(covars) <- mv_priors

prior_sample <- function() {
  mv_draws <- map(
    pfts,
    ~mvtnorm::rmvnorm(1, means[.,], covars[, , .]) %>% setNames(mv_priors)
  )
  all_draws <- pmap(
    list(
      mv_draws,
      clumping_factor = runif(length(pfts), 0, 1),
      orient_factor = runif(length(pfts), -0.5, 0.5)
    ),
    c
  )
  names(all_draws) <- pfts
  c(unlist(all_draws), residual = rgamma(1, 0.01, 0.01))
}

prior_density <- function(params) {
  residual <- params["residual"]
  res_dens <- dgamma(residual, 0.01, 0.01, log = TRUE)
  traits <- params2edr(params, prospect = FALSE)$trait.values
  mvdens <- imap_dbl(
    traits,
    ~mvtnorm::dmvnorm(.x[mv_priors], means[.y, ], covars[, , .y], log = TRUE)
  )
  cf_dens <- dunif(map_dbl(traits, "clumping_factor"), 0, 1, log = TRUE)
  of_dens <- dunif(map_dbl(traits, "orient_factor"), -0.5, 0.5, log = TRUE)
  sum(mvdens, cf_dens, of_dens, res_dens)
}

prior <- BayesianTools::createPrior(
  density = prior_density,
  sampler = prior_sample
)
test_priors <- prior$density(prior$sampler())
stopifnot(!is.na(test_priors))

############################################################
# Run PDA
############################################################
message("Initializing PDA")
PEcAn.logger::logger.setLevel("INFO")

custom_settings <- list(
  init = list(iterations = 50),
  loop = list(iterations = 50),
  other = list(
    save_progress = here("ed-outputs/multi_site_pda/progress.rds")
  )
)
samples <- invert_bt(observed, model, prior, custom_settings = custom_settings)
