library(PEcAnRTM)
library(PEcAn.ED2)
library(here)
library(purrr)
library(foreach)
library(doParallel)

# Configuration for geo
img_path <- NULL
edr_exe_path <- "/projectnb/dietzelab/ashiklom/ED2/EDR/build/ed_2.1-dbg"

#if (hostname == "ashiklom") {
  #img_path <- "~/Projects/ED2/ed2.simg"
  #edr_exe_path <- NULL
#}

pda_dir <- here("ed-outputs", "multi_site_pda_allom")
source("scripts/multi_site_pda/setup_sites.R")

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
source("scripts/multi_site_pda/priors.R")

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
