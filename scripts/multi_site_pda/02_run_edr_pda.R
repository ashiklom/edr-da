library(PEcAnRTM)
library(PEcAn.ED2)
library(redr)
library(here)
library(purrr)
library(foreach)
library(doParallel)
library(optparse)

parser <- OptionParser() %>%
  add_option("--hetero", action = "store_true", default = FALSE) %>%
  add_option("--fix_allom2", action = "store_true", default = FALSE) %>%
  add_option("--geo", action = "store_true", default = FALSE) %>%
  add_option("--ncores", action = "store", default = 1, type = "double") %>%
  add_option("--prefix", action = "store", default = "multi_site_pda_default", type = "character")

argl <- parse_args(parser)

# Configuration for geo
if (argl$geo) {
  img_path <- NULL
  edr_exe_path <- "/projectnb/dietzelab/ashiklom/ED2/EDR/build/ed_2.1-dbg"
} else {
  img_path <- "~/Projects/ED2/ed2.simg"
  edr_exe_path <- NULL
}

orig_pda_dir <- here("ed-outputs", "multi_site_pda")
pda_dir <- here("ed-outputs", argl$prefix)

message("Copying run data to: ", pda_dir)
stopifnot(file.copy(orig_pda_dir, pda_dir, recursive = TRUE))
message("Done!")

sites <- readLines(here("other_site_data", "site_list"))

observed <- load_observations(sites)
site_setup <- setup_edr_multisite(sites, pda_dir)

############################################################
# Define RTM inversion model
############################################################
message("Creating cluster with ", argl$ncores, " cores.")
cl <- makeCluster(argl$ncores)
registerDoParallel(cl)

lai_files <- file.path(dirname(site_setup), "lai_store")

# Returns a matrix of spectra
model <- function(params) {
  edr_in <- params2edr(params, prospect = TRUE, version = 5)
  walk(lai_files, ~cat("\n", file = .))   # always move to next line so file lines up with iterations
  pkgs <- c("PEcAnRTM", "PEcAn.ED2")
  exports <- c("img_path", "edr_in", "edr_exe_path")
  result <- tryCatch({
    foreach(s = site_setup, .combine = cbind, .packages = pkgs, .export = exports) %dopar% {
      outspec <- EDR(
        img_path = img_path,
        ed2in_path = s,
        spectra_list = edr_in$spectra_list,
        trait.values = edr_in$trait.values,
        verbose_error = FALSE,
        edr_exe_path = edr_exe_path
      )
      lai <- attr(outspec, "LAI")
      lai_file <- file.path(dirname(s), "lai_store")
      cat(lai, file = lai_file, append = TRUE)
      outspec
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
prior <- create_prior(fix_allom2 = argl$fix_allom2, heteroskedastic = argl$hetero)
message("Testing prior sampler and density functions.")
invalid <- check_prior(prior, error = TRUE, progress = FALSE)
message("Priors seem reliable")

message("Generating diagnostic prior plots")
pdf(here("ed-outputs", argl$prefix, "priors.pdf"))
for (j in seq_len(n_prior)) {
  hist(prior_samps[, j], main = names(p1)[j])
}
dev.off()

############################################################
# Test
############################################################
message("Done generating diagnostic prior plots")

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
# Remove old LAI files
invisible(walk(lai_files, file.remove))

############################################################
# Run PDA
############################################################
message("Initializing PDA")
PEcAn.logger::logger.setLevel("INFO")

custom_settings <- list(
  init = list(iterations = 500),
  loop = list(iterations = 200),
  other = list(
    save_progress = file.path(pda_dir, "progress.rds")
  )
)
samples <- invert_bt(observed, model, prior, custom_settings = custom_settings)
saveRDS(samples, file.path(pda_dir, "results.rds"))
