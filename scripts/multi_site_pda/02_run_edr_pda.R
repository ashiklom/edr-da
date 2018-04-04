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
  add_option("--prefix", action = "store", default = "multi_site_pda_test", type = "character")

argl <- parse_args(parser)
print(argl)

# Configuration for geo
if (argl$geo) {
  img_path <- NULL
  edr_exe_path <- "/projectnb/dietzelab/ashiklom/ED2/EDR/build/ed_2.1-opt"
} else {
  img_path <- "~/Projects/ED2/ed2.simg"
  edr_exe_path <- NULL
}

orig_pda_dir <- here("ed-outputs", "multi_site_pda")
pda_dir <- here("ed-outputs", argl$prefix)
dir.create(pda_dir, showWarnings = FALSE)

message("Copying run data to: ", pda_dir)
system2("rsync", c("-az", paste0(orig_pda_dir, "/"), pda_dir))
message("Done!")

sites <- readLines(here("other_site_data", "site_list"))

observed <- load_observations(sites)
waves <- wavelengths(observed)
site_setup <- setup_edr_multisite(sites, pda_dir)
aviris_inds <- match(colnames(observed), sites)

############################################################
# Define RTM inversion model
############################################################
message("Creating cluster with ", argl$ncores, " cores.")
cl <- makeCluster(argl$ncores)
registerDoParallel(cl)

lai_files <- map_chr(site_setup, dirname) %>% file.path(., "lai_store")
spec_files <- map_chr(site_setup, dirname) %>% file.path(., "spec_store")

# Returns a matrix of spectra
model <- function(params) {
  edr_in <- params2edr(params, prospect = TRUE, version = 5)
  walk(lai_files, ~cat("\n", file = ., append = TRUE))   # always move to next line so file lines up with iterations
  walk(spec_files, ~cat("\n", file = ., append = TRUE))
  pkgs <- c("PEcAnRTM", "PEcAn.ED2")
  exports <- c("img_path", "edr_in", "edr_exe_path")
  result <- tryCatch({
    foreach(s = site_setup, .combine = cbind, .packages = pkgs, .export = exports) %dopar% {
      outspec <- EDR(
        img_path = img_path,
        ed2in_path = s,
        spectra_list = edr_in$spectra_list,
        trait.values = edr_in$trait.values,
        verbose_error = TRUE,
        edr_exe_path = edr_exe_path
      )
      spec_file <- file.path(dirname(s), "spec_store")
      cat(outspec, sep = ",", file = spec_file, append = TRUE)
      lai <- attr(outspec, "LAI")
      lai_file <- file.path(dirname(s), "lai_store")
      cat(lai, sep = ",", file = lai_file, append = TRUE)
      outspec
    }
  }, error = function(e) {
    message("EDR failed for at least one site. Returning -1e10 and continuing.")
    NULL
  })
  # Make result match the dimensions of the aviris spectra
  if (!is.null(result)) {
    colnames(result) <- sites
    result[waves - 399, aviris_inds]
  } else {
    -1e10
  }
}

############################################################
# Set up prior
############################################################
message("Setting up prior")
prior <- create_prior(fix_allom2 = argl$fix_allom2, heteroskedastic = argl$hetero, verbose = TRUE)
message("Testing prior sampler and density functions.")
prior_samps <- check_prior(prior, error = FALSE, progress = FALSE)
if (attr(prior_samps, "n_invalid") > 15) {
  stop("More than 15 invalid priors. Something is probably wrong...")
}
message("Priors seem reliable")

message("Generating diagnostic prior plots")
pdf(here("ed-outputs", argl$prefix, "priors.pdf"))
for (j in seq_len(ncol(prior_samps))) {
  hist(prior_samps[, j], main = colnames(prior_samps)[j])
}
dev.off()

############################################################
# Test
############################################################

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

# Remove old LAI and spec files
invisible(walk(lai_files, file.remove))
invisible(walk(spec_files, file.remove))

############################################################
# Run PDA
############################################################
message("Initializing PDA")
PEcAn.logger::logger.setLevel("INFO")

custom_settings <- list(
  init = list(iterations = 500),
  loop = list(iterations = 200),
  other = list(
    save_progress = file.path(pda_dir, "progress.rds"),
    heteroskedastic = argl$hetero,
    threshold = 1.15
  )
)
samples <- invert_bt(observed, model, prior, custom_settings = custom_settings)
saveRDS(samples, file.path(pda_dir, "results.rds"))
