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
av_inds <- purrr::map(site_tags, ~which(aviris_h5[["iPLOT"]][] == .))
# Subset sites to those that have observations
# Problem sites were causing trouble with the inversion
#problem_sites <- c("BH02_site_1-25665", "IDS05_site_1-25682")
problem_sites <- c("IDS05_site_1-25682")
keep_sites <- purrr::map_lgl(av_inds, ~length(.) > 0) &
  !(all_sites %in% problem_sites)

sites <- all_sites[keep_sites]

# Grab AVIRIS data
aviris_spec_raw <- purrr::map(av_inds[keep_sites], ~aviris_h5[["spectral_data"]][, .])
aviris_wl <- aviris_h5[["wavelengths"]][]
aviris_h5$close_all()
aviris_spec <- purrr::map(aviris_spec_raw, spectra, wavelengths = aviris_wl)
names(aviris_spec) <- sites

# Resample AVIRIS data to integer wavelengths, and keep only VIS and NIR
use_wl <- seq(400, 1300, by = 10)
aviris_spec <- purrr::map(aviris_spec, resample, to = use_wl)
aviris_nspec <- purrr::map_int(aviris_spec, ncol)
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
