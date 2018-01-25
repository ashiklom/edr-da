#' Retrieve variable from ED history file
#'
#' @param prefix File path prefix
#' @param varname Name of variable in HDF5 file
#' @param edr_dir Name of EDR output subdirectory (default = "edr")
#' @param hist_prefix Prefix of history file name (default = "history-S")
get_edvar <- function(prefix, varname, edr_dir = "edr", hist_prefix = "history-S") {
  full_path <- here::here(prefix, edr_dir)
  file_name <- list.files(full_path, hist_prefix)
  hfile <- hdf5r::H5File$new(file.path(full_path, file_name))
  on.exit(hfile$close_all())
  hfile[[varname]][]
}

#' Prepare an FFT plot for testing
#'
#' Does the following:
#' - Identify css, pss, and site files
#' - Run ED at that site
#'
#' @param fft_plot FFT plot code, as character
#' @param meas_year Measurement year, as integer or character
#' @param clobber Overwrite existing ED run if present (default = TRUE)
#' @return ED run directory prefix, as character
setup_fft <- function(fft_plot, meas_year, clobber = TRUE) {
  ed_inputs_dir <- normalizePath(Sys.glob(here::here('ed-inputs/sites', paste0(fft_plot,"_*"))) )
  ed_inputs <- list.files(ed_inputs_dir)
  tmp <- strsplit(ed_inputs_dir, "ed-inputs")[[1]]
  edr_dir <- tmp[1]

  site_lat <- as.numeric(gsub('^.*lat\\s*|\\s*lon.*$', '', ed_inputs[1]))
  site_lon <- as.numeric(gsub('^.*lon\\s*|\\s*.css.*$', '', ed_inputs[1]))
  meas_year <- as.numeric(gsub('^.*FFT.\\s*|\\s*.lat.*$', '', ed_inputs[1]))
  
  # link to css, pss, site files
  lsf <- partial(
    list.files,
    path = ed_inputs_dir,
    recursive = TRUE,
    full.names = TRUE
  )
  css_in <- lsf(pattern = "*.css")
  pss_in <- lsf(pattern = "*.pss")
  site_in <- lsf(pattern = "*.site")
  
  ### create output prefix
  prefix <- paste("EDR_sim_output", fft_plot, meas_year, sep = "_")
  if (file.exists(paste0(edr_dir, prefix)) && clobber) {
    PEcAn.logger::logger.info("*** Removing previous simulation results ***")
    unlink(paste0(edr_dir, prefix), recursive = TRUE)
  }
  if (!file.exists(paste0(edr_dir, prefix))) {
    PEcAn.logger::logger.info(paste0("Running simulation in dir: ", prefix))
    
    ed2in_changes <- list(IMONTHA = 07, IDATEA = 15, IYEARA = 2006,
                          IMONTHZ = 08, IDATEZ = 31, IYEARZ = 2006)
    datetime <- ISOdate(2006, 08, 15, 16, 00, 00)
    
    genrun <- generate_run(prefix = prefix,
                          site_lat = site_lat,
                          site_lon = site_lon,
                          site_df = site_in,
                          pss_df = pss_in,
                          css_df = css_in,
                          common_inputs_dir = common_inputs_dir,
                          site_met_dir = site_met_dir,
                          ed_exe_path = ed_exe_path,
                          ed2in_changes = ed2in_changes,
                          RMDIR = TRUE)
    
    message('Running ED...')
    runed <- run_ed(prefix)
    tail(runed)
    message('Done!')
  } else {
    PEcAn.logger::logger.info("Skipping run because results present in dir: ", prefix)
  }

  list(
    prefix = prefix,
    css = css_in,
    pss = pss_in,
    site = site_in
  )
}

### Setup helper functions
pft_lookup <- tribble(
  ~pft_num, ~pft_name,
  6, "temperate.Northern_Pine",
  7, "temperate.Southern_Pine",
  8, "temperate.Late_Conifer",
  9, "temperate.Early_Hardwood",
  10, "temperate.North_Mid_Hardwood",
  11, "temperate.Late_Hardwood"
)

get_pfts <- function(css_file_path, delim = ' ', reference = pft_lookup) {
  css_data <- read.delim(css_file_path, header = T, sep = delim)
  css_pfts <- as.vector(unique(css_data["pft"]))
  filter(pft_lookup, pft_num %in% css_pfts$pft)
}

spec_ci <- function(refl_mat, ...) {
  wl <- 400:2500
  x <- c(wl, rev(wl))
  mu <- rowMeans(refl_mat)
  q1 <- apply(refl_mat, 1, quantile, 0.025)
  q2 <- apply(refl_mat, 1, quantile, 0.975)
  qv <- c(q1, rev(q2))
  polygon(x, qv, ...)
  #lines(wl, mu, col = "grey")
}

