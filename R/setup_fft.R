#' Prepare an FFT plot for testing
#'
#' Does the following:
#' - Identify css, pss, and site files
#' - Run ED at that site
#'
#' @param fft_plot FFT plot code, as character
#' @param meas_year Measurement year, as integer or character
#' @param clobber Overwrite existing ED run if present (default = TRUE)
#' @return List, containing prefix path and paths to css, pss, and site files
#' @export
setup_fft <- function(fft_plot, meas_year, clobber = TRUE) {
  ed_inputs_dir <- normalizePath(Sys.glob(here::here('ed-inputs/sites', paste0(fft_plot, "_*"))) )
  ed_inputs <- list.files(ed_inputs_dir)
  edr_dir <- strsplit(ed_inputs_dir, "ed-inputs")[[1]][1]

  site_lat <- as.numeric(gsub("^.*lat\\s*|\\s*lon.*$", "", ed_inputs[1]))
  site_lon <- as.numeric(gsub("^.*lon\\s*|\\s*.css.*$", "", ed_inputs[1]))
  meas_year <- as.numeric(gsub("^.*FFT.\\s*|\\s*.lat.*$", "", ed_inputs[1]))

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

    message("Running ED...")
    runed <- run_ed(prefix)
    tail(runed)
    message("Done!")
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
