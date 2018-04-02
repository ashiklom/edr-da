#' Download GFDL data
#'
#' @param lat.in Latitude coordinate
#' @param lon.in Longitude coordinate
#' @param outfolder Target directory for storing met
#' @param start_run Run start date
#' @param end_run Run end date (inclusive)
#' @export
get_gfdl <- function(lat.in, lon.in, outfolder, start_run, end_run) {
  full_outdir <- list.files(
    dirname(outfolder),
    basename(outfolder),
    full.names = TRUE
  )
  if (length(full_outdir > 1) && file.exists(full_outdir)) {
    message("GFDL data already exists.")
    return(tibble(file_path = full_outdir))
  }
  PEcAn.data.atmosphere::download.GFDL(
    lat.in = lat.in,
    lon.in = lon.in,
    outfolder = outfolder,
    start_date = start_run,
    end_date = end_run,
    site_id = NULL
  )
}

#' Download CRUNCEP data
#'
#' @inheritParams get_gfdl
#' @export
get_cruncep <- function(lat.in, lon.in, outfolder, start_run, end_run) {
  PEcAn.data.atmosphere::download.CRUNCEP(
    outfolder = outfolder,
    start_date = start_run,
    end_date = end_run,
    site_id = NULL,
    lat.in = lat.in,
    lon.in = lon.in
  )
}

#' Convert meteorology to ED format
#'
#' @param rundir ED runtime directory
#' @param full_name Full path and prefix of meteorology files
#' @param latitude Latitude coordinate
#' @param longitude Longitude coordinate
#' @param overwrite Logical. If `TRUE`, overwrite existing met
#' @export
met2ed <- function(rundir, full_name, latitude, longitude, overwrite = TRUE) {
  met2model.ED2(
    in.path = dirname(full_name),
    in.prefix = gsub("\\.[[:digit:]]{4}\\.nc$", "", basename(full_name)),
    start_date = start_run,
    end_date = end_run,
    outfolder = rundir,
    lat = latitude,
    lon = longitude,
    overwrite = overwrite
  )
}
