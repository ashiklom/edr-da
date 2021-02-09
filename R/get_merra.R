#' Download MERRA flux variables
#'
#' @export
#' @author Alexey Shiklomanov
get_merra_date <- function(date, latitude, longitude, outdir, overwrite = FALSE) {
  date <- as.character(date)
  dpat <- "([[:digit:]]{4})-([[:digit:]]{2})-([[:digit:]]{2})"
  year <- as.numeric(gsub(dpat, "\\1", date))
  month <- as.numeric(gsub(dpat, "\\2", date))
  day <- as.numeric(gsub(dpat, "\\3", date))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  version <- if (year >= 2011) {
    400
  } else if (year >= 2001) {
    300
  } else {
    200
  }
  base_url <- "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2"

  lat_grid <- seq(-90, 90, 0.5)
  lon_grid <- seq(-180, 180, 0.625)
  ilat <- which.min(abs(lat_grid - latitude))
  ilon <- which.min(abs(lon_grid - longitude))
  idxstring <- sprintf("[0:1:23][%d][%d]", ilat, ilon)

  # For more on MERRA variables, see:
  # - The MERRA2 readme -- https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXRAD.5.12.4/doc/MERRA2.README.pdf
  # - The MERRA2 file spec -- https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
  # Page numbers below correspond to pages in the file spec.

  # Surface flux diagnostics (pg. 33)
  merra_prod <- "M2T1NXFLX.5.12.4"
  merra_file <- "tavg1_2d_flx_Nx"
  merra_vars <- tibble::tribble(
    ~CF_name, ~MERRA_name, ~units,
    # NIRDF - Surface downwelling nearinfared diffuse flux
    "surface_diffuse_downwelling_nearinfrared_radiative_flux_in_air", "NIRDF", "W/m2",
    # NIRDR - Surface downwelling nearinfrared beam flux
    "surface_direct_downwelling_nearinfrared_radiative_flux_in_air", "NIRDR", "W/m2"
  )

  # Standard variables
  url <- glue::glue(
    "{base_url}/{merra_prod}/{year}/{sprintf('%02d', month)}/",
    "MERRA2_{version}.{merra_file}.",
    "{year}{sprintf('%02d', month)}{sprintf('%02d', day)}.nc4.nc4"
  )
  qvars <- sprintf("%s%s", merra_vars$MERRA_name, idxstring)
  qstring <- paste(qvars, collapse = ",")
  outfile <- file.path(outdir, sprintf("merra-most-lat%.2f-lon%.2f-%d-%02d-%02d.nc",
                                       latitude, longitude, year, month, day))
  if (overwrite || !file.exists(outfile)) {
    req <- httr::GET(
      paste(url, qstring, sep = "?"),
      httr::authenticate(user = "pecanproject", password = "Data4pecan3"),
      httr::write_disk(outfile, overwrite = TRUE)
    )
  }

  # Land surface forcings (pg. 39)
  merra_lfo_prod <- "M2T1NXLFO.5.12.4"
  merra_lfo_file <- "tavg1_2d_lfo_Nx"
  merra_lfo_vars <- tibble::tribble(
    ~CF_name, ~MERRA_name, ~units,
    # Surface downwelling PAR diffuse flux, PARDF
    "surface_diffuse_downwelling_photosynthetic_radiative_flux_in_air", "PARDF", "W/m2",
    # Surface downwelling PAR beam flux, PARDR
    "surface_direct_downwelling_photosynthetic_radiative_flux_in_air", "PARDR", "W/m2"
  )

  # Land forcing
  url <- glue::glue(
    "{base_url}/{merra_lfo_prod}/{year}/{sprintf('%02d', month)}/",
    "MERRA2_{version}.{merra_lfo_file}.",
    "{year}{sprintf('%02d', month)}{sprintf('%02d', day)}.nc4.nc4"
  )
  qvars <- sprintf("%s%s", merra_lfo_vars$MERRA_name, idxstring)
  qstring <- paste(qvars, collapse = ",")
  outfile <- file.path(outdir, sprintf("merra-lfo-lat%.2f-lon%.2f-%d-%02d-%02d.nc",
                                       latitude, longitude, year, month, day))
  if (overwrite || !file.exists(outfile)) {
    req <- httr::GET(
      paste(url, qstring, sep = "?"),
      httr::authenticate(user = "pecanproject", password = "Data4pecan3"),
      httr::write_disk(outfile, overwrite = TRUE)
    )
  }
}
