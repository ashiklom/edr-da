get_gfdl <- function(lat.in, lon.in, outfolder) {
  stopifnot(exists("start_run"), exists("end_run"))
  full_outdir <- list.files(
    dirname(outfolder),
    basename(outfolder),
    full.names = TRUE
  )
  if (length(full_outdir > 1) && file.exists(full_outdir)) {
    message("GFDL data already exists.")
    return(tibble(file_path = full_outdir))
  }
  download.GFDL(
    lat.in = lat.in,
    lon.in = lon.in,
    outfolder = outfolder,
    start_date = start_run,
    end_date = end_run,
    site_id = NULL
  )
}

get_cruncep <- function(lat.in, lon.in, outfolder) {
  stopifnot(exists("start_run"), exists("end_run"))
  download.CRUNCEP(
    outfolder = outfolder,
    start_date = start_run,
    end_date = end_run,
    site_id = NULL,
    lat.in = lat.in,
    lon.in = lon.in
  )
}
