#' Run EDR for a single date
#' 
#' @param date Date at which to run EDR
#' @param ed2in ED2IN list to use as ED input
#' @param trait_values Named list of trait values
#' @param output_dir Directory where to store outputs
#' @param pb Optional progress bar object
#' @inheritParams PEcAnRTM::EDR
#' @export
run_edr_date <- function(date, ed2in, trait_values,
                         output_dir = tempdir(),
                         pb = NULL,
                         img_path = NULL,
                         edr_exe_path = file.path("/projectnb/dietzelab/ashiklom",
                                                  "ED2/EDR/build",
                                                  "ed_2.1-opt")
                                                  ) {
  on.exit(if (!is.null(pb)) pb$tick())
  if (!lubridate::is.Date(date)) date <- lubridate::as_date(date)
  dtime <- ISOdatetime(
    lubridate::year(date),
    lubridate::month(date),
    lubridate::mday(date),
    12, 00, 00, tz = "UTC"
  )
  edr_ed2in <- PEcAnRTM::setup_edr(ed2in, output_dir, dtime)
  ptraits <- purrr::map(trait_values, ~.[grepl("prospect", names(.))])
  spectra_list <- purrr::map(ptraits, PEcAnRTM::prospect, version = "5")
  PEcAnRTM::EDR(
    img_path = img_path,
    ed2in_path = edr_ed2in,
    spectra_list = spectra_list,
    trait.values = trait_values,
    edr_exe_path = edr_exe_path
  )
}
