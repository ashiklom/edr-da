#' Setup EDR runs for a bunch of sites
#'
#' @param sites Character vector of site names
#' @param pda_dir EDR PDA directory (containing all sites)
#' @export
setup_edr_multisite <- function(sites, pda_dir) {
  ed2in_paths <- purrr::map(
    sites,
    ~list.files(file.path(pda_dir, .), "ED2IN", full.names = TRUE)
  )
  ed2in_sites <- purrr::map(ed2in_paths, read_ed2in)
  edr_paths <- file.path(pda_dir, sites, "edr")
  met_paths <- file.path(pda_dir, sites, "ED_MET_DRIVER_HEADER")
  stopifnot(all(file.exists(met_paths)))
  site_setup <- purrr::pmap(
    list(
      ed2in = ed2in_sites,
      output_dir = edr_paths,
      ED_MET_DRIVER_DB = met_paths
    ),
    setup_edr
  )
  site_setup
}
