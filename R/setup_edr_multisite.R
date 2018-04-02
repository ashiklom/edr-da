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
  site_setup <- purrr::map2(ed2in_sites, edr_paths, setup_edr)
  site_setup
}
