#' Refresh current site list based on ED runs and available spectra
#'
#' @param pda_dir EDR multi-site PDA directory
#' @param drop_sites Character vector of sites to drop from list.
#' @param outfile Target file for site list. Default = "other_site_data/site_list"
#' @inheritParams load_observations
#' @export
refresh_site_list <- function(pda_dir,
                              aviris_specfile = here::here("aviris/aviris.rds"),
                              drop_sites = "IDS05_site_1-25682",
                              outfile = here::here("other_site_data/site_list")) {
  stopifnot(file.exists(pda_dir))
  if (!is.null(aviris_specfile)) stopifnot(file.exists(aviris_specfile))
  all_sites <- list.files(pda_dir, "_site_")
  stopifnot(length(all_sites) > 0)

  aviris_h5 <- hdf5r::H5File$new(aviris_specfile)
  on.exit(aviris_h5$close_all())
  site_tags <- gsub("([[:alnum:]]+)_.*", "\\1", all_sites)
  av_inds <- purrr::map(site_tags, ~which(aviris_h5[["iPLOT"]][] == .))
  keep_sites <- purrr::map_lgl(av_inds, ~length(.) > 0) &
    !(all_sites %in% drop_sites)

  sites <- all_sites[keep_sites]
  writeLines(sites, outfile)
}
