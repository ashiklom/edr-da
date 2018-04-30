#' Load AVIRIS observations into a list
#'
#' @param sites Character vector of site tags
#' @param aviris_specfile Path to AVIRIS spectra H5 file. Default = "aviris/aviris.rds"
#' @param use_waves Wavelengths to use for observation. Default = 400 to 1300, by 10nm
#' @export
load_observations <- function(sites,
                              aviris_specfile = here::here("aviris/aviris.rds"),
                              use_waves = seq(400, 1300, by = 10)) {
  aviris_spec_all <- readRDS(aviris_specfile) / 10000
  aviris_plots <- colnames(aviris_spec_all)
  site_tags <- gsub("([[:alnum:]]+)_.*", "\\1", sites)
  av_inds <- aviris_plots %in% site_tags
  aviris_spec <- aviris_spec_all[, av_inds, drop = FALSE]
  colnames(aviris_spec) <- sites[match(colnames(aviris_spec), site_tags)]
  observed <- PEcAnRTM::resample(aviris_spec, use_waves)
  observed
}
