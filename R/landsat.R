#' Landsat band centers, in nm
#'
#' @export
l5bands <- c(485, 560, 660, 830, 1650, 2215)

#' @rdname l5bands
#' @export
l7bands <- c(485, 560, 660, 835, 1650, 2220)

#' @rdname l5bands
#' @export
l8bands <- c(443, 482, 561.5, 654.5, 865, 1608.5, 2200.5)

#' Convert full reflectance spectra to Landsat reflectance
#'
#' @param reflectance Vector of reflectance values
#' @return `tibble` containing wave
#' @export
spec2landsat <- function(reflectance) {
  lsats <- c("landsat5", "landsat7", "landsat8")
  lbands <- list(l5bands, l7bands, l8bands)
  bandnums <- map(lbands, seq_along)
  bandnums[[3]] <- bandnums[[3]] - 1  # Set Landsat 8 coastal band to 0
  lbandmat <- map2(bandnums, lbands, cbind)
  purrr::map(
    lsats,
    PEcAnRTM::spectral.response,
    spec = reflectance
  ) %>%
    purrr::map2(lbandmat, cbind) %>%
    purrr::map(`colnames<-`, c("value", "band", "wavelength")) %>%
    setNames(lsats) %>%
    purrr::map_dfr(tibble::as_tibble, .id = "landsat")
}
