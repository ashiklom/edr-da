#' Ribbon plot
#'
#' @param x Vector of x values
#' @param low Lower bounds of envelope
#' @param hi Upper bounds of envelope
#' @export
ribbon <- function(x, low, hi, ...) {
  stopifnot(
    length(x) == length(low),
    length(x) == length(hi)
  )
  x2 <- c(x, rev(x))
  y <- c(low, rev(hi))
  polygon(x2, y, ...)
}

#' Draw spectra confidence interval
#'
#' Based on empirical quantiles of reflectance spectra
#'
#' @param refl_mat Matrix of reflectance values (wavelengths in rows)
#' @param wl Vector of wavelengths
#' @param alpha Quantile probability (e.g. 0.05 means 95% CI)
#' @export
spec_ci <- function(refl_mat, wl = 400:2500, alpha = 0.05, ...) {
  stopifnot(
    nrow(refl_mat) == length(wl),
    alpha > 0,
    alpha < 1
  )
  a2 <- alpha / 2
  q1 <- apply(refl_mat, 1, quantile, a2, na.rm = TRUE)
  q2 <- apply(refl_mat, 1, quantile, 1 - a2, na.rm = TRUE)
  ribbon(wl, q1, q2, ...)
}
