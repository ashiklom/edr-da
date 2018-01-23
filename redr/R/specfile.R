#' Extract spectral data from HDF5 file based on criteria
#'
#' @param aviris_h5 HDF5 file object
#' @param ... Filters, similar to those used in `dplyr::filter`
#'
#' @return `spectra` object (matrix) containing AVIRIS data
#' @export
filter_specfile <- function(aviris_h5, ...) {
  dots <- rlang::quos(...)
  filters <- purrr::map(dots, ~apply_filter_specfile(aviris_h5, .))
  filter_combined <- purrr::reduce(filters, combine_logical)
  aviris_spec <- aviris_h5[["spectral_data"]][, filter_combined]
  aviris_wl <- aviris_h5[["wavelengths"]][]
  out <- PEcAnRTM::spectra(aviris_spec, aviris_wl)
  colnames(out) <- names(filter_combined)[filter_combined]
  out
}

#' Apply single filter to spectra HDF5 data file
#'
#' @inheritParams filter_specfile
#' @param quosure An `rlang` `quosure` object (e.g. a `formula`)
apply_filter_specfile <- function(aviris_h5, quosure) {
  expr <- quosure[[2]]
  fun <- match.fun(as.character(expr[[1]]))
  varname <- as.character(expr[[2]])
  arg <- rlang::eval_tidy(expr[[3]])
  full_attr <- aviris_h5[[varname]][]
  out <- fun(full_attr, arg)
  names(out) <- full_attr
  out
}

#' Combine two logical vectors
#'
#' Is identical to `&`, except that it also concatenates the names of the 
#' vectors.
#' @param x Named logical vector
#' @param y Named logical vector
#' @return Logical vector `x & y`, with concatenated names
combine_logical <- function(x, y, sep = ".", ...) {
  xy <- x & y
  names(xy) <- paste(names(x), names(y), sep = sep, ...)
  xy
}
