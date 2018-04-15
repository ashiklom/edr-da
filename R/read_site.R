#' Get site prefix from site path and filename
#'
#' @param site_path Directory containing site css and pss files
#' @param which If multiple sites present, which to choose. Either "first" or "last".
#' @return Site prefix, as character
#' @export
get_site_prefix <- function(site_path, which = "first") {
  stopifnot(which %in% c("first", "last"))
  opts <- list.files(site_path, "\\.css$")
  ind <- switch(which, `first` = 1, `last` = length(opts))
  fname <- opts[ind]
  file_prefix <- gsub("(FFT\\.[[:digit:]]{4}\\.).*", "\\1", fname)
  file.path(site_path, file_prefix)
}
