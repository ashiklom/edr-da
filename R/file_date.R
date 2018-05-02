#' Extract date from an ED output file name
#'
#' @param fname File name
#' @return Date
#' @export
file_date <- function(fname) {
  rxp <- "([[:digit:]]{4})-([[:digit:]]{2})-([[:digit:]]{2})"
  hits <- stringr::str_match(fname, rxp)
  y <- as.numeric(hits[, 2])
  m <- as.numeric(hits[, 3])
  m[m == 0] <- 1
  d <- as.numeric(hits[, 4])
  d[d == 0] <- 1
  lubridate::as_date(paste(y, m, d, sep = "-"))
}
