#' Read an `Rdata` file into a list (rather than into the global
#' environment)
#'
#' @param file Target file name
#' @return Named list containing objects in `file`.
#' @author Alexey Shiklomanov
#' @examples
#' \dontrun{
#' mylist <- load_local("file.RData")
#' }
#' @export
load_local <- function(file) {
  menv <- new.env(parent = baseenv())
  load(file, envir = menv)
  as.list(menv)
}
