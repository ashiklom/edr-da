#' Adverb to add a progress bar to a function
#'
#' @param .f Function
#' @param total Total number of progress bar ticks
#' @param pb_opts Named list of other options to be passed to `progress_bar$new`
#' @return Function with progress bar
#' @export
with_prog <- function(.f, total, pb_opts = list()) {
  defaults <- getOption("prog_default")
  if (is.null(defaults)) defaults <- list()
  pb_opts <- modifyList(c(defaults, list(total = total)), pb_opts)
  .pb <<- do.call(progress::progress_bar$new, pb_opts)
  function(...) {
    .pb$tick()
    .f(...)
  }
}
