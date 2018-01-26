#' Retrieve variable from ED history file
#'
#' @param prefix File path prefix
#' @param varname Name of variable in HDF5 file
#' @param edr_dir Name of EDR output subdirectory (default = "edr")
#' @param hist_prefix Prefix of history file name (default = "history-S")
#' @export
get_edvar <- function(prefix, varname, edr_dir = "edr", hist_prefix = "history-S") {
  full_path <- here::here(prefix, edr_dir)
  file_name <- list.files(full_path, hist_prefix)
  hfile <- hdf5r::H5File$new(file.path(full_path, file_name))
  on.exit(hfile$close_all())
  hfile[[varname]][]
}
