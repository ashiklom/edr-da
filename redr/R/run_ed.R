#' Run ED from a given directory
#' 
#' @param edr_exe_name Name of ED2 executable
#' @inheritParams exec_in_dir
#' @export
run_ed <- function(dir, edr_exe_name = 'ed_2.1') {
    exec_in_dir(dir, edr_exe_name)
}

#' Execute a system command in a directory
#'
#' @param dir Directory in which to execute command
#' @param exe Name of command to execute
#' @export
exec_in_dir <- function(dir, exe) {
    exe_string <- sprintf('(cd %s; ./%s)', dir, exe)
    system(exe_string, intern = TRUE)
}
