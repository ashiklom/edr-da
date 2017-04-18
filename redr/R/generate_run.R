#' Prepare ED run inputs in a single location
#'
#' @param prefix Input file prefix, matched by ED2IN (full filename is `prefix.lat(-)DD.Dlon(-)DD.D.{css,pss,site}`
#' @param site_lat Site latitude (decimal degrees)
#' @param site_lon Site longitude (decimal degrees)
#' @param css_df `data.frame` contianing cohort (.css) file
#' @param pss_df `data.frame` containing patch (.pss) file
#' @param site_df `data.frame` containing site (.site) file
#' @param output_dir Root directory where outputs are stored
#' @param common_inputs_dir Path to directory containing `ed-inputs/{chd,dgd}` and `OGE2old`
#' @param site_met_dir Path to directory containing `ED_MET_DRIVER_HEADER` and `{month}.h5` files
#' @param ed_exe_path Path to ED executable, which will be linked to directory
#' @param ed2in_template Path to template ED2IN file, which will be copied and modified
#' @param ed2in_changes Optional list of additional changes to make to ED2IN, in the form `TAG = "value"`
#' @param RMDIR If TRUE, remove output directory and contents before starting
#' @export
generate_run <- function(prefix, site_lat, site_lon, css_df, pss_df, site_df, 
                         output_dir, common_inputs_dir, site_met_dir, ed_exe_path,
                         ed2in_template, ed2in_changes = NULL, RMDIR = FALSE) {
    if (RMDIR) {
        unlink(output_dir, recursive = TRUE)
    }

    dir.create(output_dir, showWarnings = FALSE)
    common_inputs_link <- file.path(output_dir, 'common-inputs')
    site_met_link <- file.path(output_dir, 'site-met')
    sites_dir <- file.path(output_dir, 'site-inputs')
    ed_output_dir <- file.path(output_dir, 'outputs')


    # Symlink input files
    file.symlink(from = common_inputs_dir, to = common_inputs_link)
    file.symlink(from = site_met_dir, to = site_met_link)

    # Write PSS, CSS, and SITE files to output_dir/site_files
    dir.create(sites_dir)
    latlon_string <- sprintf('lat%.1flon%.1f', site_lat, site_lon)
    write_edfile <- function(df, extension, ...) {
        fname <- file.path(sites_dir, paste(prefix, latlon_string, extension, sep = '.'))
        append <- FALSE
        # Site file has special syntax
        if (extension == 'site') {
            site_header <- sprintf('nsite %d format 1', nrow(site_df))
            writeLines(text = site_header, con = fname)
            append <- TRUE
        }
        write.table(x = df,
                    file = fname,
                    quote = FALSE,
                    sep = "     ",
                    row.names = FALSE,
                    col.names = TRUE,
                    append = append)
    }
    write_edfile(css_df, 'css')
    write_edfile(pss_df, 'pss')
    write_edfile(site_df, 'site')

    # Write ED2IN file
    ed2in_text <- readLines(ed2in_template)
    ed2in_values <- list('VEG_DATABASE' = file.path(normalizePath(file.path(common_inputs_link, 'oge2OLD')), 'OGE2_'),
                         'THSUMS_DATABASE' = paste0(normalizePath(file.path(common_inputs_link, 'EDI/ed_inputs')), '/'),
                         'ED_MET_DRIVER_DB' = normalizePath(file.path(site_met_link, 'ED_MET_DRIVER_HEADER')),
                         'SFILIN' = file.path(normalizePath(sites_dir), paste0(prefix, '.')))
    if (!is.null(ed2in_changes)) {
        ed2in_values <- modifyList(ed2in_values, ed2in_changes)
    }
    ed2in_text <- PEcAn.ED2::ed2in_set_value_list(ed2in_values, ed2in = ed2in_text)
    ed2in_output <- file.path(output_dir, 'ED2IN')
    writeLines(ed2in_text, ed2in_output)

    # Symlink ED2 executable
    file.symlink(from = ed_exe_path, to = file.path(output_dir, 'ed_2.1'))

    # Create ED outputs directory
    dir.create(ed_output_dir)
}
