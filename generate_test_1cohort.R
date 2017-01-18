generate_test_1cohort <- function(dbh, pft, dens = 0.05) {

    # Single cohort
    c1 <- "1cohort"

    densstring <- paste0("dens", dens)
    denspath <- file.path(c1, densstring)

    dbhstring <- paste0("dbh", dbh)
    relpath <- file.path(denspath, dbhstring)

    # Create site, css, and pss files
    sitepath <- file.path(sites_dir, relpath)
    dir.create(sitepath, showWarnings = FALSE, recursive = TRUE)
    pft_number <- pfts[pft]
    sitepath_pft <- file.path(sitepath, pft)

    dir.create(sitepath_pft, showWarnings = FALSE)
    css <- css_common
    css[,"dbh"] <- dbh
    css[,"pft"] <- pft_number
    css[,"den"] <- dens
    prefix <- file.path(sitepath_pft, 
                        paste(pft, latlon.string, sep = "."))
    write.table(css, 
                file = paste0(prefix, ".css"), 
                sep = "     ", 
                row.names = FALSE, 
                col.names = TRUE, 
                quote = FALSE)
    file.symlink(site_path, paste0(prefix, ".site"))
    file.symlink(pss_path, paste0(prefix, ".pss"))

    # Add files to BETY
    new_inputs <- data.frame(name = paste0(basename(prefix), c('.site', '.css', '.pss')),
                             notes = 'EDR_testing.1cohort',
                             file_path = dirname(prefix),
                             site_id = site_id,
                             machine_id = machine_id,
                             format_id = c('site' = 10, 'css' = 15, 'patch' = 11)) %>%
    new_inputs$file_name <- new_inputs$name

    inputs <- db_merge_into(db = db, 
                            table = 'inputs',
                            values = new_inputs,
                            by = 'name',
                            id_colname = 'id')

    dbfiles <- db_merge_into(db = db,
                             table = 'dbfiles',
                             values = new_inputs,
                             by = c('file_name', 'file_path'),
                             id_colname = 'id')

    # TODO: Generate PEcAn XML
}

arg <- commandArgs(trailingOnly = TRUE)
print(arg)
dbh <- as.numeric(arg[1])
pft <- arg[2]
dens <- ifelse(is.na(arg[3]), 0.05, arg[3])

generate_test_1cohort(dbh, pft, dens)
