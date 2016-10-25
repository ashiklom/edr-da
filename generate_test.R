css_common <- data.frame(year = 2000,
                         patch = 1,
                         cohort = 1, 
                         dbh = NA, 
                         ht = 0, 
                         pft = NA, 
                         den = NA, 
                         bdead = 0,
                         balive = 0,
                         lai = -999)

sites_dir <- normalizePath("ed-inputs/sites/US-WCr/rtm")
site_path <- file.path(sites_dir, "common.site")
pss_path <- file.path(sites_dir, "common.pss")
latlon.string <- "lat45.5lon-90.5"

pfts <- c("temperate.Early_Hardwood" = 9,
          "temperate.Mid_Hardwood" = 10,
          "temperate.Late_Hardwood" = 11,
          "temperate.North_Pine" = 6,
          "temperate.Late_Conifer" = 8)

run_dir <- normalizePath("run-ed")
template_dir <- file.path(run_dir, "template")

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
                col.names = TRUE)
    file.symlink(site_path, paste0(prefix, ".site"))
    file.symlink(pss_path, paste0(prefix, ".pss"))

    # Create runtime files (ED2IN, ED executable)
    runpath <- file.path(run_dir, relpath)
    dir.create(runpath, showWarnings = FALSE, recursive = TRUE)

    runpath_pft <- file.path(runpath, pft)
    dir.create(runpath_pft, showWarnings = FALSE)
    system2("cp", c("-r", paste0(template_dir, "/*"), runpath_pft))
    ed2in_path <- file.path(runpath_pft, "ED2IN")
    fin <- readLines(ed2in_path)
    fout <- fin
    fout <- gsub("XXX_DENS_XXX", densstring, fout)
    fout <- gsub("XXX_PFT_XXX", pft, fout)
    fout <- gsub("XXX_DBH_XXX", dbhstring, fout)
    cat(fout, file = ed2in_path, sep = "\n")
}

arg <- commandArgs(trailingOnly = TRUE)
print(arg)
dbh <- as.numeric(arg[1])
pft <- arg[2]
dens <- ifelse(is.na(arg[3]), 0.05, arg[3])

generate_test_1cohort(dbh, pft, dens)
