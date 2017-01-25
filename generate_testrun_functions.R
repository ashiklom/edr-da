#--------------------------------------------------------------------------------------------------#
## !! could spend some time here to generalize into a single function where number of
## !! cohorts is an input

## Generate test functions - 1 cohort 
generate_test_1cohort <- function(path, site, num.cohorts=1, dbh, pft, dens = 0.05) {
  
  sites_dir <- path
  #print(paste0("Site input path: ",sites_dir))
  
  # Single cohort
  num.cohorts <- paste0(num.cohorts,"cohort")
  c1 <- num.cohorts
  
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
  css <- css_common.1cohort
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
  
  # Create runtime files (ED2IN, ED executable)
  runpath <- file.path(run_dir, relpath)
  dir.create(runpath, showWarnings = FALSE, recursive = TRUE)
  
  runpath_pft <- file.path(runpath, pft)
  dir.create(runpath_pft, showWarnings = FALSE)
  system2("cp", c("-r", paste0(template_dir, "/*"), runpath_pft))
  ed2in_path <- file.path(runpath_pft, "ED2IN")
  fin <- readLines(ed2in_path)
  fout <- fin
  fout <- gsub("XXX_SITE_XXX", site, fout)
  fout <- gsub("XXX_NUMCOHORTS_XXX", num.cohorts, fout)
  fout <- gsub("XXX_DENS_XXX", densstring, fout)
  fout <- gsub("XXX_PFT_XXX", pft, fout)
  fout <- gsub("XXX_DBH_XXX", dbhstring, fout)
  cat(fout, file = ed2in_path, sep = "\n")
}

## Generate test functions - multi-cohort
generate_test_multicohort <- function(path, site, num.cohorts=3, dbh, pft, stand_type, dens = 0.05) {
  
  sites_dir <- path
  #print(paste0("Site input path: ",sites_dir))
  
  # Multi cohort
  num.cohorts <- paste0(num.cohorts,"cohort")
  c1 <- num.cohorts
  
  densstring <- paste0("dens", dens)
  denspath <- file.path(c1, densstring)
  
  dbhstring <- paste0("dbh", dbh)
  relpath <- file.path(denspath, dbhstring)
  
  # Create site, css, and pss files
  sitepath <- file.path(sites_dir, relpath)
  print(paste0("Sitepath: ",sitepath))
  dir.create(sitepath, showWarnings = FALSE, recursive = TRUE)
  pft_number <- pfts[pft]
  sitepath_pft <- file.path(sitepath, paste(stand_type))
  
  dir.create(sitepath_pft, showWarnings = FALSE)
  css <- css_common.multicohort
  css[,"dbh"] <- dbh
  css[,"pft"] <- pft_number
  css[,"den"] <- dens
  prefix <- file.path(sitepath_pft, 
                      paste("EMLH", latlon.string, sep = ".")) # hacky!
  print(prefix)
  
  write.table(css, 
              file = paste0(prefix, ".css"), 
              sep = "     ", 
              row.names = FALSE, 
              col.names = TRUE, 
              quote = FALSE)
  file.symlink(site_path, paste0(prefix, ".site"))
  file.symlink(pss_path, paste0(prefix, ".pss"))
  
  # Create runtime files (ED2IN, ED executable)
  runpath <- file.path(run_dir, relpath)
  print(runpath)
  dir.create(runpath, showWarnings = FALSE, recursive = TRUE)
  
  runpath_pft <- file.path(runpath, paste(stand_type)) # hacky
  dir.create(runpath_pft, showWarnings = FALSE)
  system2("cp", c("-r", paste0(template_dir, "/*"), runpath_pft))
  ed2in_path <- file.path(runpath_pft, "ED2IN")
  fin <- readLines(ed2in_path)
  fout <- fin
  fout <- gsub("XXX_SITE_XXX", site, fout)
  fout <- gsub("XXX_NUMCOHORTS_XXX", num.cohorts, fout)
  fout <- gsub("XXX_DENS_XXX", densstring, fout)
  #fout <- gsub("XXX_PFT_XXX", pft, fout)
  fout <- gsub("XXX_PFT_XXX", stand_type, fout)
  fout <- gsub("XXX_DBH_XXX", dbhstring, fout)
  cat(fout, file = ed2in_path, sep = "\n")
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## EOF