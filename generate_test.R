#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## !! could spend some time here to generalize into a single function where number of
## !! cohorts is an input

## Generate test functions - 1 cohort 
generate_test_1cohort <- function(path, dbh, pft, dens = 0.05) {
  
  sites_dir <- path
  print(paste0("Path: ",sites_dir))
  
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
  fout <- gsub("XXX_DENS_XXX", densstring, fout)
  fout <- gsub("XXX_PFT_XXX", pft, fout)
  fout <- gsub("XXX_DBH_XXX", dbhstring, fout)
  cat(fout, file = ed2in_path, sep = "\n")
}

## Generate test functions - multi-cohort
# placeholder
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Create 1 cohort test runs
css_common.1cohort <- data.frame(year = 2000,
                         patch = 1,
                         cohort = 1, 
                         dbh = as.numeric(NA), 
                         ht = 0, 
                         pft = as.numeric(NA), 
                         den = as.numeric(NA), 
                         bdead = 0,
                         balive = 0,
                         lai = -999)

## Create multi-cohort test runs
css_common.multicohort <- data.frame(year = 2000,
                         patch = 1,
                         cohort = 3,
                         dbh = as.numeric(NA),
                         ht = 0,
                         pft = as.numeric(NA),
                         den = as.numeric(NA),
                         bdead = 0,
                         balive = 0,
                         lai = -999)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Sites
test.sites <- c("US-WCr","US-Syv")

# US-WCr
sites_dir <- normalizePath(paste0("ed-inputs/sites/",test.sites[1],"/rtm"))
site_path <- file.path(sites_dir, "common.site")
pss_path <- file.path(sites_dir, "common.pss")
latlon.string <- "lat45.5lon-90.5"

pfts <- c("temperate.Early_Hardwood" = 9,
          "temperate.North_Mid_Hardwood" = 10,
          "temperate.Late_Hardwood" = 11)

# US-Syv
#sites_dir <- normalizePath(paste0("ed-inputs/sites/",test.sites[2],"/rtm"))
#site_path <- file.path(sites_dir, "common.site")
#pss_path <- file.path(sites_dir, "common.pss")
#latlon.string <- "lat46.5lon-89.5"
#pfts <- c("temperate.Early_Hardwood" = 9,
#          "temperate.North_Mid_Hardwood" = 10,
#          "temperate.Late_Hardwood" = 11,
#          "temperate.Northern_Pine" = 6,
#          "temperate.Late_Conifer" = 8)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Generate test runs and ED2 output
run_dir <- normalizePath("run-ed")
template_dir <- file.path(run_dir, "template")

arg <- commandArgs(trailingOnly = TRUE)
print(arg)
dbh <- as.numeric(arg[1])
pft <- arg[2]
dens <- ifelse(is.na(arg[3]), 0.05, arg[3])

generate_test_1cohort(sites_dir,dbh, pft, dens)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## EOF