#--------------------------------------------------------------------------------------------------#
source('testruns/generate_testrun_functions.R')
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Create multi-cohort test runs
css_common.multicohort <- data.frame(year = 2000,
                                     patch = 1,
                                     cohort = 1:3,
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
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Generate test runs and ED2 output
run_dir <- normalizePath("run-ed")
template_dir <- file.path(run_dir, "template")

arg <- commandArgs(trailingOnly = TRUE)
print(arg)
dbh <- as.numeric(arg[1])
pft <- unlist(strsplit(arg[2], split=" "))
print(paste0("PFTs: ",pfts[pft]))
dens <- ifelse(is.na(arg[3]), 0.05, arg[3])
stand_type <- as.character(arg[4])
print(stand_type)
generate_test_multicohort(sites_dir, test.sites[1], length(css_common.multicohort$cohort),
                          dbh, pft, stand_type, dens=dens)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF