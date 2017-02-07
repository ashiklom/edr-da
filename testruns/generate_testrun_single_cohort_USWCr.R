#--------------------------------------------------------------------------------------------------#
source('testruns/generate_testrun_functions.R')
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
pft <- arg[2]
dens <- ifelse(is.na(arg[3]), 0.05, arg[3])
generate_test_1cohort(sites_dir, test.sites[1], num.cohorts=1, dbh, pft, dens)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF