library(dbhelpers)
library(dplyr)

db <- src_postgres(dbname = 'bety',
                   user = 'bety',
                   password = 'bety')

site_id <- 676          # Willow Creek (US-WCr)
machine_id <- 77        # My laptop machine ID

css_common <- data.frame(year = 2000,
                         patch = 1,
                         cohort = 1, 
                         dbh = as.numeric(NA), 
                         ht = 0, 
                         pft = as.numeric(NA), 
                         den = as.numeric(NA), 
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
