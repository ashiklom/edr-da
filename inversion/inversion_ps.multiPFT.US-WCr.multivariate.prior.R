#--------------------------------------------------------------------------------------------------#
#
#
# S. Serbin & A. Shiklomanov
#--------------------------------------------------------------------------------------------------#


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?

## Load functions
source("common.R")
library(mvtnorm)

## Load priors
load(file = normalizePath('priors/stan_priors_sun.RData')) # for sun exposed leaves only
priors <- priors_sun$means
PEcAn.utils::logger.info(" *** Single PFT ***  Running with sun exposed priors only")

#load(file = normalizePath('priors/prior_all.RData')) # based on leaves from sun/shaded positions
#priors <- priors_all$means
#PEcAn.utils::logger.info(" *** Single PFT ***  Running with priors using all leaves (sun and shade)")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Set user email address
email_add <- "sserbin@bnl.gov"

## Setup tag
dttag <- strftime(Sys.time(), "%Y%m%d_%H%M%S")

## Define PFT and canopy structure
pft <- c("temperate.Early_Hardwood","temperate.North_Mid_Hardwood","temperate.Late_Hardwood")
#pft <- list("temperate.Early_Hardwood","temperate.North_Mid_Hardwood","temperate.Late_Hardwood")
dens <- 0.015
dbh <- 20 # 20, 30 or 40
lai <- getvar("LAI_CO", dbh, pft)
num.cohorts <- 3
multi.pft <- "EMLH"

data_dir <- normalizePath(paste0('../run-ed/',num.cohorts,'cohort/dens',dens,'/dbh',dbh,'/',multi.pft))
paths <- list(ed2in = file.path(data_dir, 'ED2IN'),
              history = file.path(data_dir, 'outputs'))

## Setup PDA options
nchains <- 3
ngibbs.max <- 100000
ngibbs.min <- 500
ngibbs.step <- 1000

## Setup output
main_out <- paste("PDA", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
if (! file.exists(main_out)) dir.create(main_out,recursive=TRUE)
PEcAn.utils::logger.info(paste0("Running inversion in dir: ",main_out))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Setup PROSPECT for psuedo data
prospect_ver <- 5
pp <- list()
pp[["temperate.Early_Hardwood"]] <- c("N" = 1.4, "Cab" = 30, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)
pp[["temperate.North_Mid_Hardwood"]] <- c("N" = 1.8, "Cab" = 45, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)
pp[["temperate.Late_Hardwood"]] <- c("N" = 1.95, "Cab" = 65, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)
spectra_list <- list()
spectra_list[["temperate.Early_Hardwood"]] <- prospect(pp$temperate.Early_Hardwood, prospect_ver, TRUE)
spectra_list[["temperate.North_Mid_Hardwood"]] <- prospect(pp$temperate.North_Mid_Hardwood, prospect_ver, TRUE)
spectra_list[["temperate.Late_Hardwood"]] <- prospect(pp$temperate.Late_Hardwood, prospect_ver, TRUE)

## Setup EDR wavelengths and run date
par.wl = 400:2499
nir.wl = 2500
datetime <- ISOdate(2004, 07, 01, 16, 00, 00)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Output directory function
outdir_path <- function(runID) {
  #paste("inversion_prior", dttag, runID, sep = ".")
  paste0(main_out,"/inversion_prior.", dttag, ".",runID)
  
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Setup the output directories
run_first <- function(inputs) {
  outdir <- outdir_path(inputs$runID)
  dir.create(outdir, showWarnings = FALSE)
  try_link <- link_ed(outdir)
  
  trait.values <- list()
  for (i in seq_along(pft)) {
    trait.values[[pft[i]]] <- list()
  }
  
  albedo <- EDR(paths = paths,
                spectra_list = spectra_list,
                par.wl = par.wl,
                nir.wl = nir.wl,
                datetime = datetime,
                #trait.values = list(temperate.Late_Hardwood = list()),  
                trait.values = trait.values,
                output.path = outdir)
  return(albedo)
}
# Create initial output directory
first_run <- run_first(list(runID = 0))
#--------------------------------------------------------------------------------------------------#




