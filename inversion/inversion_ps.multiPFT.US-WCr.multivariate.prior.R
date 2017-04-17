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
PEcAn.utils::logger.info(" *** Multi PFT ***  Running with sun exposed priors only")

edr.exe.name <- 'ed_2.1'

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
#pft <- c("temperate.Late_Hardwood","temperate.North_Mid_Hardwood")
pft <- list("temperate.Early_Hardwood","temperate.North_Mid_Hardwood","temperate.Late_Hardwood")
dens <- 0.015
dbh <- 20 # 20, 30 or 40
lai <- getvar("LAI_CO", dbh, pft)
num.cohorts <- 3
multi.pft <- "EMLH"

data_dir <- normalizePath(paste0('../run-ed/',num.cohorts,'cohort/dens',dens,'/dbh',dbh,'/',multi.pft))
paths <- list(ed2in = file.path(data_dir, 'ED2IN'),
              history = file.path(data_dir, 'outputs'))



#--------------------------------------------------------------------------------------------------#
## Setup PROSPECT for psuedo data
prospect_ver <- 5
pp <- list()
pp[["temperate.Early_Hardwood"]] <- c("N" = 1.4, "Cab" = 30, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)
pp[["temperate.Late_Hardwood"]] <- c("N" = 1.95, "Cab" = 65, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)
pp[["temperate.North_Mid_Hardwood"]] <- c("N" = 1.8, "Cab" = 45, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)

spectra_list <- list()
spectra_list[["temperate.Early_Hardwood"]] <- prospect(pp$temperate.Early_Hardwood, prospect_ver, TRUE)
spectra_list[["temperate.Late_Hardwood"]] <- prospect(pp$temperate.Late_Hardwood, prospect_ver, TRUE)
spectra_list[["temperate.North_Mid_Hardwood"]] <- prospect(pp$temperate.North_Mid_Hardwood, prospect_ver, TRUE)


## Setup EDR wavelengths and run date
par.wl = 400:2499
nir.wl = 2500
datetime <- ISOdate(2004, 07, 01, 16, 00, 00)
#--------------------------------------------------------------------------------------------------#

## Setup output
main_out <- paste("PDA", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
if (! file.exists(main_out)) dir.create(main_out,recursive=TRUE)
PEcAn.utils::logger.info(paste0("Running inversion in dir: ",main_out))


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
    trait.values[[pft[[i]]]] <- list()
  }
  
  albedo <- EDR(paths = paths,
                spectra_list = spectra_list,
                par.wl = par.wl,
                nir.wl = nir.wl,
                datetime = datetime,
                trait.values = trait.values,
                edr.exe.name = edr.exe.name,
                output.path = outdir)
  return(albedo)
}
# Create initial output directory
first_run <- run_first(list(runID = 0))
#--------------------------------------------------------------------------------------------------#

print(head(first_run))
print(tail(first_run))
print(range(first_run))

extract_params <- function(params) {
    list(prospect.params = 
            list('temperate.Early_Hardwood' = params[1:5],
                'temperate.North_Mid_Hardwood' = params[6:10],
                'temperate.Late_Hardwood' = params[11:15]),
         trait.values = 
             list('temperate.Early_Hardwood' = list('orient_factor' = params[16],
                                                    'clumping_factor' = params[17]),
                  'temperate.North_Mid_Hardwood' = list('orient_factor' = params[18],
                                                        'clumping_factor' = params[19]),
                  'temperate.Late_Hardwood' = list('orient_factor' = params[20],
                                                   'clumping_factor' = params[21])
                  )
         )
}

invert_model <- function(param, runID = 0) {

    outdir <- outdir_path(runID)
    paths_run <- list(ed2in = NA, history = outdir)

    # Parse parameters
    pars_list <- extract_params(param)
    spectra_list <- lapply(pars_list$prospect.params, prospect, version = 5, include.wl = TRUE)
    trait.values <- pars_list$trait.values

    albedo <- EDR(spectra_list = spectra_list,
                  trait.values = trait.values,
                  paths = paths_run,
                  par.wl = par.wl,
                  nir.wl = nir.wl,
                  datetime = datetime,
                  edr.exe.name = "ed_2.1",
                  output.path = outdir, 
                  change.history.time = FALSE)

    # Create quick figure
    #waves <- seq(400,2500,1)
    #png(paste(outdir,"/",'simulated_albedo.png',sep="/"),width=4900, height =2700,res=400)
    #par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
    #plot(waves,unlist(albedo)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
         #cex.axis=1.5, cex.lab=1.7,col="black")
    #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
           #"grey80")
    #lines(waves,unlist(albedo)*100,lwd=3, col="black")
    #dev.off()
    
    return(albedo)
}

# Simulate observations
inits <- c(unlist(pp[c('temperate.Early_Hardwood', 
                       'temperate.North_Mid_Hardwood', 
                       'temperate.Late_Hardwood')]),
           0.5, 0.5, # Early orient, clumping
           0.5, 0.5, # Mid orient, clumping
           0.5, 0.5) # Late orient, clumping

obs <- invert_model(inits, runID = 0) + generate.noise()

prior <- function(params) {
    do_prosp_prior <- function(params, pft) {
        dmvnorm(params,
                priors_sun$means$M[1:5, pft],
                priors_sun$means$Sigma[1:5, 1:5, pft],
                log = TRUE)
    }

    prospect_prior <- do_prosp_prior(params[1:5], 'temperate.Early_Hardwood') +
        do_prosp_prior(params[6:10], 'temperate.North_Mid_Hardwood') +
        do_prosp_prior(params[11:15], 'temperate.Late_Hardwood')

    traits_prior <- 
        # Early Hardwood
        dunif(params[16], 0, 1, log = TRUE) +   # Orient
        dunif(params[17], -1, 1, log = TRUE) +  # Clumping
        # North Mid Hardwood
        dunif(params[18], 0, 1, log = TRUE) +   # Orient
        dunif(params[19], -1, 1, log = TRUE) +  # Clumping
        # Late Hardwood
        dunif(params[20], 0, 1, log = TRUE) +   # Orient
        dunif(params[21], -1, 1, log = TRUE)    # Clumping

    return(prospect_prior + traits_prior)
}

init_function <- function() {
    c(1.1, 10, 2, 0.01, 0.01,
      1.1, 10, 2, 0.01, 0.01,
      1.1, 10, 2, 0.01, 0.01,
      0, 0,
      0, 0,
      0, 0)
}

# Test that prior function generates meaningful values
prior(inits)
prior(init_function())

param.mins <- c(1, 0, 0, 0, 0,
                1, 0, 0, 0, 0,
                1, 0, 0, 0, 0,
                0, 0,
                0, 0,
                0, 0)
param.maxs <- c(Inf, Inf, Inf, Inf, Inf,
                Inf, Inf, Inf, Inf, Inf,
                Inf, Inf, Inf, Inf, Inf,
                1, 1,
                1, 1,
                1, 1)

## Setup PDA options
invert.options <- list(model = invert_model,
                       run_first = run_first,
                       inits.function = init_function,
                       prior.function = prior,
                       param.mins = param.mins,
                       param.maxs = param.maxs,
                       nchains = 5,
                       ngibbs.step = 2000,
                       adapt = 1000,
                       do.lsq = FALSE)

samples <- invert.auto(observed = obs,
                       invert.options = invert.options,
                       parallel = TRUE,
                       save.samples = file.path(main_out, 'test_samples_ongoing.rds'))

saveRDS(samples, file = file.path(main_out, 'PDA_samples_output.rds')

samples.bt <- PEcAn.assim.batch::autoburnin(samples$samples)
samples.bt <- PEcAn.assim.batch::makeMCMCList(samples.bt)
par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
png(paste0(main_out,"/",paste("trace", runtag, "png", sep = ".")), width = 1500, height = 1600, res=150)
plot(samples.bt)
dev.off()

rawsamps <- do.call(rbind, samples.bt)
par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
png(paste0(main_out,"/",paste("pairs", runtag, "png", sep = ".")), width = 1500, height = 1600, res=150)
pairs(rawsamps)
dev.off()

par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
png(paste0(main_out,"/",paste("deviance", runtag, "png", sep = ".")), width = 1500, height = 1600, res=150)
plot(PEcAn.assim.batch::makeMCMCList(input.pda.data$deviance))
dev.off()

par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
png(paste0(main_out,"/",paste("n_eff", runtag, "png", sep = ".")), width = 1500, height = 1600, res=150)
plot(PEcAn.assim.batch::makeMCMCList(input.pda.data$n_eff_list))
dev.off()
