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
#load(file = normalizePath('priors/prior_all.RData')) # based on leaves from sun/shaded positions
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Set user email address
email_add <- "sserbin@bnl.gov"

## Setup tag
dttag <- strftime(Sys.time(), "%Y%m%d_%H%M%S")

## Define PFT and canopy structure
pft <- "temperate.Late_Hardwood"
dens <- 0.05
dbh <- 40 # 20, 30 or 40
lai <- getvar("LAI_CO", dbh, pft)

data_dir <- normalizePath(paste0('../run-ed/1cohort/dens',dens,'/dbh',dbh,'/',pft))
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
pp <- c("N" = 1.4, "Cab" = 30, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)
spectra_list <- list()
spectra_list[[pft]] <- prospect(pp, 5, TRUE)

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
  trait.values[[pft]] <- list()
  
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


#--------------------------------------------------------------------------------------------------#
# Setup model function for inversion. Need to update if adding additional traits
invert_model <- function(param, runID = 0) {
  
  outdir <- outdir_path(runID)
  paths_run <- list(ed2in = NA, history = outdir)
  
  prospect_params <- param[1:5]
  sla <- param[6]
  orient_factor <- param[7]
  spectra_list <- list()
  trait.values <- list()
  spectra_list[[pft]] <- prospect(pp, 5, TRUE)
  trait.values[[pft]] <- list(SLA = sla, orient_factor = orient_factor)
  
  albedo <- EDR(spectra_list = spectra_list,
                trait.values = trait.values,
                paths = paths_run,
                par.wl = par.wl,
                nir.wl = nir.wl,
                datetime = datetime,
                edr.exe.name = "ed_2.1-opt", # OK to change this from ed_2.1 to ed_2.1-opt??
                output.path = outdir, 
                change.history.time = FALSE)
  
  # Create quick figure
  waves <- seq(400,2500,1)
  png(paste(outdir,"/",'simulated_albedo.png',sep="/"),width=4900, height =2700,res=400)
  par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
  plot(waves,unlist(albedo)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
       cex.axis=1.5, cex.lab=1.7,col="black")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
         "grey80")
  lines(waves,unlist(albedo)*100,lwd=3, col="black")
  dev.off()
  
  return(albedo)
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Create pseudo data observation
sim_obs_params <- c(pp, "sla" = 18.85, "orient_factor" = 0.25)
#"clumping_factor" = 0.75,
#"b1Bl_large" = 0.05,
#"b2Bl_large" = 1.45)

obs <- invert_model(sim_obs_params) + generate.noise()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Setup prior function
priors <- priors_sun$means
prior_function <- function(params) {
  pft_priors <- dmvnorm(c(params[1:6]),
                           mean = priors$M[,paste(pft)], sigma = priors$Sigma[,,paste(pft)], log = TRUE)
  orient_factor_late <- dunif(params[7], -0.5, 0.5)
  return(sum(pft_priors+orient_factor_late,na.rm=T))
}
prior_function(rep(0,7))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Parameter ranges
param.mins <- c(N = 1, Cab = 1, Car = 0, Cw = 0.0001, Cm = 0.0001, sla = 5, orient_factor = -0.5)
param.maxs <- c(N = 6, Cab = 200, Car = 50, Cw = 0.99, Cm = 0.99, sla = Inf, orient_factor = 0.5)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Initialization function - Initial conditions
init_function <- function() return(c('N' = 1.5, 'Cab' = 30, 'Car' = 5, 'Cw' = 0.05, 'Cm' = 0.05, 'sla' = 20, 'orient_factor' = 0.2))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Complete inversion options
invert.options <- list(model = invert_model, 
                       run_first = run_first,
                       nchains = nchains,
                       inits.function = init_function,
                       prior.function = prior_function,
                       ngibbs.max = ngibbs.max,
                       ngibbs.min = ngibbs.min,
                       ngibbs.step = ngibbs.step,
                       param.mins = param.mins,
                       param.maxs = param.maxs,
                       adapt = 100,
                       adj_min = 0.1,
                       target = 0.234)
invert.options$do.lsq <- FALSE # TRUE/FALSE
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Run
runtag <- paste(paste0("dbh", dbh), pft, sep = ".")
fname <- paste("samples", runtag, "rds", sep = ".")
fname_prog <- paste("prog_samples", runtag, "rds", sep = ".")
logfile <- paste("output_prior", dttag, "log", sep = ".")
samples <- invert.auto(observed = obs, 
                       invert.options = invert.options,
                       parallel = TRUE,
                       parallel.output = paste0(main_out,"/",logfile),
                       save.samples = paste0(main_out,"/",fname_prog))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Output results
saveRDS(samples, file = paste0(main_out,"/",fname))
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

## send an email that the job is done
PEcAn.utils::sendmail(paste0(email_add),paste0(email_add),
                      paste0('PDA inversion in: ',main_out,' is complete'),
                      paste0("Log file can be found in: ", main_out,"/",logfile))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF