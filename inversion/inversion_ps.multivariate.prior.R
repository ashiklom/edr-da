#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?

## Load functions
source("common.R")
library(mvtnorm)

## Load priors
load(file = normalizePath('priors/stan_priors_sun.RData'))
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
## Setup PROSPECT
pp <- c(1.4, 30, 8, 0.01, 0.01)
spectra_list <- list(temperate.Late_Hardwood = prospect(pp, 5, TRUE))

#paths <- getpaths(dbh, pft)
par.wl = 400:2499
nir.wl = 2500
datetime <- ISOdate(2004, 07, 01, 16, 00, 00)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
outdir_path <- function(runID) {
  #paste("inversion_prior", dttag, runID, sep = ".")
  paste0(main_out,"/inversion_prior.", dttag, ".",runID)
  
}

get_trait_values <- function(param) {
  orient_factor <- param[1]
  #clumping_factor <- param[2]
  #sla <- param[3]
  #b1Bl_large <- param[4]
  #b2Bl_large <- param[5]
  
  trait.values <- list()
  trait.values[[pft]] <- list(orient_factor = orient_factor)
  #trait.values[[pft]] <- list(orient_factor = orient_factor,
  #                            clumping_factor = clumping_factor,
  #                            SLA = sla)
  #b1Bl_large = b1Bl_large,
  #b2Bl_large = b2Bl_large)
  return(trait.values)
}

# setup the output directories
run_first <- function(inputs) {
  outdir <- outdir_path(inputs$runID)
  dir.create(outdir, showWarnings = FALSE)
  try_link <- link_ed(outdir)
  
  albedo <- EDR(paths = paths,
                spectra_list = spectra_list,
                par.wl = par.wl,
                nir.wl = nir.wl,
                datetime = datetime,
                trait.values = list(temperate.Late_Hardwood = list()),  # hacky, need to allow this to be set by pft above
                output.path = outdir)
  return(albedo)
}


invert_model <- function(param, runID = 0) {
  
  outdir <- outdir_path(runID)
  paths_run <- list(ed2in = NA, history = outdir)
  
  trait.values <- get_trait_values(param) 
  
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

# Simulate observation
inits <- c("orient_factor" = 0.25)
#           "clumping_factor" = 0.75,
#           "sla" = 18.85)
#"b1Bl_large" = 0.05,
#"b2Bl_large" = 1.45)

alb <- run_first(list(runID = 0))
obs <- invert_model(inits) + generate.noise()

## Set initial conditions
inits <- c("orient_factor" = 0.15)
prior_def <- list(orient_factor = list("unif", list(-0.5, 0.5)))

init_function <- function() {
  out <- sapply(prior_def, 
                function(x) do.call(get(paste0("r", x[[1]])),
                                    c(1, x[[2]])))
  return(out)
}



priors <- priors_sun$means
#prior_function <- function(params) {
#  early_hardwood <- dmvnorm(c(params[1:6]), 
#                            mean = priors$M[,'temperate.Early_Hardwood'], sigma = priors$Sigma[,,'temperate.Early_Hardwood'], log = TRUE)
#  mid_hardwood <- dmvnorm(c(params[7:12]), 
#                          mean = priors$M[,'temperate.North_Mid_Hardwood'], sigma = priors$Sigma[,,'temperate.North_Mid_Hardwood'], log = TRUE)
#  orient_factor_early <- dunif(params[13], -0.5, 0.5)
# orient_factor_mid <- dunif(params[14], -0.5, 0.5)
  
  #  mid_hardwood <- dmvnorm(params[6:10], params[12]), mean = ..., log = TRUE)
  #...
  #return(early_hardwood + mid_hardwood + late_hardwood)
  #return(early_hardwood)
  #return(early_hardwood + mid_hardwood)
#  return(sum(early_hardwood, mid_hardwood, orient_factor_early, orient_factor_mid))
#}
#prior_function(rep(0,14))

prior_function <- function(params) {
  late_hardwood <- dmvnorm(c(params[1:6]),
                           mean = priors$M[,'temperate.Late_Hardwood'], sigma = priors$Sigma[,,'temperate.Late_Hardwood'], log = TRUE)
  orient_factor_late <- dunif(params[7], -0.5, 0.5)
  return(sum(late_hardwood+orient_factor_late,na.rm=T))
}
prior_function(rep(0,7))                         
                           
                           
param.mins <- c(orient_factor = -0.5)
param.maxs <- c(orient_factor = 0.5)

invert.options <- list(model = invert_model, 
                       run_first = run_first,
                       nchains = nchains,
                       inits.function = init_function,
                       prior.function = prior_function(rep(0,7)),
                       ngibbs.max = ngibbs.max,
                       ngibbs.min = ngibbs.min,
                       ngibbs.step = ngibbs.step,
                       param.mins = param.mins,
                       param.maxs = param.maxs,
                       adapt = 100,
                       adj_min = 0.1,
                       target = 0.234)
invert.options$do.lsq <- FALSE # TRUE/FALSE

runtag <- paste(paste0("dbh", dbh), pft, sep = ".")
fname <- paste("samples", runtag, "rds", sep = ".")
fname_prog <- paste("prog_samples", runtag, "rds", sep = ".")
logfile <- paste("output_prior", dttag, "log", sep = ".")
samples <- invert.auto(observed = obs, 
                       invert.options = invert.options,
                       parallel = TRUE,
                       parallel.output = paste0(main_out,"/",logfile),
                       save.samples = paste0(main_out,"/",fname_prog))



