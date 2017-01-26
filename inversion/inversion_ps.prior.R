#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?

## Load functions
source("common.R")
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

## Link to ed_2.1-opt executable
#file.copy(from = edr_exe_path,
#          to = file.path(main_out, 'ed_2.1-opt'), 
#          overwrite = TRUE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Setup PROSPECT
#prospect.param <- c('N' = 1.4,
#                    'Cab' = 40,
#                    'Car' = 10,
#                    'Cw' = 0.01,
#                    'Cm' = 0.01)

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
    clumping_factor <- param[2]
    sla <- param[3]
    #b1Bl_large <- param[4]
    #b2Bl_large <- param[5]

    trait.values <- list()
    trait.values[[pft]] <- list(orient_factor = orient_factor,
                                clumping_factor = clumping_factor,
                                SLA = sla)
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
                  trait.values = list(),
                  output.path = outdir)
    return(albedo)
}


invert_model <- function(param, runID = 0) {

    outdir <- outdir_path(runID)
    paths_run <- list(ed2in = NA, history = outdir)

    trait.values <- get_trait_values(param) 

    albedo <- EDR.prospect(spectra_list = spectra_list,
                           trait.values = trait.values,
                           paths = paths_run,
                           par.wl = par.wl,
                           nir.wl = nir.wl,
                           datetime = datetime,
                           edr.exe.name = "ed_2.1",
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
inits <- c("orient_factor" = 0.25,
           "clumping_factor" = 0.75,
           "sla" = 18.85)
           #"b1Bl_large" = 0.05,
           #"b2Bl_large" = 1.45)

alb <- run_first(list(runID = 0))
obs <- invert_model(inits) + generate.noise()

## Set initial conditions
inits <- c("orient_factor" = 0.15,
           "clumping_factor" = 0.9,
           "sla" = 40)
            #"b1Bl_large" = 0.05,
            #"b2Bl_large" = 1.45)

prior_def <- list(orient_factor = list("unif", list(-0.5, 0.5)),
                  clumping_factor = list("unif", list(0, 1)),
                  sla = list("norm", list(18.904, 0.108)))
                  #b1Bl_large = list("lnorm", list(log(0.05), 0.1)),
                  #b2Bl_large = list("lnorm", list(log(1.45), 0.025)))

prior <- function(params) {
    out <- 0
    for (p in seq_along(prior_def)) {
        func <- get(paste0("d", prior_def[[p]][[1]]))
        out <- out + do.call(func, c(unname(params[p]), 
                                     prior_def[[p]][[2]],
                                     log = TRUE))
    }
    return(out)
}


init_function <- function() {
    out <- sapply(prior_def, 
                  function(x) do.call(get(paste0("r", x[[1]])),
                                      c(1, x[[2]])))
    return(out)
}

prior(inits)
prior(init_function())

param.mins <- c(orient_factor = -0.5,
                clumping_factor = 0,
                sla = 0)
                #b1Bl_large = 0,
                #b2Bl_large = 0)

param.maxs <- c(orient_factor = 0.5,
                clumping_factor = 1,
                sla = Inf)
                #b1Bl_large = Inf,
                #b2Bl_large = 1.488) # Hangs if any higher


invert.options <- list(model = invert_model, 
                       run_first = run_first,
                       nchains = nchains,
                       inits.function = init_function,
                       prior.function = prior,
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