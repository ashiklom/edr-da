####################################################################################################
#          
#   Generate forward EDR simulations to compare against AVIRIS observations
#   Author: Shawn P. Serbin
#
#
#
#
#
#    --- Last updated:  08.29.2017 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
source('config.R')
PEcAn.logger::logger.setLevel("INFO")

## Load libraries
library(redr)
library(magrittr)

## Setup options
plot_albedo <- TRUE #TRUE/FALSE
generate_summary_figs <- TRUE #TRUE/FALSE
hidden <- FALSE  #TRUE/FALSE  # output to a hidden (e.g. .forward_sim_output) folder?
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (hidden) {
  prefix <- '.EDR_sim_output'
} else {
  prefix <- paste("EDR_sim_output", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
}
PEcAn.logger::logger.info(paste0("Running simulation in dir: ",prefix))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Run ED2
fft_plot <- "NC22"
css_in <- '/data/sserbin/Modeling/edr-da/ed-inputs/sites/US-WCr/NASA_FFT/US-WCr.Inv2010.NC22.lat45.5lon-90.5.css'
pss_in <- '/data/sserbin/Modeling/edr-da/ed-inputs/sites/US-WCr/NASA_FFT/US-WCr.Inv2010.NC22.lat45.5lon-90.5.pss'
site_in <- '/data/sserbin/Modeling/edr-da/ed-inputs/sites/US-WCr/NASA_FFT/US-WCr.Inv2010.NC22.lat45.5lon-90.5.site'
ed2in_changes <- list(IMONTHA = 07, IDATEA = 15, IYEARA = 2006,
                      IMONTHZ = 08, IDATEZ = 31, IYEARZ = 2006)
datetime <- ISOdate(2006, 08, 26, 16, 00, 00)
genrun <- generate_run(prefix = prefix,
                       site_lat = site_lat,
                       site_lon = site_lon,
                       site_df = site_in,
                       pss_df = pss_in,
                       css_df = css_in,
                       common_inputs_dir = common_inputs_dir,
                       site_met_dir = site_met_dir,
                       ed_exe_path = ed_exe_path,
                       ed2in_changes = ed2in_changes,
                       RMDIR = TRUE)

message('Running ED...')
runed <- run_ed(prefix)
tail(runed)
message('Done!')

edr_setup <- setup_edr(prefix, edr_exe_path = edr_exe_path)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Setup priors
load('priors/sunlit_meanjp.RData')

# Params vector
# Temperate.Early_Hardwood: 1:8
#   1:5 -- PROSPECT params
#   6 -- sla
#   7 -- clumping_factor
#   8 -- orient_factor
# Temperate.North_Mid_Hardwood: 9:16
# Temperate.Late_Hardwood: 17:24
pft_end <- c(temperate.Early_Hardwood = 8,
             temperate.North_Mid_Hardwood = 16,
             temperate.Late_Hardwood = 24)

param_sub <- function(i, params) {
  param_seq <- (pft_end[i] - 7):pft_end[i]
  param_sub <- params[param_seq]
  return(param_sub)
}

inits_function <- function() {
  samples <- numeric()
  for (i in seq_along(pft_end)) {
    pft <- names(pft_end)[i]
    # PROSPECT prior
    samples <- c(samples, mvtnorm::rmvnorm(1, means[pft,], covars[pft,,]))
    # ED priors
    #samples <- c(samples, runif(1, 0, 1), runif(1, -0.5, 0.5)) # clumping and orient factor
    samples <- c(samples, runif(1, 0.8, 0.87), runif(1, 0.15, 0.2)) # clumping and orient factor
  }
  
  name_vec <- rep(c("N","Cab","Car","Cw","Cm","leaf_mass_per_area","clumping_factor","orient_factor"),length(pft_end))
  samples[which(name_vec=="leaf_mass_per_area")] <- 1/samples[which(name_vec=="leaf_mass_per_area")]*1000 # convert to SLA
  names(samples) <- rep(c('N', 'Cab', 'Car', 'Cw', 'Cm', 'SLA', 'clumping_factor', 'orient_factor'),3)
  return(samples)
}

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Setup EDR
vec2list <- function(params, ...) {
  spectra_list <- list()
  ed_list <- list()
  for (i in seq_along(pft_end)) {
    pft <- names(pft_end)[i]
    param_sub <- param_sub(i, params)
    spectra_list[[pft]] <- PEcAnRTM::prospect(param_sub[1:5], 5, TRUE)
    ed_list[[pft]] <- list(SLA = param_sub[6],
                           clumping_factor = param_sub[7],
                           orient_factor = param_sub[8])
  }
  outlist <- list(spectra_list = spectra_list, trait.values = ed_list, ...)
  return(outlist)
}

testrun <- inits_function() %>%
  vec2list(datetime = datetime, par.wl = 400:2499, nir.wl = 2500) %>%
  run_edr(prefix, edr_args = .)

model <- function(params) {
  edr_dir <- 'edr'
  args_list <- vec2list(params,
                        paths = list(ed2in = NA, history = file.path(prefix, 'outputs')),
                        par.wl = 400:2499,
                        nir.wl = 2500,
                        edr.exe.name = 'edr 2> /dev/null',
                        datetime = datetime,
                        change.history.time = FALSE)
  albedo <- run_edr(prefix, args_list, edr_dir)
  
  if (plot_albedo) {
    # Create quick figure
    waves <- seq(400,2500,1)
    png(file.path(prefix,edr_dir,'simulated_albedo.png'),width=4900, height =2700,res=400)
    par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
    plot(waves,unlist(albedo)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",
         ylab="Reflectance (%)",
         cex.axis=1.5, cex.lab=1.7,col="black")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey80")
    lines(waves,unlist(albedo)*100,lwd=3, col="black")
    dev.off()
  }
  
  return(albedo)
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Run EDR
# Test that the model works
run_model <- model(inits_function())
head(run_model)



#spectra_list <- list(temperate.Early_Hardwood = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.01), 4, TRUE),
#                     temperate.North_Mid_Hardwood = PEcAnRTM::prospect(c(1.6, 50, 0.01, 0.012), 4, TRUE),
#                     temperate.Late_Hardwood = PEcAnRTM::prospect(c(1.8, 60, 0.01, 0.012), 4, TRUE))
#trait.values <- list(temperate.Early_Hardwood =  list(orient_factor = 0.1, clumping_factor = 0.86, SLA=15),
#                     temperate.North_Mid_Hardwood = list(orient_factor = 0.1, clumping_factor = 0.86, SLA=7),
#                     temperate.Late_Hardwood = list(orient_factor = 0.1, clumping_factor = 0.86, SLA=15))
#edr_args <- list(spectra_list = spectra_list, trait.values = trait.values,par.wl = 400:2499,nir.wl = 2500,
#                 datetime = datetime)
#run_model <- run_edr(dir=prefix,edr_args=edr_args)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Generate comparison figure
aviris.data <- read.csv('/data/sserbin/Modeling/edr-da/aviris/NASA_FFT/allPlotsAggregated_NIT_A_SpecJoined.csv',header = T,sep=",")
aviris.spectra <- as.matrix(droplevels(aviris.data[,paste0("band_",seq(1,224,1))]))
spectra <- aviris.spectra[which(aviris.data$iPLOT==fft_plot),]
av.waves <- read.csv('/data/sserbin/Modeling/edr-da/aviris/NASA_FFT/aviris_c_wavelength.csv',header = T,sep=",")

waves <- seq(400,2500,1)
png(file.path(prefix,'edr','EDR_vs_AVIRIS_Albedo.png'),width=4900, height =2700,res=400)
par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
plot(waves,unlist(run_model)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
     cex.axis=1.5, cex.lab=1.7,col="black")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
       "grey80")
lines(waves,unlist(run_model)*100,lwd=3, col="green4")
lines(unlist(av.waves[,1]), spectra*0.01,lwd=2)
legend("topleft",legend=c("Modeled","Measured"), lwd=3,col=c("green4","black"))
dev.off()

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF