####################################################################################################
#
#   Generate forward EDR simulations to compare against AVIRIS observations
#   Author: Shawn P. Serbin
#
#
#
#
#
#    --- Last updated:  10.03.2017 By Shawn P. Serbin <sserbin@bnl.gov>
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
fft_plot <- "GR08"  # NC22, BH02
#css_in <- '/data/sserbin/Modeling/edr-da/ed-inputs/sites/US-WCr/NASA_FFT/US-WCr.Inv2010.NC22.lat45.5lon-90.5.css'
#pss_in <- '/data/sserbin/Modeling/edr-da/ed-inputs/sites/US-WCr/NASA_FFT/US-WCr.Inv2010.NC22.lat45.5lon-90.5.pss'
#site_in <- '/data/sserbin/Modeling/edr-da/ed-inputs/sites/US-WCr/NASA_FFT/US-WCr.Inv2010.NC22.lat45.5lon-90.5.site'

#css_in <- 'ed-inputs/sites/BH02_site_1-25665/FFT.2008.lat43.3724lon-89.9071.css'
#pss_in <- 'ed-inputs/sites/BH02_site_1-25665/FFT.2008.lat43.3724lon-89.9071.pss'
#site_in <- 'ed-inputs/sites/BH02_site_1-25665/FFT.2008.lat43.3724lon-89.9071.site'

css_in <- 'ed-inputs/sites/GR08_site_1-25681/FFT.2008.lat39.695lon-78.4661.css'
pss_in <- 'ed-inputs/sites/GR08_site_1-25681/FFT.2008.lat39.695lon-78.4661.pss'
site_in <- 'ed-inputs/sites/GR08_site_1-25681/FFT.2008.lat39.695lon-78.4661.site'


site_lat <- 39.695 #45.5 #43.3724
site_lon <- -78.4661 #-90.5 #-89.9071

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
#load('priors/sunlit_meanjp.RData') # old prior
load('priors/mvtraits_priors.RData')

# Params vector
# Temperate.Early_Hardwood: 1:8
#   1:5 -- PROSPECT params
#   6 -- sla
#   7 -- clumping_factor
#   8 -- orient_factor
# Temperate.North_Mid_Hardwood: 9:16
# Temperate.Late_Hardwood: 17:24
#pft_end <- c(temperate.North_Mid_Hardwood = 8,
#             temperate.Late_Hardwood = 16,
#             temperate.Late_Conifer = 24)
pft_end <- c(temperate.North_Mid_Hardwood = 8,
             temperate.Late_Hardwood = 16,
             temperate.Northern_Pine = 24)


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
    #samples <- c(samples, mvtnorm::rmvnorm(1, means[pft,], covars[pft,,]))
    draw <- -1
    while (any(draw < 0)) {
      draw <- mvtnorm::rmvnorm(1, means[pft,], covars[,,pft])
    }
    samples <- c(samples, draw)
    # ED priors
    samples <- c(samples, runif(1, 0, 1), runif(1, -0.5, 0.5)) # clumping and orient factor
    #samples <- c(samples, runif(1, 0.8, 0.87), runif(1, 0.01, 0.15)) # clumping and orient factor
    #samples <- c(samples, runif(1, 0.55, 0.98), runif(1, -0.5, 0.5)) # clumping and orient factor
  }
  name_vec <- rep(c("N","Cab","Car","Cw","Cm","SLA","clumping_factor","orient_factor"),length(pft_end))
  names(samples) <- name_vec
  samples <- c(samples, residual = rlnorm(1, log(0.001), 2.5))
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
# single run
#run_model <- model(inits_function())
#head(run_model)

# run X times
x <- 100
edr_albedo <- array(data=NA, c(x,2101))
for (i in seq_along(1:dim(edr_albedo)[1])) {
  print(paste0(" Simulation: ",i))
  edr_albedo[i,] <- model(inits_function())
  #tryCatch(edr_albedo[i,] <- model(inits_function()), finally = print("Hello"))
}

edr_albedo_mean <- colMeans(edr_albedo, na.rm=T)
edr_albedo_quant <- apply(edr_albedo,2,quantile,na.rm=T,probs=c(0,0.05,0.95,1))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Generate comparison figure
if (generate_summary_figs) {
  waves <- seq(400,2500,1)
  aviris.data <- read.csv('/data/sserbin/Modeling/edr-da/aviris/NASA_FFT/allPlotsAggregated_NIT_A_SpecJoined.csv',
                          header = T,sep=",")
  aviris.spectra <- as.matrix(droplevels(aviris.data[,paste0("band_",seq(1,224,1))]))
  spectra <- aviris.spectra[which(aviris.data$iPLOT==fft_plot),]
  av.waves <- read.csv('/data/sserbin/Modeling/edr-da/aviris/NASA_FFT/aviris_c_wavelength.csv',header = T,sep=",")
  #spectra_quant <- apply(spectra,2,quantile,na.rm=T,probs=c(0,1))
  #spectra_mean <- colMeans(spectra)
  png(file.path(prefix,'edr','EDR_vs_AVIRIS_Albedo.png'),width=4900, height =2700,res=400)
  par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
  plot(waves,unlist(edr_albedo_mean)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
       cex.axis=1.5, cex.lab=1.7,col="black")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
         "grey80")
  polygon(c(waves ,rev(waves)),c(edr_albedo_quant[3,]*100, rev(edr_albedo_quant[2,]*100)),
          col="#99CC99",border=NA)
  lines(waves,unlist(edr_albedo_mean)*100,lwd=3, col="green4")
  if (is.null(dim(spectra)[1])) {
    lines(unlist(av.waves[,1]), as.vector(spectra)*0.01,lwd=2)
  } else {
    for (i in seq_along(1:dim(spectra)[1])) {
      lines(unlist(av.waves[,1]), spectra[i,]*0.01,lwd=2)
    }
  }
  legend("topleft",legend=c("Modeled","Measured"), lwd=3,col=c("green4","black"))
  dev.off()
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF
