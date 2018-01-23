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
## Setup parameter prior info
load('priors/mvtraits_priors.RData')

# EDR params
num_of_edr_params <- 3    # length of edr param vector, e.g. 3 for SLA, of, and cf
#cf <- "runif(1, 0, 1)"  # update to allow to change dist options here
#of <- "runif(1, -0.5, 0.5)"
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup plots to simulate and other options

fft_plots <- c("BH02","NC22", "GR08") # all
#fft_plots <- c("BH02", "NC22")
#fft_plots <- c("NC22")
years <- 2008   # update so that we can select specific years to run
clobber <- TRUE

edr_iterations <- 5
#--------------------------------------------------------------------------------------------------#

# end of user input
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup helper functions
pft_lookup <- list(pft_num = c(6,7,8,9,10,11), 
                   pft_names = c("temperate.Northern_Pine","temperate.Southern_Pine",
                               "temperate.Late_Conifer","temperate.Early_Hardwood",
                               "temperate.North_Mid_Hardwood","temperate.Late_Hardwood"))
get_pfts <- function(css_file_path, delim = ' ', reference = pft_lookup) {
  css_data <- read.delim(css_file_path, header = T, sep = delim)
  css_pfts <- as.vector(unique(css_data["pft"]))
  index <- which(pft_lookup$pft_num %in% css_pfts$pft)
  return(data.frame(pft_num = pft_lookup$pft_num[index], pft_name = pft_lookup$pft_names[index]))
}
#get_pfts(css_in)

param_sub <- function(i, params) {
  param_seq <- (pft_ends[i] - 7):pft_ends[i]
  param_sub <- params[param_seq]
  return(param_sub)
}

inits_function <- function() {  # lots of hard coded param assumptions here.....
  samples <- numeric()
  for (i in seq_along(1:edr_run_pfts_length)) {
    #pft <- names(pft_end)[i]
    pft <- as.character(edr_run_pfts$pft_name[i])
    
    # PROSPECT prior
    #samples <- c(samples, mvtnorm::rmvnorm(1, means[pft,], covars[pft,,]))
    draw <- -1
    while (any(draw < 0)) {
      draw <- mvtnorm::rmvnorm(1, means[pft,], covars[,,pft])
    }
    samples <- c(samples, draw)
    # ED priors
    samples <- c(samples, runif(1, 0, 1), runif(1, -0.5, 0.5)) # clumping and orient factor
  }
  name_vec <- rep(c("N","Cab","Car","Cw","Cm","SLA","clumping_factor","orient_factor"),edr_run_pfts_length)
  names(samples) <- name_vec
  samples <- c(samples, residual = rlnorm(1, log(0.001), 2.5))
  return(samples)
}

vec2list <- function(params, ...) {  # lots of hard-coded params here
  spectra_list <- list()
  ed_list <- list()
  for (i in seq_along(1:edr_run_pfts_length)) {
    pft <- as.character(edr_run_pfts$pft_name[i])
    param_sub <- param_sub(i, params)
    spectra_list[[pft]] <- PEcAnRTM::prospect(param_sub[1:5], 5, TRUE)
    ed_list[[pft]] <- list(SLA = param_sub[6],
                           clumping_factor = param_sub[7],
                           orient_factor = param_sub[8])
  }
  outlist <- list(spectra_list = spectra_list, trait.values = ed_list, ...)
  return(outlist)
}

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
# run sites
for (i in seq_along(1:length(fft_plots))) {
  PEcAn.logger::logger.info(paste0("Current Plot:  ",fft_plots[i]))
  
  ed_inputs_dir <- normalizePath(Sys.glob(file.path('ed-inputs/sites', paste0(fft_plots[i],"_*"))) )
  ed_inputs <- list.files(ed_inputs_dir)
  tmp <- strsplit(ed_inputs_dir, "ed-inputs")[[1]]
  edr_dir <- tmp[1]

  site_lat <- as.numeric(gsub('^.*lat\\s*|\\s*lon.*$', '', ed_inputs[1]))
  site_lon <- as.numeric(gsub('^.*lon\\s*|\\s*.css.*$', '', ed_inputs[1]))
  meas_year <- as.numeric(gsub('^.*FFT.\\s*|\\s*.lat.*$', '', ed_inputs[1]))
  
  # link to css, pss, site files
  css_in <- list.files(ed_inputs_dir,pattern = "*.css", recursive = TRUE, full.names = T)
  pss_in <- list.files(ed_inputs_dir,pattern = "*.pss", recursive = TRUE, full.names = T)
  site_in <- list.files(ed_inputs_dir,pattern = "*.site", recursive = TRUE, full.names = T)
  
  ### create output prefix
  prefix <- paste("EDR_sim_output", fft_plots[i], meas_year, sep = "_")
  if (clobber) PEcAn.logger::logger.info("*** Removing previous simulation results ***")
  if (file.exists(paste0(edr_dir,prefix)) && clobber) unlink(paste0(edr_dir,prefix), recursive = T)
  PEcAn.logger::logger.info(paste0("Running simulation in dir: ",prefix))
  
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
  
  
  ## Setup EDR
  edr_setup <- setup_edr(prefix, edr_exe_path = edr_exe_path)
  
  ## Setup EDR PFTs
  edr_run_pfts <- get_pfts(css_in)
  edr_run_pfts_length <- dim(edr_run_pfts)[1]
  pft_vec_length <- 5+edr_run_pfts_length
  pft_ends <- seq(pft_vec_length,edr_run_pfts_length*pft_vec_length, pft_vec_length)

  testrun <- inits_function() %>%
    vec2list(datetime = datetime, par.wl = 400:2499, nir.wl = 2500) %>%
    run_edr(prefix, edr_args = .)
  
  
  ## Run EDR
  #edr_iterations <- 5
  edr_albedo <- array(data=NA, c(edr_iterations,2101))
  for (j in seq_along(1:dim(edr_albedo)[1])) {
    print(paste0(" Simulation: ",j))
    edr_albedo[j,] <- model(inits_function())
  }
  
  edr_albedo_mean <- colMeans(edr_albedo, na.rm=T)
  edr_albedo_quant <- apply(edr_albedo,2,quantile,na.rm=T,probs=c(0,0.05,0.95,1))
  
  
  ## Generate comparison figure
  if (generate_summary_figs) {
    waves <- seq(400,2500,1)
    aviris.data <- read.csv('/data/sserbin/Modeling/edr-da/aviris/NASA_FFT/allPlotsAggregated_NIT_A_SpecJoined.csv',
                            header = T,sep=",")
    aviris.spectra <- as.matrix(droplevels(aviris.data[,paste0("band_",seq(1,224,1))]))
    spectra <- aviris.spectra[which(aviris.data$iPLOT==fft_plots[i]),]
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
  
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF