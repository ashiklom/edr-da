#--------------------------------------------------------------------------------------------------#
source("config.R")
source("common.R")
PEcAn.logger::logger.setLevel("INFO")

## Load libraries
library(redr)
library(ncdf4)
library(tidyverse)
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
### Setup field plots to simulate and other options
#fft_plots <- c("BH02","NC22", "GR08", "IDS05", "IDS10") # all
#fft_plots <- c("BH02","NC22","BI01","BI02","BI03","GR08", "OF04", "OF05", "PB03", "PB09","PB10",
#               "PB13","SF01")
#fft_plots <- c("BH02","OF04", "OF05", "PB03", "PB09","PB10","PB13","SF01")
# fft_plots <- c("PB09","PB10","PB13","SF01")
#fft_plots <- c("PB03") # plot missing from AVIRIS data
#fft_plots <- c("PB09","OF05")
#fft_plots <- c("PB03")

#fft_plots <- c("BH02", "NC22","BI01","BI02","BI03","GR08", "OF04", "OF05", "PB03", "PB09","PB10",
#               "PB13","SF01")
fft_plots <- c("BH02")

years <- 2008   # update so that we can select specific years to run...NOT YET USED
clobber <- TRUE

edr_iterations <- 5

## In-situ LAI data
obs_LAI <- read_csv(here::here('other_site_data/NASA_FFT_LAI_FPAR_Data.csv'))

## AVIRIS data
# move up reading in AVIRIS data to memory here
#--------------------------------------------------------------------------------------------------#

# end of user input
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

model <- function(params) {
  edr_dir <- 'edr'
  paths_list <- list(
    ed2in = NA,
    history = here::here(prefix, "outputs"),
    soil_reflect = soil_refl_file
  )
  args_list <- vec2list(params,
                        paths = paths_list,
                        par.wl = 400:800,
                        nir.wl = 801:2500,
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
for (i in seq_along(fft_plots)) {
  PEcAn.logger::logger.info(paste0("Current Plot:  ", fft_plots[i]))

  setup <- setup_fft(fft_plots[i], meas_year, clobber = FALSE)
  
  ## Setup EDR
  edr_setup <- setup_edr(setup$prefix, edr_exe_path = edr_exe_path)
  
  ## Setup EDR PFTs
  edr_run_pfts <- get_pfts(setup$css)
  PEcAn.logger::logger.info(paste0("Running with :  ", edr_run_pfts$pft_name))
  edr_run_pfts_length <- dim(edr_run_pfts)[1]
  pft_vec_length <- 5 + num_of_edr_params
  pft_ends <- seq(pft_vec_length, edr_run_pfts_length * pft_vec_length, pft_vec_length)

  testrun <- inits_function() %>%
    vec2list(datetime = datetime, par.wl = 400:2499, nir.wl = 2500) %>%
    run_edr(setup$prefix, edr_args = .)
  
  ## Grab initial LAI
  LAI <- sum(get_edvar(setup$prefix, "LAI_CO"))
  PEcAn.logger::logger.info(paste0("Initial LAI:  ", LAI))
  
  ## Run EDR
  #history_files <- list.files(file.path(prefix,"edr"),pattern="history-S-*.*h5")
  edr_albedo <- array(data = NA, c(edr_iterations, 2101))
  LAI <- rep(NA, edr_iterations)
  for (j in seq_len(dim(edr_albedo)[1])) {
    print(paste0(" Simulation: ", j))
    param_vector <- inits_function()
    edr_albedo[j,] <- model(param_vector)
    LAI[j] <- sum(get_edvar(prefix, "LAI_CO"))
  }
  
  ## create outputs
  # LAI 
  output_LAI <- data.frame(Iteration = seq(1,edr_iterations,1), LAI)
  names(output_LAI) <- c("Iteration", "ED2_LAI")
  write.csv(output_LAI, file = file.path(paste0(edr_dir,prefix),"edr","EDR_LAI.csv"), row.names = FALSE)
  # Albedo
  output_albedo <- data.frame(Iteration = seq(1,edr_iterations,1), edr_albedo)
  names(output_albedo) <- c("Iteration", paste0("Wave_",seq(400,2500,1)))
  #output_albedo[1:3, 1:10]
  write.csv(output_albedo, file = file.path(paste0(edr_dir,prefix),"edr","EDR_Albedo.csv"), row.names = FALSE)

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
    if (length(which(aviris.data$iPLOT==fft_plots[i]))<1 ) {
      png(file.path(prefix,'edr','EDR_vs_AVIRIS_Albedo.png'),width=4900, height =2700,res=400)
      par(mfrow=c(1,1), mar=c(4.3,6,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
      plot(waves,unlist(edr_albedo_mean)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
           cex.axis=1.5, cex.lab=1.7,col="black", main=fft_plots[i])
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
             "grey80")
      polygon(c(waves ,rev(waves)),c(edr_albedo_quant[3,]*100, rev(edr_albedo_quant[2,]*100)),
              col="#99CC99",border=NA)
      lines(waves,unlist(edr_albedo_mean)*100,lwd=3, col="green4")
      legend("topleft",legend=c("Modeled"), lwd=3,col=c("green4"))
      dev.off()
    } else {
      png(file.path(prefix,'edr','EDR_vs_AVIRIS_Albedo.png'),width=4900, height =2700,res=400)
      par(mfrow=c(1,1), mar=c(4.3,6,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
      plot(waves,unlist(edr_albedo_mean)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
           cex.axis=1.5, cex.lab=1.7,col="black", main=fft_plots[i])
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
             "grey80")
      polygon(c(waves ,rev(waves)),c(edr_albedo_quant[3,]*100, rev(edr_albedo_quant[2,]*100)),
              col="#99CC99",border=NA)
      lines(waves,unlist(edr_albedo_mean)*100,lwd=3, col="green4")
      if (is.null(dim(spectra)[1])) {
        lines(unlist(av.waves[,1]), as.vector(spectra)*0.01,lwd=2)
      } else {
        for (k in seq_along(1:dim(spectra)[1])) {
          lines(unlist(av.waves[,1]), spectra[k,]*0.01,lwd=2)
        }
      }
      legend("topleft",legend=c("Modeled","Measured"), lwd=3,col=c("green4","black"))
      dev.off()
    }

    
    
    # LAI histogram
    if (!any(obs_LAI$Site_Plot==fft_plots[i])==TRUE) {
      # need to check if LAI data availible for comparison first
      png(file.path(prefix,'edr','EDR_LAI_histogram.png'),width=4200, height =2700,res=400)
      EDR_LAI_range <- range(LAI, na.rm=TRUE)
      hist(LAI, xlim = c(EDR_LAI_range[1]-0.5, EDR_LAI_range[2]+0.5))
      dev.off()
    } else {
      plot_LAI <- obs_LAI[which(obs_LAI$Site_Plot==fft_plots[i]),] 
      plot_LAI_mn <- mean(plot_LAI$LAI_Alpha_GammaE_CLX_10_60_m2_m2, na.rm=TRUE)
      EDR_LAI_range <- range(LAI, na.rm=TRUE)
      if (plot_LAI_mn < EDR_LAI_range[1]) {
        hist_range <- c(plot_LAI_mn-0.5, EDR_LAI_range[2]+0.5)
      } else if (plot_LAI_mn > EDR_LAI_range[2]) {
        hist_range <- c(EDR_LAI_range[1]-0.5, plot_LAI_mn+0.5)
      } else {
        hist_range <- EDR_LAI_range
      }
      png(file.path(prefix,'edr','EDR_LAI_histogram.png'),width=4200, height =2700,res=400)
      hist(LAI, xlim = hist_range)
      abline(v=plot_LAI_mn,col="red", lwd=3)
      dev.off()
    }

  }
  
  # cleanup
  rm(i)
  
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF
