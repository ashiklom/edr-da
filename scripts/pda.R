# Parameter data assimilation
# Author: Alexey Shiklomanov
library(redr)
library(hdf5r)
library(here)

source("config.R")

fft_plots <- c("BH02","NC22","BI01","BI02","BI03","GR08", "OF04", "OF05", "PB03", "PB09","PB10",
               "PB13","SF01")
load("priors/mvtraits_priors.RData")

n_sim <- 500

meas_year <- 2008
fft_plot <- fft_plots[1]
aviris_specfile <- here("aviris/aviris.h5")

run_pda <- function(fft_plot,
                    meas_year = 2008,
                    aviris_specfile = here("aviris/aviris.h5")) {

  aviris_h5 <- hdf5r::H5File$new(aviris_specfile)
  ## Aviris spectra
  aviris_refl <- filter_specfile(
    aviris_h5,
    iPLOT == fft_plot
  ) / 10000
  aviris_h5$close_all()
  aviris_interp_wl <- round(wavelengths(aviris_refl))
  aviris_keep <- aviris_interp_wl > 400 &
    rowSums(is.na(aviris_refl)) == 0 &
    aviris_interp_wl <= 1300   # Use only VIS and NIR spectra
  observed <- aviris_refl[aviris_keep,]
  aviris_ind <- aviris_interp_wl[aviris_keep] - 399
  observation_operator <- function(x) x[aviris_ind]
  
  setup <- setup_fft(fft_plot, meas_year = 2008, clobber = FALSE)

  custom_settings <- list(
    #init = list(iterations = 25),
    #loop = list(iterations = 25),
    #other = list(max_iter = 50)
  )

  samples <- run_pda_edr(setup, observed, observation_operator,
                         means, covars, custom_settings = custom_settings)

  results_dir <- here("pda_results")
  dir.create(results_dir, showWarnings = FALSE)
  fname <- paste("PDA", "EDR", fft_plot, meas_year, "rds", sep = ".")
  saveRDS(samples, file.path(results_dir, fname))
}

walk(fft_plots, possibly(run_pda, NULL))
