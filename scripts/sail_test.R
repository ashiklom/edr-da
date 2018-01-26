library(tidyverse)
library(redr)
library(PEcAnRTM)
library(hdf5r)
library(here)

source("config.R")

fft_plots <- c("BH02","NC22","BI01","BI02","BI03","GR08", "OF04", "OF05", "PB03", "PB09","PB10",
               "PB13","SF01")
meas_year <- 2008

n_sim <- 500

load("priors/mvtraits_priors.RData")
obs_LAI_dat <- readr::read_csv(here::here("other_site_data/NASA_FFT_LAI_FPAR_Data.csv"))
aviris_specfile <- hdf5r::H5File$new("aviris/aviris.h5")

sail_simulation <- function(fft_plot,
                            meas_year = 2008,
                            aviris_specfile = system.file("aviris/aviris.h5",
                                                          package = "redr")) {
  obs_LAI <- obs_LAI_dat %>%
    dplyr::filter(Site_Plot == fft_plot, Collection_Year == meas_year) %>%
    dplyr::pull(LAI_Alpha_GammaE_CLX_10_60_m2_m2) %>%
    mean(na.rm = TRUE)

  aviris_h5 <- hdf5r::H5File$new(aviris_specfile)
  on.exit(aviris_h5$close_all())

  ## Aviris spectra
  aviris_refl <- filter_specfile(
    aviris_h5,
    iPLOT == fft_plot
  ) / 10000

  plot(0, 0, type = "n", xlim = c(400, 2500), ylim = c(0, 0.8),
       main = paste(fft_plot, meas_year),
       xlab = "Wavelength (nm)",
       ylab = "Reflectance (0-1)")
  spec_ci(refl_store, col = "green")
  if (all(dim(aviris_refl) > 0)) {
    matplot(aviris_refl, lty = "solid", col = "black", add = TRUE)
  }
}

pdf("sail_sims.pdf")
walk(fft_plots, possibly(sail_simulation, NULL))
dev.off()
