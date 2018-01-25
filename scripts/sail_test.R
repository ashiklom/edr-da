library(tidyverse)
library(redr)
library(PEcAnRTM)
library(hdf5r)
library(here)

source("config.R")
source("common.R")

fft_plots <- c("BH02","NC22","BI01","BI02","BI03","GR08", "OF04", "OF05", "PB03", "PB09","PB10",
               "PB13","SF01")
meas_year <- 2008

n_sim <- 500

load("priors/mvtraits_priors.RData")
obs_LAI_dat <- readr::read_csv(here::here("other_site_data/NASA_FFT_LAI_FPAR_Data.csv"))
aviris_specfile <- hdf5r::H5File$new("aviris/aviris.h5")

draw_params <- function(pft, lai) {
  prosp_params <- -1
  while (any(prosp_params < 0)) {
    prosp_params <- mvtnorm::rmvnorm(1, means[pft, ], covars[,,pft])
  }
  sail_params <- c(
    N = prosp_params[1],
    Cab = prosp_params[2],
    Car = prosp_params[3],
    Cbrown = 0,
    Cw = prosp_params[4],
    Cm = prosp_params[5],
    LIDFa = runif(1, -0.9, 0.9),   # Default: -0.35
    LIDFb = runif(1, -0.9, 0.9),  # Default: -0.15
    TypeLIDF = 1,
    LAI = lai,
    q = runif(1, 0.01, 0.1),     # Hot spot, default = 0.01
    tts = 30,     # Solar zenith
    tto = 0,      # Observer zenith
    psi = 0,      # Sun-sensor azimuth
    psoil = 0.5   # Soil moisture
  )
  sail_params
}

sail_simulation <- function(fft_plot) {
  print(fft_plot)
  obs_LAI <- obs_LAI_dat %>%
    dplyr::filter(Site_Plot == fft_plot, Collection_Year == meas_year) %>%
    dplyr::pull(LAI_Alpha_GammaE_CLX_10_60_m2_m2) %>%
    mean(na.rm = TRUE)

  ## Aviris spectra
  aviris_refl <- filter_specfile(
    aviris_specfile,
    iPLOT == fft_plot
  ) / 10000

  setup <- setup_fft(fft_plot, meas_year, clobber = FALSE)
  edr_setup <- redr::setup_edr(setup$prefix, edr_exe_path = edr_exe_path)

  ed_LAI_pft <- get_edvar(setup$prefix, "LAI_CO")
  ed_LAI_total <- sum(ed_LAI_pft)

  ed_pft_table <- get_pfts(setup$css)
  ed_pft_co <- get_edvar(setup$prefix, "PFT")
  ed_pft_top <- ed_pft_table %>%
    filter(pft_num == ed_pft_co[1]) %>%
    pull(pft_name)

  param_store <- matrix(0, 15, n_sim)
  refl_store <- matrix(0, 2101, n_sim)

  for (i in seq_len(n_sim)) {
    param_store[, i] <- draw_params(ed_pft_top, ed_LAI_total)
    refl_store[, i] <- PEcAnRTM::pro4sail(param_store[, i])[, 1]
  }

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
