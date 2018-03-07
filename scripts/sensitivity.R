library(redr)
library(hdf5r)
library(here)

source("config.R")

data(fft_plots)

n_sim <- 500

load("priors/mvtraits_priors.RData")
obs_LAI_dat <- read_csv(here("other_site_data/NASA_FFT_LAI_FPAR_Data.csv"))
aviris_specfile <- H5File$new("aviris/aviris.h5")

run_simulation <- function(fft_plot,
                           meas_year,
                           aviris_specfile = here("aviris/aviris.h5")) {

  obs_LAI <- obs_LAI_dat %>%
    dplyr::filter(Site_Plot == fft_plot, Collection_Year == meas_year) %>%
    dplyr::pull(LAI_Alpha_GammaE_CLX_10_60_m2_m2) %>%
    mean(na.rm = TRUE)

  aviris_h5 <- hdf5r::H5File$new(aviris_specfile)
  on.exit(aviris_h5$close_all())

  ## Aviris spectra
  aviris_refl <- filter_specfile(
    aviris_h5,
    iPLOT == !!fft_plot
  ) / 10000

  setup <- setup_fft(fft_plot, clobber = TRUE)

  edr_result <- run_simulation_edr(setup, means, covars, n_sim = 100)
  sail_result <- run_simulation_sail(setup, means, covars)

  plot(0, 0, type = "n", xlim = c(400, 2500), ylim = c(0, 0.8),
       main = paste(fft_plot, meas_year),
       xlab = "Wavelength (nm)",
       ylab = "Reflectance (0-1)")
  spec_ci(sail_result$refl_store, col = "green")
  spec_ci(edr_result$refl_store, col = "purple")
  if (all(dim(aviris_refl) > 0)) {
    matplot(aviris_refl, lty = "solid", col = "black", add = TRUE)
  }

}


#pdf("ed_sail_sims.pdf")
walk2(fft_plots$plot_code, fft_plots$meas_year, run_simulation)
#dev.off()
