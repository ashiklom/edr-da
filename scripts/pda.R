# Parameter data assimilation
# Author: Alexey Shiklomanov
library(redr)
library(hdf5r)
library(here)
library(PEcAn.logger)

logger.setLevel("INFO")   # This prevents the annoying "history" file debug messages

source("config.R")

data(fft_plots)
fft_plots$plot_code <- as.character(fft_plots$plot_code)

load("priors/mvtraits_priors.RData")

cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) == 0) {
  ind <- 1
} else {
  ind <- as.numeric(cmd_args)
}

stopifnot(
  all(!is.na(ind)),
  all(ind <= nrow(fft_plots))
)

fft_plots_sub <- fft_plots[ind, ]
message(
  "Running FFT plot ", fft_plots_sub$plot_code,
  " for year ", fft_plots_sub$meas_year,
  " (index ", ind, ")"
)

run_pda <- function(fft_plot,
                    meas_year,
                    aviris_specfile = here("aviris/aviris.h5")) {

  aviris_h5 <- hdf5r::H5File$new(aviris_specfile)
  ## Aviris spectra
  aviris_refl <- filter_specfile(
    aviris_h5,
    iPLOT == !!fft_plot
  ) / 10000
  aviris_h5$close_all()
  if (any(dim(aviris_refl) == 0)) {
    stop("Missing AVIRIS spectra. Not running PDA.")
  }
  aviris_interp_wl <- round(wavelengths(aviris_refl))
  aviris_keep <- aviris_interp_wl > 400 &
    rowSums(is.na(aviris_refl)) == 0 &
    aviris_interp_wl <= 1300   # Use only VIS and NIR spectra
  observed <- aviris_refl[aviris_keep,]
  aviris_ind <- aviris_interp_wl[aviris_keep] - 399
  observation_operator <- function(x) x[aviris_ind]
  
  setup <- setup_fft(fft_plot, clobber = TRUE)

  results_dir <- here("pda_results")
  progress_dir <- here("pda_progress")
  dir.create(results_dir, showWarnings = FALSE)
  dir.create(progress_dir, showWarnings = FALSE)
  fname <- paste("PDA", "EDR", fft_plot, meas_year, "rds", sep = ".")

  custom_settings <- list(
    init = list(iterations = 500),
    loop = list(iterations = 500),
    other = list(save_progress = file.path(progress_dir, fname))
  )

  samples <- run_pda_edr(setup, observed, observation_operator,
                         means, covars, custom_settings = custom_settings)

  saveRDS(samples, file.path(results_dir, fname))
}

walk2(fft_plots_sub$plot_code, fft_plots_sub$meas_year, run_pda)
