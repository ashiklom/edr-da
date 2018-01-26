#' Run SAIL simulation
#'
#' Returns the following:
#' - `refl_store` -- Stored reflectance, wavelength x simulation
#' - `param_store` -- Drawn parameters, param x simulation
#' - `ed_LAI_total` -- ED LAI summed across cohorts, and used as SAIL input
#' - `ed_LAI_pft` -- ED LAI for each cohort
#' - `ed_pft_co` -- PFT designation of each cohort
#'
#' @param setup FFT setup object (see [setup_fft()])
#' @param n_sim Number of simulations
#' @inheritParams sample_params_prior
#' @return List of results (see description)
run_simulation_sail <- function(setup, prospect_means, prospect_covar,
                                n_sim = 500) {
  edr_setup <- setup_edr(setup$prefix, edr_exe_path = edr_exe_path)
  stopifnot(edr_setup)
  ed_LAI_pft <- get_edvar(setup$prefix, "LAI_CO")
  ed_LAI_total <- sum(ed_LAI_pft)
  ed_pft_table <- get_pfts(setup$css)
  ed_pft_co <- get_edvar(setup$prefix, "PFT")
  ed_pft_top <- ed_pft_table %>%
    filter(pft_num == ed_pft_co[1]) %>%
    pull(pft_name)
  test_sample <- sample_params_sail(ed_pft_top, ed_LAI_total, means, covars)
  param_store <- matrix(0, 15, n_sim)
  rownames(param_store) <- names(test_sample)
  refl_store <- matrix(0, 2101, n_sim)
  pb <- txtProgressBar()
  for (i in seq_len(n_sim)) {
    setTxtProgressBar(pb, i / seq_len(n_sim))
    param_store[, i] <- sample_params_sail(ed_pft_top, ed_LAI_total, means, covars)
    refl_store[, i] <- PEcAnRTM::pro4sail(param_store[, i])[, 1]
  }
  close(pb)
  list(
    refl_store = refl_store,
    param_store = param_store,
    ed_LAI_total = ed_LAI_total,
    ed_LAI_pft = ed_LAI_pft,
    ed_pft_co = ed_pft_co
  )
}

