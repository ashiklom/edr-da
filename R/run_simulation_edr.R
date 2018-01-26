#' Run EDR simulation
#'
#' @inheritParams run_simulation_sail
#' @param num_of_edr_params Number of EDR params (default = 3)
#' @return List of outputs. See [run_simulation_sail()].
#' @export
run_simulation_edr <- function(setup,
                               prospect_means,
                               prospect_covar,
                               n_sim = 500,
                               num_of_edr_params = 3) {
  edr_setup <- setup_edr(setup$prefix, edr_exe_path = edr_exe_path)
  stopifnot(edr_setup)
  edr_run_pfts <- get_pfts(setup$css)
  PEcAn.logger::logger.info(paste0("Running with :  ", edr_run_pfts$pft_name))
  edr_run_pfts_length <- dim(edr_run_pfts)[1]
  pft_vec_length <- 5 + num_of_edr_params
  pft_ends <- seq(pft_vec_length, edr_run_pfts_length * pft_vec_length, pft_vec_length)
  test_params <- sample_params_edr(edr_run_pfts$pft_name, prospect_means, prospect_covar)
  testrun <- test_params %>%
    vec2list(datetime = datetime, par.wl = 400:2499, nir.wl = 2500) %>%
    run_edr(setup$prefix, edr_args = .)
  ed_LAI_pft <- get_edvar(setup$prefix, "LAI_CO")
  ed_LAI_total <- sum(ed_LAI_pft)
  ed_pft_co <- get_edvar(setup$prefix, "PFT")
  PEcAn.logger::logger.info(paste0("Initial LAI:  ", LAI))
  refl_store <- matrix(numeric(), 2101, n_sim)
  param_store <- matrix(numeric(), length(test_params), n_sim)
  pb <- txtProgressBar()
  for (i in seq_len(n_sim)) {
    param_vector <- sample_params_edr(
      edr_run_pfts$pft_name,
      prospect_means,
      prospect_covar
    )
    param_store[i, ] <- param_vector
    refl_store[i, ] <- edr_model(param_vector)
    setTxtProgressBar(pb, i / seq_len(n_sim))
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

#' Subset parameters for ith PFT
#'
#' @param i PFT index
#' @param params Full vector of parameter values
#' @param pft_ends Vector of PFT index ends
param_sub <- function(i, params, pft_ends) {
  param_seq <- (pft_ends[i] - 7):pft_ends[i]
  params[param_seq]
}

#' Convert EDR parameter vector to EDR input list
#'
#' @param params Vector of EDR parameters (see [sample_params_edr()])
#' @export
vec2list <- function(params, pft_vec, ...) {  # lots of hard-coded params here
  spectra_list <- list()
  ed_list <- list()
  for (i in seq_along(pft_vec)) {
    pft <- as.character(pft_vec[i])
    param_sub <- param_sub(i, params)
    spectra_list[[pft]] <- PEcAnRTM::prospect(param_sub[1:5], 5)
    ed_list[[pft]] <- list(
      SLA = param_sub[6],
      clumping_factor = param_sub[7],
      orient_factor = param_sub[8]
    )
  }
  list(spectra_list = spectra_list, trait.values = ed_list, ...)
}

