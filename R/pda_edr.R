#' Run EDR parameter data assimilation
#'
#' @inheritParams run_simulation_edr
#' @param observed Observation vector or matrix (see [PEcAnRTM::invert_bt()])
#' @param observation_operator Function for aligning EDR output with observation
#' @param custom_settings Custom settings list for [PEcAnRTM::invert_bt()]
#' @return Output of [PEcAnRTM::invert_bt()]
#' @export
run_pda_edr <- function(setup,
                        observed,
                        observation_operator,
                        prospect_means,
                        prospect_covar,
                        num_of_edr_params = 3,
                        custom_settings = list()) {
  stopifnot(is.function(observation_operator))
  edr_setup <- setup_edr(setup$prefix, edr_exe_path = edr_exe_path)
  edr_run_pfts <- get_pfts(setup$css)
  PEcAn.logger::logger.info("Running with :  ", edr_run_pfts$pft_name)
  edr_run_pfts_length <- dim(edr_run_pfts)[1]
  pft_vec_length <- 5 + num_of_edr_params
  pft_ends <- seq(pft_vec_length, edr_run_pfts_length * pft_vec_length, pft_vec_length)
  test_params <- sample_params_edr(edr_run_pfts$pft_name, prospect_means, prospect_covar)
  itest <- 1
  while (itest <= 20) {
    testrun <- try({
      test_params %>%
        vec2list(
          pft_vec = edr_run_pfts$pft_name,
          pft_ends = pft_ends,
          datetime = datetime,
          par.wl = 400:2499,
          nir.wl = 2500
        ) %>%
        run_edr(setup$prefix, edr_args = .)
    })
    if (class(testrun) == "try-error") {
      itest <- itest + 1
      message("Bad parameter draw for test run. Trying again with attempt #", itest)
    } else {
      break
    }
  }
  if (itest > 20) {
    stop("Could not run EDR after 20 attempts.")
  }
  ed_LAI_pft <- get_edvar(setup$prefix, "LAI_CO")
  ed_LAI_total <- sum(ed_LAI_pft)
  ed_pft_co <- get_edvar(setup$prefix, "PFT")

  invert_model <- function(params) {
    raw_out <- edr_model(
      params,
      setup$prefix,
      edr_run_pfts$pft_name,
      pft_ends
    )
    observation_operator(raw_out)
  }

  prior_function <- function(params) {
    prior <- numeric(1)
    for (i in seq_along(pft_ends)) {
      param_sub <- param_sub(i, params, pft_ends)
      pft <- edr_run_pfts$pft_name[i]
      # PROSPECT prior
      prior <- prior + mvtnorm::dmvnorm(param_sub[1:6], means[pft,], covars[,,pft], log = TRUE)
      # ED priors
      prior <- prior + dunif(param_sub[7], 0, 1, TRUE) + dunif(param_sub[8], -0.5, 0.5, TRUE)
    }
    # Residual standard deviation
    prior <- prior + dlnorm(params[length(params)], log(0.001), 2.5, TRUE)
    unname(prior)
  }

  prior <- BayesianTools::createPrior(
    density = prior_function,
    sampler = function() sample_params_edr(edr_run_pfts$pft_name, prospect_means, prospect_covar)
  )

  samples <- invert_bt(observed, invert_model, prior, custom_settings)
  samples
}
