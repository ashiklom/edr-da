#' Create prior functions
#'
#' @param fix_allom2 Logical. If `TRUE`, fix allometry exponent to prior mean. 
#' If `FALSE`, use full multivariate prior.
#' @param heteroskedastic Logical. If `TRUE`, use heteroskedastic error model 
#' (with residual slope and intercept). If `FALSE`, use scalar residual.
#' @param verbose Logical. If `TRUE`, print verbose error messages for invalid 
#' prior densities
#' @export
create_prior <- function(fix_allom2 = TRUE,
                         heteroskedastic = TRUE,
                         verbose = TRUE,
                         nsite = 1,
                         limits = FALSE,
                         best = TRUE,
                         param_names = NULL,
                         site_specific_var = FALSE) {
  lower <- NULL
  upper <- NULL
  if (limits) {
    npft <- 5
    pft_lower <- c(
      1, 0, 0, 0, 0, # prospect
      0, # SLA
      0, if (fix_allom2) NULL else -50, # leaf allom
      0, if (fix_allom2) NULL else -50, # wood allom
      0, -1 # clumping, orient
    )
    pft_upper <- c(
      10, 1000, 500, 0.1, 0.1, # prospect
      10000, # SLA
      50, if (!fix_allom2) 50, # leaf allom
      50, if (!fix_allom2) 50, # wood allom
      1, 0.6 # clumping, orient
    )
    resid_lower <- if (heteroskedastic) c(0, 0) else 0
    resid_upper <- if (heteroskedastic) c(100, 100) else 100
    if (site_specific_var) {
      resid_lower <- rep(resid_lower, each = nsite)
      resid_upper <- rep(resid_upper, each = nsite)
    }
    lower <- c(rep(pft_lower, npft), rep(0, nsite), resid_lower)
    upper <- c(rep(pft_upper, npft), rep(1, nsite), resid_upper)
  }
  resid_best <- if (heteroskedastic) c(0.01, 0.01)
  if (site_specific_var) resid_best <- rep(resid_best, each = nsite)
  bestval <- NULL
  if (best) {
    pft_best <- c(
      1.4, 40, 10, 0.01, 0.01,
      40,
      1, if (!fix_allom2) 1,
      1, if (!fix_allom2) 1,
      0.5, 0
    )
    bestval <- c(
      rep(pft_best, npft),
      rep(0.5, nsite),
      resid_best
    )
  }
  BayesianTools::createPrior(
    density = create_prior_density(
      fix_allom2 = fix_allom2,
      heteroskedastic = heteroskedastic,
      verbose = verbose,
      param_names = param_names,
      site_specific_var = site_specific_var
    ),
    sampler = create_prior_sampler(fix_allom2, heteroskedastic, nsite, site_specific_var = site_specific_var),
    lower = lower, upper = upper, best = bestval
  )
}

#' @rdname create_prior
#' @export
create_prior_sampler <- function(fix_allom2 = TRUE, heteroskedastic = TRUE, nsite = 1,
                                 site_specific_var = FALSE) {
  function(n = 1) {
    out <- numeric()
    for (i in seq_len(n)) {
      prospect_params <- rprospect()
      alloms <- if (fix_allom2) rallom1() else rallom2()
      walloms <- if (fix_allom2) rwallom1() else rwallom2()
      if (fix_allom2) allom_names <- allom_names[1]
      if (fix_allom2) wallom_names <- wallom_names[1]
      alloms <- purrr::map(alloms, setNames, allom_names)
      walloms <- purrr::map(walloms, setNames, wallom_names)
      cf <- rclumping() %>% purrr::map(setNames, "clumping_factor")
      of <- rorient() %>% purrr::map(setNames, "orient_factor")
      # HACK: Hard-coded uniform prior on soil moisture
      soil <- runif(nsite, 0, 1)
      names(soil) <- paste0("sitesoil_", seq_len(nsite))
      nresid <- if (site_specific_var) nsite else 1
      resid <- if (heteroskedastic) rresidual2(nresid) else rresidual(nresid)
      if (!site_specific_var) {
        names(resid) <- gsub("[[:digit:]]+", "", names(resid))
      }
      curr_params <- purrr::pmap(
        list(
          prospect_params,
          alloms,
          walloms,
          cf,
          of
        ),
        c
      ) %>% unlist()
      out <- c(out, curr_params, soil, resid)
    }
    out
  }
}

#' @rdname create_prior
#' @export
create_prior_density <- function(fix_allom2 = TRUE, heteroskedastic = TRUE, verbose = TRUE,
                                 param_names = NULL, site_specific_var = FALSE) {
  function(params) {
    if (!is.null(param_names)) {
      stopifnot(length(param_names) == length(params))
      names(params) <-  param_names
    }
    # HACK: Assume a uniform 0-1 prior for site soil, so no effect on likelihood
    soil_params <- params[grepl("sitesoil", names(params))]
    if (any(!is.finite(soil_params) | soil_params < 0 | soil_params > 1)) {
      if (verbose) message("Invalid soil prior")
      return(-Inf)
    }
    # Residual
    if (heteroskedastic) {
      residual_slope <- params[grep("residual_slope", names(params))]
      residual_intercept <- params[grep("residual_intercept", names(params))]
      ld_resid <-
        sum(dexp(
          residual_slope,
          prior_residual2$slope,
          log = TRUE
        )) +
        sum(dexp(
          residual_intercept,
          prior_residual2$intercept,
          log = TRUE
        ))
    } else {
      residuals <- params[grep("residual", names(params))]
      ld_resid <- sum(dgamma(
        residuals,
        prior_residual[1],
        prior_residual[2],
        log = TRUE
      ))
    }
    if (!is.finite(ld_resid)) {
      if (verbose) message("Invalid residual prior")
      print(params)
      return(-Inf)
    }
    # Drop the residuals and soil parameters from the remaining parameter vector
    params <- params[!grepl("residual|sitesoil", names(params))]
    traits <- PEcAnRTM::params2edr(params, prospect = FALSE)$trait.values
    ld_allom <- if (fix_allom2) dallom1(traits) else dallom2(traits)
    if (!all(is.finite(ld_allom))) {
      if (verbose) {
        message("Invalid allometry prior")
        print(purrr::map_dbl(traits, allom_names[1]))
        print(purrr::map_dbl(traits, allom_names[2]))
      }
      return(-Inf)
    }
    ld_wallom <- if (fix_allom2) dwallom1(traits) else dwallom2(traits)
    if (!all(is.finite(ld_wallom))) {
      if (verbose) {
        message("Invalid wood allometry prior")
        print(purrr::map_dbl(traits, wallom_names[1]))
        print(purrr::map_dbl(traits, wallom_names[2]))
      }
      return(-Inf)
    }
    ld_prosp <- dprospect(traits)
    if (!all(is.finite(ld_prosp))) {
      if (verbose) message("Invalid prospect prior")
      return(-Inf)
    }
    ld_clumping <- dclumping(traits)
    if (!all(is.finite(ld_clumping))) {
      if (verbose) message("Invalid clumping prior")
      return(-Inf)
    }
    ld_orient <- dorient(traits)
    if (!all(is.finite(ld_orient))) {
      if (verbose) message("Invalid orient prior")
      return(-Inf)
    }
    sum(ld_resid, ld_allom, ld_wallom, ld_prosp, ld_clumping, ld_orient)
  }
}


#' Check prior for validity
#'
#' @param prior BayesianTools prior object
#' @param n_test Number of samples to test
#' @param error Logical. If `TRUE`, throw an error if priors are invalid. If 
#' `FALSE`, throw a warning and store the invalid parameters.
#' @param progress Logical. If `TRUE`, print a progress bar.
#' @export
check_prior <- function(prior, n_test = 500, error = FALSE, progress = TRUE) {
  stopifnot(inherits(prior, "prior"))
  samp <- prior$sampler()
  n_prior <- length(samp)
  prior_samps <- matrix(0, n_test, n_prior)
  colnames(prior_samps) <- names(samp)
  n_invalid <- 0
  if (progress) {
    pb <- txtProgressBar()
    on.exit(close(pb))
  }
  for (i in seq_len(n_test)) {
    if (progress) setTxtProgressBar(pb, i / n_test)
    test_params <- prior$sampler()
    prior_samps[i, ] <- test_params
    test_priors <- tryCatch(
      prior$density(test_params),
      error = function(e) {
        if (error) {
          print(test_params)
          stop(e)
        } else {
          message("Hit the following error at index: ", i, " : ", e)
          NULL
        }
      }
    )
    if (is.null(test_priors)) n_invalid <- n_invalid + 1
  }
  #stopifnot(is.numeric(test_priors), is.finite(test_priors))
  attr(prior_samps, "n_invalid") <- n_invalid
  invisible(prior_samps)
}
