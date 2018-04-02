#' Create prior functions
#'
#' @param fix_allom2 Logical. If `TRUE`, fix allometry exponent to prior mean. 
#' If `FALSE`, use full multivariate prior.
#' @param heteroskedastic Logical. If `TRUE`, use heteroskedastic error model 
#' (with residual slope and intercept). If `FALSE`, use scalar residual.
#' @export
create_prior <- function(fix_allom2 = TRUE, heteroskedastic = TRUE) {
  BayesianTools::createPrior(
    density = create_prior_density(fix_allom2, heteroskedastic),
    sampler = create_prior_sampler(fix_allom2, heteroskedastic)
  )
}

#' @rdname create_prior
#' @export
create_prior_sampler <- function(fix_allom2 = TRUE, heteroskedastic = TRUE) {
  function(n = 1) {
    out <- numeric()
    for (i in seq_len(n)) {
      prospect_params <- rprospect()
      alloms <- if (fix_allom2) rallom1() else rallom2()
      if (fix_allom2) allom_names <- allom_names[1]
      alloms <- purrr::map(alloms, setNames, allom_names)
      cf <- rclumping() %>% purrr::map(setNames, "clumping_factor")
      of <- rorient() %>% purrr::map(setNames, "orient_factor")
      resid <- if (heteroskedastic) rresidual2() else rresidual()
      curr_params <- purrr::pmap(
        list(
          prospect_params,
          alloms,
          cf,
          of
        ),
        c
      ) %>% unlist()
      out <- c(out, curr_params, resid)
    }
    out
  }
}

#' @rdname create_prior
#' @export
create_prior_density <- function(fix_allom2 = TRUE, heteroskedastic = TRUE) {
  function(params) {
    ld_resid <- if (heteroskedastic) dresidual2(params) else dresidual(params)
    traits <- PEcAnRTM::params2edr(params, prospect = FALSE)$trait.values
    ld_allom <- if (fix_allom2) dallom1(traits) else dallom2(traits)
    ld_prosp <- dprospect(traits)
    ld_clumping <- dclumping(traits)
    ld_orient <- dorient(traits)
    sum(ld_resid, ld_allom, ld_prosp, ld_clumping, ld_orient)
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
  invalid <- list()
  n_invalid <- 0
  if (progress) pb <- txtProgressBar()
  for (i in seq_len(n_test)) {
    if (progress) setTxtProgressBar(pb, i / n_test)
    test_params <- prior$sampler()
    test_priors <- prior$density(test_params)
    prior_samps[i, ] <- test_params
    if (!is.numeric(test_priors) || !is.finite(test_priors)) {
      msg <- paste0("Problem with priors at index ", i)
      if (error) {
        print(test_params)
        stop(msg)
      } else {
        warning(msg)
        n_invalid <- n_invalid + 1
        invalid[[n_invalid]] <- test_params
      }
    }
    #stopifnot(is.numeric(test_priors), is.finite(test_priors))
  }
  invisible(invalid)
}
