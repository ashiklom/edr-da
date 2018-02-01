#' Draw parameters for EDR simulation
#'
#' @param pft_vec Vector of PFTs for which to draw parameters
#' @inheritParams sample_params_prior
#' @export
sample_params_edr <- function(pft_vec, prospect_means, prospect_covars) {
  samples <- numeric()
  for (pft in pft_vec) {
    # PROSPECT prior
    prior_params <- sample_params_prior(pft, prospect_means, prospect_covars)
    draw <- c(
      prior_params,
      clumping_factor = runif(1, 0, 1),
      orient_factor = runif(1, -0.5, 0.5)
    )
    names(draw) <- paste(pft, names(draw), sep = ".")
    samples <- c(samples, draw)
  }
  samples <- c(samples, residual = rlnorm(1, log(0.001), 2.5))
  samples
}
