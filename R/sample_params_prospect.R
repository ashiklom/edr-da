#' Draw parameters from prior
#'
#' @param pft Name of PFT. Must be in `prospect_means`
#' @param prospect_means Matrix of PROSPECT parameter prior means
#' @param prospect_covars Array of PROSPECT parameter prior covariance matrices
#' @export
sample_params_prior <- function(pft, prospect_means, prospect_covars) {
  params <- -1
  while (any(params < 0)) {
    params <- mvtnorm::rmvnorm(1, prospect_means[pft, ], prospect_covars[,,pft])[1, ]
  }
  params
}
