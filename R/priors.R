#' Multivariate allometry prior (varying both base and exponent)
#'
#' @param traits Traits list
#' @export
rallom2 <- function() {
  purrr::map2(
    allom_mu,
    allom_Sigma,
    ~mvtnorm::rmvnorm(1, .x, .y)[1, ]
  ) %>%
    purrr::map(~c(exp(.[1]), .[2]))
}

#' @rdname rallom2
#' @export
dallom2 <- function(traits, log = TRUE) {
  b1Bl <- purrr::map_dbl(traits, allom_names[1]) %>% log()
  b2Bl <- purrr::map_dbl(traits, allom_names[2])
  pfts <- names(traits)
  purrr::pmap_dbl(
    list(b1Bl, b2Bl, pfts),
    ~mvtnorm::dmvnorm(c(..1, ..2), allom_mu[[..3]], allom_Sigma[[..3]], log = log)
  )
}

#' Allometry mean vector
#'
#' @name allom_mu
#' @export allom_mu
NULL

#' Allometry variance-covariance matrix
#'
#' @name allom_Sigma
#' @export allom_Sigma
NULL

#' Univariate allometry prior (varying just the base)
#'
#' @inheritParams rallom2
#' @export
rallom1 <- function() {
  out <- rnorm(npfts, b1Bl_means, b1Bl_sds) %>% exp()
  names(out) <- paste(pfts, allom_names[1], sep = ".")
  out
}

#' @rdname rallom1
#' @export
dallom1 <- function(traits, log = TRUE) {
  allom_vals <- purrr::map_dbl(traits, allom_names[1]) %>% log()
  ld <- dnorm(allom_vals, b1Bl_means, b1Bl_sds, log = log)
  ld
}

#' @export
prior_clumping <- c(0, 1)

#' Priors on other EDR parameters
#'
#' @inheritParams rallom2
#' @export
rclumping <- function() {
  out <- runif(npfts, prior_clumping[1], prior_clumping[2])
  names(out) <- paste(pfts, "clumping_factor", sep = ".")
  out
}

#' @rdname rclumping
#' @export
dclumping <- function(traits, log = TRUE) {
  x <- purrr::map_dbl(traits, "clumping_factor")
  out <- dunif(x, prior_clumping[1], prior_clumping[2], log = log)
  names(out) <- pfts
  out
}

#' @export
prior_orient <- c(6, 4) * 3

#' @rdname rclumping
#' @export
rorient <- function() {
  out <- 2 * rbeta(length(pfts), prior_orient[1], prior_orient[2]) - 1
  names(out) <- paste(pfts, allom_names[1], sep = ".")
  out
}

#' @rdname rclumping
#' @export
dorient <- function(traits, log = TRUE) {
  x <- purrr::map_dbl(traits, "orient_factor")
  out <- dbeta((x + 1) / 2, prior_orient[1], prior_orient[2], log = log)
  names(out) <- pfts
  out
}

#' @export
prior_residual <- c(0.01, 0.01)

#' @rdname rclumping
#' @export
rresidual <- function() {
  rgamma(1, prior_residual[1], prior_residual[2]) %>%
    setNames("residual")
}

#' @rdname rclumping
#' @export
dresidual <- function(params, log = TRUE) {
  x <- params["residual"]
  dgamma(x, prior_residual[1], prior_residual[2], log = log) %>%
    setNames("residual")
}

#' @export
prior_residual2 <- list(
  intercept = c(0, 0.5),
  slope = c(0, 10)
)

#' @rdname rclumping
#' @export
rresidual2 <- function() {
  rint <- rnorm(1, prior_residual2$intercept[1], prior_residual2$intercept[2])
  rslope <- rnorm(1, prior_residual2$slope[1], prior_residual2$slope[2])
  c("residual_intercept" = rint, "residual_slope" = rslope)
}

#' @rdname rclumping
#' @export
dresidual2 <- function(params, log = TRUE) {
  rint <- dnorm(
    params["residual_intercept"],
    prior_residual2$intercept[1],
    prior_residual2$intercept[2],
    log = log
  )
  rslope <- dnorm(
    params["residual_slope"],
    prior_residual2$slope[1],
    prior_residual2$slope[2],
    log = log
  )
  c("residual_intercept" = rint, "residual_slope" = rslope)
}

#' @rdname rclumping
#' @export
rprospect <- function() {
  out <- purrr::map(
    pfts,
    ~rmvnorm_positive(prospect_means[., ], prospect_covars[, , .])
  )
  names(out) <- pfts
  out
}

#' @rdname rclumping
#' @export
dprospect <- function(traits, log = TRUE) {
  out <- purrr::imap_dbl(
    traits,
    ~mvtnorm::dmvnorm(.x[prospect_names], prospect_means[.y, ], prospect_covars[, , .y], log = TRUE)
  )
  names(out) <- pfts
  out
}

#' Draw only positive values from multivariate prior
rmvnorm_positive <- function(mu, Sigma) {
  draw <- -1
  while(any(draw < 0)) {
    draw <- mvtnorm::rmvnorm(1, mu, Sigma)[1,]
  }
  draw
}
