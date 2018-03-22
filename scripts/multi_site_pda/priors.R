load(here("priors/mvtraits_priors.RData"))
pfts <- rownames(means)
pfts <- pfts[!grepl("Southern_Pine", pfts)]

# Fix row and column names of means and covariances
mv_priors <- c(
  "prospect_N", "prospect_Cab", "prospect_Car",
  "prospect_Cw", "prospect_Cm", "SLA"
)
colnames(means) <- rownames(covars) <- colnames(covars) <- mv_priors

# Draw only positive values from multivariate prior
rmvnorm_positive <- function(mu, Sigma) {
  draw <- -1
  while(any(draw < 0)) {
    draw <- mvtnorm::rmvnorm(1, mu, Sigma)
  }
  draw
}

# Load allometry priors
allom_stats <- readRDS(here("priors/allometry_stats.rds")) %>%
  map(list(18, "statistics")) %>%
  .[pfts]

allom_names <- c("b1Bl", "b2Bl")
allom_mu <- map(allom_stats, ~.[c("mu0", "mu1"), "Mean"] %>% setNames(allom_names))
allom_Sigma <- map(allom_stats, ~matrix(.[c("tau11", "tau12", "tau12", "tau22"), "Mean"], 2, 2))

# Define priors for other parameters
prior_clumping <- c(0, 1)
rclumping <- function(n) runif(length(pfts), prior_clumping[1], prior_clumping[2])
dclumping <- function(x, log = TRUE) dunif(x, prior_clumping[1], prior_clumping[2], log = log)

prior_orient <- c(6, 4) * 3
rorient <- function(n) 2 * rbeta(n, prior_orient[1], prior_orient[2]) - 1
dorient <- function(x, log = TRUE) dbeta((x + 1) / 2, prior_orient[1], prior_orient[2], log = log)

prior_residual <- c(0.01, 0.01)

prior_sample <- function() {
  mv_draws <- -1
  mv_draws <- map(
    pfts,
    ~rmvnorm_positive(means[.,], covars[, , .]) %>% setNames(mv_priors)
  )
  allom_draws <- map2(
    allom_mu,
    allom_Sigma,
    ~mvtnorm::rmvnorm(1, .x, .y)[1, ]
  ) %>%
    map(~c(exp(.[1]), .[2]))
  all_draws <- pmap(
    list(
      mv_draws,
      allom_draws,
      clumping_factor = rclumping(length(pfts)),
      orient_factor = rorient(length(pfts))
    ),
    c
  )
  names(all_draws) <- pfts
  c(unlist(all_draws), residual = rgamma(1, prior_residual[1], prior_residual[2]))
}

prior_density <- function(params) {
  residual <- params["residual"]
  res_dens <- dgamma(residual, prior_residual[1], prior_residual[2], log = TRUE)
  if (!is.finite(res_dens)) {
    #message("Residual density is not finite.")
    #print(residual)
    return(-Inf)
  }
  traits <- params2edr(params, prospect = FALSE)$trait.values
  mvdens <- imap_dbl(
    traits,
    ~mvtnorm::dmvnorm(.x[mv_priors], means[.y, ], covars[, , .y], log = TRUE)
  )
  if (any(!is.finite(mvdens))) {
    #message("Multivariate density is not finite.")
    #print(params)
    return(-Inf)
  }
  allom_dens <- imap_dbl(
    traits,
    ~mvtnorm::dmvnorm(c(log(.x[allom_names[1]]), .x[allom_names[2]]),
                      allom_mu[[.y]], allom_Sigma[[.y]], log = TRUE)
  )
  if (any(!is.finite(allom_dens))) {
    return(-Inf)
  }
  cf_dens <- dclumping(map_dbl(traits, "clumping_factor"), log = TRUE)
  if (any(!is.finite(cf_dens))) {
    #message("Clumping factor density not finite.")
    #print(map_dbl(traits, "clumping_factor")[!is.finite(cf_dens)])
    return(-Inf)
  }
  of_dens <- dorient(map_dbl(traits, "orient_factor"), log = TRUE)
  if (any(!is.finite(of_dens))) {
    #message("Orient factor density not finite.")
    #print(map_dbl(traits, "orient_factor")[!is.finite(of_dens)])
    return(-Inf)
  }
  logdens <- sum(mvdens, allom_dens, cf_dens, of_dens, res_dens)
  if (!is.finite(logdens)) {
    #message("Log density is not finite.")
    return(-Inf)
  }
  logdens
}

prior <- BayesianTools::createPrior(
  density = prior_density,
  sampler = prior_sample
)
