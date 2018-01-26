#' Draw a random sample of PRO4SAIL parameters
#'
#' @param lai Value of leaf area index
#' @param custom_sail Named list of sail parameters to replace defaults
#' @inheritParams sample_params_prior
#' @export
sample_params_sail <- function(pft, lai, prospect_means, prospect_covars,
                               custom_sail = list()) {
  prosp_params <- sample_params_prior(pft, prospect_means, prospect_covars)
  sail_params <- c(
    N = prosp_params["N.mu"],
    Cab = prosp_params["Cab.mu"],
    Car = prosp_params["Car.mu"],
    Cbrown = 0,
    Cw = prosp_params["Cw.mu"],
    Cm = prosp_params["Cm.mu"],
    LIDFa = runif(1, -0.9, 0.9),   # Default: -0.35
    LIDFb = runif(1, -0.9, 0.9),  # Default: -0.15
    TypeLIDF = 1,
    LAI = lai,
    q = runif(1, 0.01, 0.1),     # Hot spot, default = 0.01
    tts = 30,     # Solar zenith
    tto = 0,      # Observer zenith
    psi = 0,      # Sun-sensor azimuth
    psoil = 0.5   # Soil moisture
  )
  if (length(custom_sail) > 0) {
    stopifnot(
      !is.null(names(custom_sail)),
      all(names(custom_sail) %in% names(sail_params))
    )
    for (i in seq_along(custom_sail)) {
      if (is.function(custom_sail[i])) {
        value <- custom_sail[[i]]()
      } else {
        value <- custom_sail[[i]]
      }
      sail_params[names(custom_sail)[i]] <- value
    }
  }
  sail_params
}

