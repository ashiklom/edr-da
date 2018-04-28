#' Retrieve posterior or prior confidence interval for a site
#'
#' @param site Site code
#' @param pda_dir Directory of PDA runs
#' @param burnin Number of samples to burnin
#' @param alpha P value for quantiles
#' @export
get_post_ci <- function(site, pda_dir, burnin, alpha = 0.05) {
  sample_spec <- read_spectra_history(site, pda_dir, burnin)
  summarize_specmat(sample_spec, alpha = alpha, waves = 400:2500)
}

#' Read spectra history
#'
#' @param site Site code
#' @param pda_dir Directory of PDA runs
#' @param burnin Number of samples to burnin
#' @return Matrix of spectra for each iteration
#' @export
read_spectra_history <- function(site, pda_dir, burnin) {
  specfile_path <- file.path(pda_dir, site, "edr", "spec_store")
  stopifnot(file.exists(specfile_path))
  nlines <- wcl(specfile_path)
  skip <- ifelse(nlines < burnin, 0, burnin)
  sample_spec_raw <- data.table::fread(specfile_path, header = FALSE,
                                       blank.lines.skip = TRUE, skip = skip,
                                       showProgress = FALSE)
  t(as.matrix(sample_spec_raw))
}

#' Count the number of lines in a file
#'
#' @param file Path to file
#' @return Number of lines in file
#' @export
wcl <- function(file) {
  string <- system2("wc", c("-l", file), stdout = TRUE)
  string2 <- gsub(" .*", "", string)
  as.numeric(string2)
}

#' @rdname get_post_ci
#' @export
get_prior_ci <- function(site, pda_dir, burnin, alpha = 0.05) {
  fname <- paste("prior_sim", site, "rds", sep = ".")
  specfile_path <- file.path(pda_dir, fname)
  sample_spec <- readRDS(specfile_path)
  waves <- PEcAnRTM::wavelengths(sample_spec)
  summarize_specmat(sample_spec, alpha = alpha, waves = waves)
}

#' Summarize a matrix of spectra
#'
#' @param mat Matrix of spectra, waves x samples
#' @param waves Wavelengths. Default = 400:2500
#' @inheritParams get_post_ci
#' @param ... Additional arguments to [PEcAnRTM::spectra]
#' @export
summarize_specmat <- function(mat, alpha = 0.05, waves = 400:2500, ...) {
  mu <- rowMeans(mat, na.rm = TRUE)
  lo <- apply(mat, 1, quantile, alpha / 2, na.rm = TRUE)
  hi <- apply(mat, 1, quantile, 1 - (alpha / 2), na.rm = TRUE)
  out <- cbind(mu = mu, lo = lo, hi = hi)
  PEcAnRTM::spectra(out, waves, ...)
}

#' Add spectra confidence interval ribbon to plot
#'
#' @param spectra_ci Spectra confidence limits object
#' @param ... Additional arguments to `polygon`
#' @export
add_spectra_ci <- function(spectra_ci, ...) {
  stopifnot(
    PEcAnRTM::is_spectra(spectra_ci),
    all(c("mu", "lo", "hi") %in% colnames(spectra_ci))
  )
  x <- PEcAnRTM::wavelengths(spectra_ci)
  xp <- c(x, rev(x))
  yp <- c(spectra_ci[, "lo"], rev(spectra_ci[, "hi"]))
  polygon(xp, yp, ...)
}

#' Plot prior, posterior, and observations for a given site
#'
#' @inheritParams get_post_ci
#' @param prior_color Color for prior band
#' @param post_color Color for posterior band
#' @param obs_color Color for observation line
#' @param obs_lwd Linewidth of observation line
#' @param pb Progress bar object (default = `NULL`)
#' @export
plot_site_spectra <- function(site, pda_dir, burnin, 
                              prior_color = alpha("blue", 0.2),
                              post_color = alpha("red", 0.5),
                              obs_color = "black",
                              obs_lwd = 2,
                              pb = NULL,
                              ...) {
  if (!is.null(pb)) pb$tick()
  prior_ci <- get_prior_ci(site, pda_dir, argl$burnin, ...)
  post_ci <- get_post_ci(site, pda_dir, argl$burnin, ...)
  obs <- load_observations(site)
  has_obs <- ncol(obs) > 0
  if (!has_obs) {
    message("No observations for site: ", site)
  }
  xrange <- c(400, 2500)
  yrange <- range(prior_ci, post_ci, obs, na.rm = TRUE)
  plot(0, 0, type = "n", xlab = "Wavelength", ylab = "Reflectance", main = site,
       xlim = xrange, ylim = yrange)
  add_spectra_ci(prior_ci, col = prior_color)
  add_spectra_ci(post_ci, col = post_color)
  if (has_obs) {
    matlines(PEcAnRTM::wavelengths(obs), obs, col = obs_color, lwd = obs_lwd, lty = "solid")
  }
  legend(
    "topright",
    legend = c("prior", "posterior", "observation"),
    lwd = c(7, 7, 2),
    col = c(prior_color, post_color, obs_color)
  )
}
