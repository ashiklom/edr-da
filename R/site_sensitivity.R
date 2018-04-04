#' Perform a single EDR simulation at a single site
#'
#' @param params Parameter vector
#' @param site String identifying the site code
#' @param pda_dir Path to PDA run output
#' @param img_path Path to EDR singularity image
#' @param edr_exe_path Path to EDR executable
#' @export
run_site <- function(params,
                     site,
                     pda_dir,
                     img_path = getOption("redr.img_path"),
                     edr_exe_path = getOption("redr.edr_exe_path")
                     ) {
  stopifnot(!(is.null(img_path) && is.null(edr_exe_path)))
  ed2in_path <- file.path(pda_dir, site, "edr", "ED2IN")
  stopifnot(file.exists(ed2in_path))
  edr_in <- PEcAnRTM::params2edr(params, prospect = TRUE, version = 5)
  EDR(
    img_path = img_path,
    ed2in_path = ed2in_path,
    spectra_list = edr_in$spectra_list,
    trait.values = edr_in$trait.values,
    verbose_error = FALSE,
    edr_exe_path = edr_exe_path
  )
}

#' Perform site sensitivity analysis over a set of samples for a given site
#'
#' @param samples_mat Matrix of samples. Each iteration is a row, each column is a parameter.
#' @param site Site to simulate
#' @param n Number of simulations to do. Default = 500.
#' @export
site_sensitivity <- function(samples_mat, site, pda_dir, n = 500, ...) {
  out <- matrix(0, 2101, n)
  pb <- progress::progress_bar$new(total = n)
  pb$tick(0)
  for (i in seq_len(n)) {
    pb$tick()
    spec <- NULL
    isim <- 0
    while (is.null(spec) && isim <= 20) {
      isim <- isim + 1
      r <- sample(nrow(samples_mat), 1)
      p <- samples_mat[r, ]
      spec <- tryCatch(
        run_site(p, site, pda_dir, ...),
        error = function(e) {message("Error in run. Trying again."); NULL}
      )
    }
    if (isim > 10) stop("More than 20 bad simulations in a row. Stopping because something is wrong.")
    out[, i] <- spec
  }
  PEcAnRTM::spectra(out, 400:2500, "R")
}

#' Plot sensitivity confidence bands
#'
#' @param result Matrix of resulting spectra
#' @param obs Matrix of observations
#' @param ... Additional arguments to `plot`
#' @export
plot_sens <- function(result, obs, ...) {
  rmu <- rowMeans(result, na.rm = TRUE)
  rlo <- apply(result, 1, quantile, 0.025, na.rm = TRUE)
  rhi <- apply(result, 1, quantile, 0.975, na.rm = TRUE)
  y_lims <- range(rmu, rlo, rhi, obs)
  plot(0, 0, type = "n", xlim = c(400, 2500), ylim = y_lims, ...)
  polygon(c(400:2500, 2500:400), c(rlo, rev(rhi)), col = "green")
  lines(400:2500, rmu, col = "green4")
  matplot(obs_list[[1]], col = "black", lty = "solid", add = TRUE)
}
