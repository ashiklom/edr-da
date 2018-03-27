run_site <- function(params, site) {
  edr_in <- params2edr(params, prospect = TRUE, version = 5)
  EDR(
    img_path = img_path,
    ed2in_path = site,
    spectra_list = edr_in$spectra_list,
    trait.values = edr_in$trait.values,
    verbose_error = FALSE,
    edr_exe_path = edr_exe_path
  )
}

site_sensitivity <- function(samples, sites, burnin = 10000, n = 500) {
  samples_mat_full <- BayesianTools::getSample(samples, start = burnin)
  out_list <- list()
  for (si in seq_along(sites)) {
    s <- sites[si]
    message("Simulating site ", si, " of ", length(sites), ": ", names(s))
    out <- matrix(0, 2101, n)
    pb <- dplyr::progress_estimated(n)
    for (i in seq_len(n)) {
      pb$tick()$print()
      spec <- NULL
      while (is.null(spec)) {
        r <- sample(nrow(samples_mat_full), 1)
        p <- samples_mat_full[r, ]
        spec <- tryCatch(
          run_site(p, s),
          error = function(e) NULL
        )
      }
      out[, i] <- spec
    }
    out_list[[s]] <- out
  }
  out_list
}

get_obs <- function(site, use_wl = seq(400, 1300, by = 10)) {
  aviris_specfile <- here("aviris/aviris.h5")
  aviris_h5 <- hdf5r::H5File$new(aviris_specfile)
  site_tags <- gsub("([[:alnum:]]+)_.*", "\\1", site)
  av_inds <- purrr::map(site_tags, ~which(aviris_h5[["iPLOT"]][] == .))
  aviris_spec_raw <- purrr::map(av_inds, ~aviris_h5[["spectral_data"]][, .])
  aviris_wl <- aviris_h5[["wavelengths"]][]
  aviris_h5$close_all()
  aviris_spec <- purrr::map(aviris_spec_raw, spectra, wavelengths = aviris_wl)
  aviris_spec <- purrr::map(aviris_spec, resample, to = use_wl)
  observed <- do.call(cbind, aviris_spec) / 10000
  colnames(observed) <- rep(site, ncol(observed))
  observed
}

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
