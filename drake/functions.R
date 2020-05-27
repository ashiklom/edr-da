last_result_file <- function(basedir = "multi_site_pda_results") {
  info <- fs::dir_info(basedir, recurse = TRUE, glob = "*.rds")
  info %>%
    dplyr::arrange(change_time) %>%
    tail(1) %>%
    dplyr::pull(path)
}

tidy_prior <- function() {
  sites <- readLines(here::here("other_site_data", "site_list"))
  nsite <- length(sites)
  param_names <- readLines(here::here("param_names.txt"))
  prior <- create_prior(
    nsite = nsite,
    heteroskedastic = FALSE,
    limits = TRUE,
    param_names = param_names
  )
  prior_draws <- purrr::map(
    seq_len(2000),
    ~prior$sampler()
  ) %>% purrr::invoke(.f = rbind)
  tidy_param_matrix(prior_draws, "prior")
}

tidy_param_matrix <- function(mat, type) {
  tibble::as_tibble(mat) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    tidyr::pivot_longer(-id, names_to = "param", values_to = "value") %>%
    dplyr::mutate(type = !!type) %>%
    split_params("param")
}

pft_posterior_plot <- function(tidy_priors, tidy_posteriors) {
  tidy_prior_sub <- tidy_priors %>%
    dplyr::filter(
      !(variable == "b1Bl" & value > 0.75),
      !(variable == "b1Bw" & value > 1.5),
      !is.na(pft)
    )
  clrs <- c("prior" = "gray70", "posterior" = "black")
 
  ggplot() +
    aes(x = forcats::fct_inorder(pft), y = value,
        fill = type, color = type) +
    geom_violin(data = tidy_prior_sub) +
    geom_violin(data = dplyr::filter(tidy_posteriors, !is.na(pft))) +
    facet_wrap(
      vars(forcats::fct_inorder(variable)),
      scales = "free_y",
      ncol = 2
    ) +
    scale_color_manual(values = clrs, aesthetics = c("color", "fill")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}

soil_moisture_plot <- function(tidy_posteriors, site_structure_data) {
  site_structure <- site_structure_data %>%
    dplyr::mutate(x = forcats::fct_reorder(site_name, frac_evergreen_wtd))

  site_posterior <- tidy_posteriors %>%
    dplyr::filter(grepl("sitesoil", variable)) %>%
    dplyr::inner_join(site_structure, c("variable" = "site_tag"))

  last_hw_site <- site_structure %>%
    dplyr::filter(frac_evergreen_wtd <= 0.5) %>%
    dplyr::arrange(dplyr::desc(frac_evergreen_wtd)) %>%
    dplyr::slice(1) %>%
    dplyr::pull(x)

  site_posterior_summary <- site_posterior %>%
    dplyr::group_by(EG = frac_evergreen_wtd > 0.5) %>%
    dplyr::summarize(Mean = mean(value),
                     SD = sd(value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      x = as.numeric(last_hw_site) + 0.5 + c(-3, 3),
      lab = sprintf("%.2f (%.2f)", Mean, SD)
    )

  ggplot(site_posterior) +
    aes(x = x, y = value,
        fill = frac_evergreen_wtd, color = frac_evergreen_wtd) +
    geom_violin() +
    geom_vline(xintercept = as.numeric(last_hw_site) + 0.5,
               linetype = "dashed") +
    geom_text(aes(x = x, y = 0, label = lab), data = site_posterior_summary,
              inherit.aes = FALSE) +
    scale_color_viridis_c(
      aesthetics = c("color", "fill"),
      guide = guide_colorbar(title = "Weighted evergreen fraction",
                             direction = "horizontal",
                             title.position = "top")
    ) +
    labs(x = "Site code", y = "Soil moisture fraction (0 - 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      legend.position = c(0.98, 1),
      legend.justification = c(1, 1),
      legend.background = element_blank()
    )
}

full_site_info <- function(site_list_file, site_dir) {
  selected_sites <- readLines(site_list_file)
  site_files <- selected_sites %>% purrr::map_chr(
    ~head(list.files(file.path(site_dir, .x), "css$", full.names = TRUE), 1)
  ) %>%
    setNames(selected_sites)
  site_df <- purrr::map_dfr(site_files, read.table,
                            header = TRUE, .id = "site") %>%
    tibble::as_tibble() %>%
    dplyr::select(site, year = time, dbh, ipft = pft, nplant = n) %>%
    dplyr::mutate(
      ipft = get_ipft(ipft),
      hite = dbh2h(dbh, ipft)
    ) %>%
    dplyr::group_by(site, year) %>%
    dplyr::arrange(desc(hite)) %>%
    dplyr::mutate(cohort = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pft = factor(ipft, 1:5, c(
      "Early_Hardwood", "North_Mid_Hardwood", "Late_Hardwood",
      "Northern_Pine", "Late_Conifer"
    )))
  site_df
}

predict_lai <- function(site_details, tidy_posteriors, max_samples = 5000) {
  tidy_params_dt <- tidy_posteriors %>%
    dplyr::filter(variable %in% c("b1Bl", "SLA", "clumping_factor"))
  b2Bl <- purrr::map_dbl(allom_mu, "b2Bl")
  names(b2Bl) <- gsub("temperate\\.", "", names(b2Bl))
  params_structure <- tidy_params_dt %>%
    tidyr::pivot_wider(
      names_from = "variable",
      values_from = "value"
    ) %>%
    dplyr::mutate(b2Bl = b2Bl[pft])
  nsamp <- min(max_samples, nrow(params_structure))
  params_structure_sub <- params_structure %>%
    dplyr::sample_n(nsamp, replace = FALSE)
  site_lai_samples <- params_structure_sub %>%
    dplyr::left_join(site_details, "pft") %>%
    dplyr::mutate(
      bleaf = size2bl(dbh, b1Bl, b2Bl),
      lai = nplant * bleaf * SLA,
      elai = lai * clumping_factor
    )
  site_lai_samples
}

summarize_lai_samples <- function(site_lai_samples) {
  site_lai_samples %>%
    dplyr::group_by(site, year, pft, ipft, hite, dbh, nplant, cohort) %>%
    dplyr::summarize(
      lai_mean = mean(lai),
      lai_sd = sd(lai),
      lai_lo = quantile(lai, 0.025),
      lai_hi = quantile(lai, 0.975),
      elai_mean = mean(elai),
      elai_sd = sd(elai),
      elai_lo = quantile(elai, 0.025),
      elai_hi = quantile(elai, 0.975)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(site, year) %>%
    dplyr::arrange(hite) %>%
    dplyr::mutate(
      cum_lai = cumsum(lai_mean),
      cum_elai = cumsum(elai_mean)
    ) %>%
    dplyr::ungroup()
}

tidy_site_spec <- function(site) {
  rawobs <- load_observations(site)
  colnames(rawobs) <- as.character(seq_len(NCOL(rawobs)))
  dfobs <- as.data.frame(rawobs, row.names = PEcAnRTM::wavelengths(rawobs))
  tibble::as_tibble(dfobs, rownames = "wavelength") %>%
    dplyr::mutate(wavelength = as.numeric(wavelength)) %>%
    tidyr::pivot_longer(-wavelength, names_to = "iobs", values_to = "observed")
}

spec_error_all_f <- function(observed_predicted, sail_predictions) {
  # Sort sites by aggregate bias
  plot_dat <- observed_predicted %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(site_mean_bias = mean(bias)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(site_f = forcats::fct_reorder(site, site_mean_bias))
  sail_sub <- sail_predictions %>%
    dplyr::semi_join(observed_predicted, "wavelength") %>%
    dplyr::mutate(site_f = factor(site, levels(plot_dat[["site_f"]])))
  sail_avg <- sail_sub %>%
    tidyr::pivot_wider(names_from = "stream", values_from = "value") %>%
    # Same configuration as EDR -- assume incident radiation is 90% direct, 10%
    # diffuse
    dplyr::mutate(value = 0.9 * dhr + 0.1 * bhr)
  ggplot(plot_dat) +
    aes(x = wavelength) +
    geom_ribbon(aes(ymin = pmax(albedo_r_q025, 0),
                    ymax = pmin(albedo_r_q975, 1)),
                fill = "gray80") +
    geom_ribbon(aes(ymin = pmax(albedo_q025, 0),
                    ymax = pmin(albedo_q975, 1)),
                fill = "green3") +
    geom_line(aes(y = albedo_mean), color = "green4", size = 1) +
    geom_line(aes(y = observed, group = iobs)) +
    geom_line(aes(y = value), color = "red", data = sail_avg) +
    facet_wrap(vars(site_f), scales = "fixed", ncol = 6) +
    labs(x = "Wavelength (nm)", y = "Reflectance (0 - 1)") +
    theme_bw()
  # TODO Facet text inside plots --- use `geom_text`
}

spec_error_aggregate_f <- function(observed_predicted) {
  ggplot(observed_predicted) +
    aes(x = wavelength, y = bias, group = interaction(iobs, site)) +
    geom_line(alpha = 0.2) +
    geom_hline(yintercept = 0, color = "red") +
    labs(x = "Wavelength (nm)",
         y = expression("Predicted (mean)" - "observed reflectance")) +
    theme_bw()
}

site_spec_dbh_plot <- function(site, observed_predicted, site_details) {
  spec_sub <- observed_predicted %>%
    dplyr::filter(site == !!site)

  pspec <- ggplot(spec_sub) +
    aes(x = wavelength) +
    geom_ribbon(aes(ymin = pmax(albedo_r_q025, 0), ymax = pmin(albedo_r_q975, 1)),
                fill = "gray70") +
    geom_ribbon(aes(ymin = pmax(albedo_q025, 0),
                    ymax = pmin(albedo_q975, 1)),
                fill = "green3") +
    geom_line(aes(y = albedo_mean), color = "green4", size = 1) +
    geom_line(aes(y = observed, group = iobs)) +
    labs(x = "Wavelength (nm)", y = "Reflectance (0 - 1)",
         title = site) +
    coord_cartesian(ylim = c(0, 1.0)) +
    theme_bw()

  dbh_dat <- site_details %>%
    dplyr::filter(site == !!site)

  pft_colors <- RColorBrewer::brewer.pal(5, "Set1")
  names(pft_colors) <- levels(dbh_dat$pft)

  pdbh <- ggplot(dbh_dat) +
    aes(x = dbh, fill = pft) +
    geom_histogram(binwidth = 5) +
    coord_cartesian(xlim = c(0, 100)) +
    labs(x = "DBH (cm)", y = "Count", fill = "PFT") +
    scale_fill_manual(values = pft_colors, drop = FALSE) +
    theme_bw() +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.background = element_blank()
    )
  pspec + pdbh
}

calc_ndvi <- function(dat, vcol) {
  dat %>%
    dplyr::filter(wavelength %in% c(690, 800)) %>%
    tidyr::pivot_wider(
      names_from = "wavelength",
      values_from = all_of(vcol)
    ) %>%
    dplyr::rename(nir = `800`, red = `690`) %>%
    dplyr::mutate(ndvi = (nir - red) / (nir + red))
}

calc_ndvi_bysite <- function(observed_spectra, predicted_spectra,
                             site_structure) {
  obs_ndvi <- calc_ndvi(observed_spectra, "observed")
  pred_ndvi <- predicted_spectra %>%
    dplyr::select(wavelength, site, albedo_mean) %>%
    calc_ndvi("albedo_mean")
  obs_ndvi %>%
    dplyr::inner_join(pred_ndvi, "site", suffix = c("_obs", "_pred")) %>%
    dplyr::left_join(site_structure, c("site" = "site_name"))
}

ndvi_dbh_plot <- function(both_ndvi) {
  ggplot(both_ndvi) +
    aes(x = mean_dbh) +
    geom_point(aes(y = ndvi_obs, shape = "observed")) +
    geom_smooth(aes(y = ndvi_obs, linetype = "observed"),
                method = "lm", se = FALSE, color = "black") +
    geom_point(aes(y = ndvi_pred, shape = "predicted")) +
    geom_smooth(aes(y = ndvi_pred, linetype = "predicted"),
                method = "lm", se = FALSE, color = "black") +
    labs(x = "Mean DBH (cm)", y = "NDVI") +
    scale_shape_manual(values = c("observed" = 19, predicted = 3)) +
    theme_bw()
}

lai_predicted_observed_plot <- function(site_lai_total, lai_observed) {
  plot_dat <- dplyr::inner_join(site_lai_total, lai_observed, "site")
  ggplot(plot_dat) +
    aes(x = lai_mean, xmin = lai_lo, xmax = lai_hi,
        y = obs_LAI, ymin = obs_LAI_lo, ymax = obs_LAI_hi) +
    geom_pointrange() +
    geom_errorbarh() +
    geom_abline(linetype = "dashed") +
    labs(x = "Predicted LAI", y = "Observed LAI") +
    theme_bw()
}

tidy_sail_predictions <- function(site_details, site_lai_total,
                                  tidy_posteriors) {
  site_tags <- tibble::tibble(
    site = readLines(site_list_file),
    site_tag = paste0("sitesoil_", seq_along((site)))
  )

  sail_input <- site_details %>%
    group_by(site) %>%
    arrange(desc(hite), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(site_lai_total, c("site", "year")) %>%
    left_join(site_tags, "site")

  posterior_means <- tidy_posteriors %>%
    group_by(pft, variable) %>%
    summarize(value = mean(value)) %>%
    ungroup()

  soil_means <- posterior_means %>%
    filter(grepl("sitesoil", variable)) %>%
    select(site_tag = variable, soil_moisture = value)

  posterior_params <- posterior_means %>%
    filter(grepl("prospect_", variable) | variable == "orient_factor") %>%
    tidyr::pivot_wider(names_from = "variable", values_from = "value")

  sail_input2 <- sail_input %>%
    left_join(posterior_params, "pft") %>%
    left_join(soil_means, "site_tag") %>%
    mutate(leaf_theta = acos((1 + orient_factor) / 2))

  sail_output <- sail_input2 %>%
    mutate(
      sail_result = purrr::pmap(list(
        prospect_N, prospect_Cab, prospect_Car, prospect_Cw, prospect_Cm,
        leaf_theta, elai_mean, soil_moisture
      ), ~PEcAnRTM::pro4sail(c(..1, ..2, ..3, 0, ..4, ..5,  # PROSPECt
                               ..6, 0, 2,  # Leaf angle
                               ..7, 0,     # LAI, hot spot
                               0, 0, 0,    # Angles
                               ..8)))
    )

  sail_output_proc <- sail_output %>%
    select(site, sail_result) %>%
    mutate(
      wavelength = purrr::map(sail_result, PEcAnRTM::wavelengths),
      bhr = purrr::map(sail_result, ~.x[, "bi-hemispherical"]) %>% purrr::map(as.numeric),
      dhr = purrr::map(sail_result, ~.x[, "directional_hemispherical"]) %>% purrr::map(as.numeric),
      hdr = purrr::map(sail_result, ~.x[, "hemispherical_directional"]) %>% purrr::map(as.numeric),
      bdr = purrr::map(sail_result, ~.x[, "bi-directional"]) %>% purrr::map(as.numeric)
    ) %>%
    select(-sail_result) %>%
    tidyr::unnest(wavelength:bdr)

  sail_output_proc %>%
    tidyr::pivot_longer(bhr:bdr, names_to = "stream", values_to = "value")
}
