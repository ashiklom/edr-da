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

pft_posterior_plot <- function(tidy_priors, tidy_posteriors, ncol = 2) {
  lvls <- c("prospect_N", "prospect_Cab", "prospect_Car",
            "prospect_Cw", "prospect_Cm",
            "SLA",
            "b1Bl", "b1Bw",
            "clumping_factor", "orient_factor")
  lbls <- c("'# mesophyll layers'",
            "Chlorophyll ~ (mu * g ~ cm^-2)",
            "Carotenoids ~ (mu * g ~ cm^-2)",
            "'Water' ~ (g ~ cm^-2)",
            "'Dry matter' ~ (g ~ cm^-2)",
            "'Specific leaf area' ~ (kg ~ m^-2)",
            "'Leaf biomass allometry'", "'Wood biomass allometry'",
            "'Canopy clumping' ~ ('0, 1')", "'Leaf orientation' ~ ('-1, 1')")
  tidy_prior_sub <- tidy_priors %>%
    dplyr::filter(
      # Clipped because priors are much wider than posteriors
      !(variable == "b1Bl" & value > 0.3),
      !(variable == "b1Bw" & value > 0.4),
      !is.na(pft)
    ) %>%
    dplyr::mutate(variable = factor(variable, lvls, lbls))
  clrs <- c("prior" = "gray70", "posterior" = "black")

  tidy_posterior2 <- tidy_posteriors %>%
    dplyr::filter(!is.na(pft)) %>%
    dplyr::mutate(variable = factor(variable, lvls, lbls))

  ggplot() +
    aes(x = forcats::fct_inorder(pft), y = value,
        fill = type, color = type) +
    geom_violin(data = tidy_prior_sub) +
    geom_violin(data = tidy_posterior2) +
    facet_wrap(vars(variable), scales = "free_y", ncol = ncol,
               labeller = label_parsed) +
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
      guide = guide_colorbar(title = "Weighted evergreen fraction")
    ) +
    labs(x = "Site code", y = "Soil moisture fraction (0 - 1)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
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

spec_error_all_f <- function(observed_predicted, sail_predictions, ncol = 6) {
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
                    ymax = pmin(albedo_r_q975, 1),
                    fill = "95% PI")) +
    geom_ribbon(aes(ymin = pmax(albedo_q025, 0),
                    ymax = pmin(albedo_q975, 1),
                    fill = "95% CI")) +
    geom_line(aes(y = albedo_mean, color = "EDR"), size = 1) +
    geom_line(aes(y = observed, group = iobs, color = "AVIRIS")) +
    geom_line(aes(y = value, color = "SAIL"), data = sail_avg) +
    facet_wrap(vars(site_f), scales = "fixed", ncol = ncol) +
    labs(x = "Wavelength (nm)", y = "Reflectance (0 - 1)") +
    scale_fill_manual(
      name = "",
      values = c("95% PI" = "gray80",
                 "95% CI" = "green3")
    ) +
    scale_color_manual(
      name = "",
      values = c("EDR" = "green4",
                 "AVIRIS" = "black",
                 "SAIL" = "red")
    ) +
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

site_spec_dbh_plot <- function(site, observed_predicted, site_details,
                               spec_additions = NULL,
                               dbh_additions = NULL) {
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
  if (!is.null(spec_additions)) {
    pspec <- Reduce("+", c(list(pspec), spec_additions))
  }

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

  if (!is.null(dbh_additions)) {
    pdbh <- Reduce("+", c(list(pdbh), dbh_additions))
  }
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
  fit <- lm(obs_LAI ~ lai_mean, data = plot_dat)
  sfit <- summary(fit)
  eqn <- paste(
    sprintf("y = %.2fx + %.2f", fit$coefficients[2], fit$coefficients[1]),
    sprintf("R2 = %.2f, p = %.3f", sfit$r.squared,
            sfit$coefficients["lai_mean", "Pr(>|t|)"]),
    sep = "\n"
  )
  xx <- max(plot_dat$lai_hi)
  yy <- max(plot_dat$obs_LAI_hi)
  ggplot(plot_dat) +
    aes(x = lai_mean, xmin = lai_lo, xmax = lai_hi,
        y = obs_LAI, ymin = obs_LAI_lo, ymax = obs_LAI_hi) +
    geom_pointrange() +
    geom_errorbarh() +
    geom_abline(aes(linetype = "1:1", color = "1:1",
                    intercept = 0, slope = 1)) +
    geom_abline(aes(linetype = "Regression", color = "Regression",
                    intercept = fit$coefficients[1], slope = fit$coefficients[2])) +
    scale_linetype_manual(values = c("1:1" = "dashed", "Regression" = "solid"),
                          name = "") +
    scale_color_manual(values = c("1:1" = "black", "Regression" = "red"),
                       name = "") +
    annotate("text", x = xx, y = yy, hjust = 1, vjust = 1,
             label = eqn) +
    labs(x = "Predicted LAI", y = "Observed LAI") +
    theme_bw() +
    theme(legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.background = element_blank())
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

edr_sensitivity_defaults <- list(
  N = 1.4, Cab = 40, Car = 10, Cw = 0.01, Cm = 0.01,
  lai = 3, cai = 1,
  clumping_factor = 1,
  orient_factor = 0,
  direct_sky_frac = 0.8,
  pft = 1,
  czen = 0.85,
  wai = 0,
  soil_moisture = 0.5
)

sail_sensitivity_defaults <- c(
  edr_sensitivity_defaults[c("N", "Cab", "Car", "Cw", "Cm",
                             "soil_moisture")],
  list(
    Cbrown = 0,
    LAI = edr_sensitivity_defaults[["lai"]] * edr_sensitivity_defaults[["clumping_factor"]], #nolint
    hot_spot = 0,
    solar_zenith = acos(edr_sensitivity_defaults[["czen"]]) * 180/pi,
    LIDFa = edr_sensitivity_defaults[["orient_factor"]],
    LIDFb = 0
  )
)


do_sens <- function(value, variable, fun, .dots = list()) {
  stopifnot(is.list(.dots))
  varlist <- list()
  varlist[[variable]] <- value
  # Recycle PFT-specific variables if necessary
  nval <- length(value)
  if (nval > 1) {
    for (v in c("pft", "lai", "wai", "cai")) {
      .dots[[v]] <- rep(.dots[[v]], nval)
    }
  }
  arglist <- modifyList(.dots, varlist)
  do.call(fun, arglist)
}

tidy_albedo <- function(result_list, values) {
  stopifnot(length(result_list) == length(values))
  values_df <- tibble::tibble(
    variable = paste0("V", seq_along(values)),
    var_value = values
  )
  names(result_list) <- values_df[["variable"]]
  albedo_dfw <- purrr::map_dfc(result_list, "albedo")
  albedo_dfw[["wavelength"]] <- seq(400, 2500)
  albedo_long <- tidyr::pivot_longer(albedo_dfw, -wavelength,
                                     names_to = "variable", values_to = "value")
  dplyr::left_join(albedo_long, values_df, by = "variable")
}

tidy_sail <- function(result_list, values) {
  stopifnot(length(result_list) == length(values))
  results_dfl <- purrr::map(result_list, tibble::as_tibble) %>%
    purrr::map(function(x) {x$wavelength <- seq(400, 2500); x})
  values_df <- tibble::tibble(
    variable = paste0("V", seq_along(values)),
    var_value = values,
    saildata = results_dfl
  )
  tidyr::unnest(values_df, saildata)
}

sensitivity_plot <- function(values, varname, label,
                             defaults = edr_sensitivity_defaults,
                             ...) {
    sens <- purrr::map(
      values, do_sens,
      fun = edr_r,
      variable = varname,
      .dots = modifyList(defaults, list(...))
    ) %>%
      tidy_albedo(values)
    plt <- ggplot(sens) +
      aes(x = wavelength, y = value, color = var_value,
          group = var_value) +
      geom_line() +
      scale_color_viridis_c() +
      labs(x = "Wavelength (nm)", y = "Albedo [0,1]",
           color = label) +
      theme_bw() +
      theme(
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank()
      )
    ggsave(path(figdir, paste0("edr-sensitivity-", varname, ".png")), plt,
           width = 4, height = 4, dpi = 300)
}
