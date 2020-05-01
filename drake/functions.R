last_result_file <- function() {
  info <- fs::dir_info(
    "multi_site_pda_results",
    recurse = TRUE,
    glob = "*.rds"
  )
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

soil_moisture_plot <- function(tidy_posteriors, site_structure) {
  site_posterior <- tidy_posteriors %>%
    dplyr::filter(grepl("sitesoil", variable)) %>%
    dplyr::inner_join(site_structure, c("variable" = "site_tag"))

  ggplot(site_posterior) +
    aes(x = forcats::fct_reorder(site_name, frac_evergreen_wtd), y = value,
        fill = frac_evergreen_wtd, color = frac_evergreen_wtd) +
    geom_violin() +
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

spec_error_all_f <- function(observed_predicted) {
  # Sort sites by aggregate bias
  plot_dat <- observed_predicted %>%
    group_by(site) %>%
    mutate(site_mean_bias = mean(bias)) %>%
    ungroup() %>%
    mutate(site_f = forcats::fct_reorder(site, site_mean_bias))
  ggplot(plot_dat) +
    aes(x = wavelength) +
    geom_ribbon(aes(ymin = pmax(albedo_r_q025, 0), ymax = albedo_r_q975),
                fill = "gray70") +
    ## geom_ribbon(aes(ymin = pmax(albedo_q025, 0), ymax = albedo_q975),
    ##             fill = "gray50") +
    geom_line(aes(y = observed, group = iobs)) +
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
    geom_ribbon(aes(ymin = pmax(albedo_r_q025, 0), ymax = albedo_r_q975),
                fill = "gray70") +
    geom_line(aes(y = observed, group = iobs)) +
    labs(x = "Wavelength (nm)", y = "Reflectance (0 - 1)",
         title = site) +
    coord_cartesian(ylim = c(0, 0.7)) +
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
