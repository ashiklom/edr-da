plan <- drake_plan(
  pda_result_file = last_result_file(),
  posterior_matrix = {
    samples_bt <- readRDS(pda_result_file)
    posterior_matrix <- BayesianTools::getSample(
      samples_bt,
      thin = "auto",
      start = 15000
    )
    colnames(posterior_matrix) <- readLines(file_in("param_names.txt"))
    posterior_matrix
  },
  tidy_posteriors = target(
    tidy_param_matrix(posterior_matrix, "posterior"),
    format = "fst_tbl"
  ),
  tidy_priors = target(tidy_prior(), format = "fst_tbl"),
  prior_posterior = ggsave(
    file_out(!!path(figdir, "posterior-pft.png")),
    pft_posterior_plot(tidy_priors, tidy_posteriors),
    width = 6, height = 7, dpi = 300
  ),
  site_structure_data = file_in(!!site_structure_file) %>%
    read_csv(col_types = "cnnnnnnclc") %>%
    mutate(site_tag = paste0("sitesoil_", row_number())),
  soil_moisture_plt = ggsave(
    file_out(!!path(figdir, "posterior-soil.png")),
    soil_moisture_plot(tidy_posteriors, site_structure_data),
    width = 15, height = 2.9, dpi = 300
  ),
  site_details = full_site_info(
    file_in(!!site_list_file),
    file_in("sites")
  ),
  site_lai_samples = target(
    predict_lai(site_details, tidy_posteriors),
    format = "fst_tbl"
  ),
  site_lai_summary = summarize_lai_samples(site_lai_samples),
  site_lai_total = site_lai_summary %>%
    group_by(site, year) %>%
    summarize(
      lai_mean = sum(lai_mean),
      lai_lo = sum(lai_lo),
      lai_hi = sum(lai_hi),
      elai_mean = sum(elai_mean),
      elai_lo = sum(elai_lo),
      elai_hi = sum(elai_hi)
    ) %>%
    ungroup(),
  lai_observed = file_in(!!fft_lai_file) %>%
    read_csv(col_types = cols(Site = "c", Site_Plot = "c", Subplot = "c",
                              .default = "n")) %>%
    group_by(site = Site_Plot) %>%
    summarize(
      obs_LAI = mean(LAI_CLX_10_60_m2_m2, na.rm = TRUE),
      obs_LAI_SD = sd(LAI_CLX_10_60_m2_m2, na.rm = TRUE),
      obs_LAI_lo = obs_LAI - obs_LAI_SD,
      obs_LAI_hi = obs_LAI + obs_LAI_SD
    ) %>%
    filter_if(is.numeric, all_vars(. > 0)),
  inversion_site_list = readLines(file_in(!!site_list_file)),
  predicted_spectra = target(
    predict_site_spectra(
      posterior_matrix,
      inversion_site_list,
      nsamp = 500, dedup = TRUE, progress = FALSE
    ) %>%
      dplyr::mutate(site = inversion_site_list),
    dynamic = map(inversion_site_list)
  ),
  observed_spectra = target(
    tidy_site_spec(inversion_site_list) %>%
      dplyr::mutate(site = inversion_site_list),
    dynamic = map(inversion_site_list)
  ),
  observed_predicted = observed_spectra %>%
    inner_join(predicted_spectra, c("site", "wavelength")) %>%
    mutate(bias = albedo_mean - observed),
  spec_error_all = ggsave(
    file_out(!!path(figdir, "spec-error-all.png")),
    spec_error_all_f(observed_predicted),
    width = 10, height = 14, dpi = 300
  ),
  spec_error_aggregate = ggsave(
    file_out(!!path(figdir, "spec-error-aggregate.png")),
    spec_error_aggregate_f(observed_predicted),
    width = 5, height = 4, dpi = 300
  ),
  overestimate_sites = ggsave(
    file_out(!!path(figdir, "overestimate-sites.png")),
    lapply(
      c("IDS35", "OF04", "SF03", "OF01", "GR08", "IDS34", "OF05"),
      site_spec_dbh_plot,
      observed_predicted = observed_predicted,
      site_details = site_details
    ) %>%
      wrap_plots(guides = "collect", ncol = 2) +
      guide_area(),
    width = 8, height = 10, dpi = 300
  )
)
