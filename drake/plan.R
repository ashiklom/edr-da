plan <- drake_plan(
  posterior_matrix = {
    samples_bt <- readRDS(file_in(!!pda_result_file))
    posterior_matrix <- BayesianTools::getSample(
      samples_bt,
      thin = "auto",
      start = !!pda_start
    )
    colnames(posterior_matrix) <- readLines(file_in(!!param_names_file))
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
  lai_pred_obs_plot = ggsave(
    file_out(!!path(figdir, "lai-pred-obs.png")),
    lai_predicted_observed_plot(site_lai_total, lai_observed),
    width = 6, height = 5, dpi = 300
  ),
  inversion_site_list = readLines(file_in(!!site_list_file)),
  predicted_spectra = target(
    predict_site_spectra(
      posterior_matrix,
      inversion_site_list,
      nsamp = 500, dedup = TRUE, progress = FALSE,
      run_config = !!run_config
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
    spec_error_all_f(observed_predicted, sail_predictions),
    width = 10, height = 14, dpi = 300
  ),
  spec_error_aggregate = ggsave(
    file_out(!!path(figdir, "spec-error-aggregate.png")),
    spec_error_aggregate_f(observed_predicted),
    width = 5, height = 4, dpi = 300
  ),
  small_tree_sites = {
    sites <- site_details %>%
      dplyr::group_by(site) %>%
      dplyr::summarize(max_dbh = max(dbh)) %>%
      dplyr::filter(max_dbh < 27.5) %>%
      dplyr::pull(site)
    plt <- lapply(sites, site_spec_dbh_plot,
                  observed_predicted = observed_predicted,
                  site_details = site_details) %>%
      wrap_plots(guides = "collect", ncol = 3) +
      guide_area()
    ggsave(
      file_out(!!path(figdir, "small-tree-sites.png")),
      plt,
      width = 12, height = 8, dpi = 300
    )
  },
  med_tree_sites = {
    sites <- site_details %>%
      dplyr::group_by(site) %>%
      dplyr::summarize(max_dbh = max(dbh)) %>%
      dplyr::filter(max_dbh >= quantile(max_dbh, 0.25) &
                      max_dbh <= quantile(max_dbh, 0.75)) %>%
      dplyr::pull(site)
    plt <- lapply(sites, site_spec_dbh_plot,
                  observed_predicted = observed_predicted,
                  site_details = site_details) %>%
      wrap_plots(guides = "collect", ncol = 3) +
      guide_area()
    ggsave(
      file_out(!!path(figdir, "med-tree-sites.png")),
      plt,
      width = 12, height = 8, dpi = 300
    )
  },
  big_tree_sites = {
    sites <- site_details %>%
      dplyr::group_by(site) %>%
      dplyr::summarize(max_dbh = max(dbh)) %>%
      dplyr::filter(max_dbh > quantile(max_dbh, 0.75)) %>%
      dplyr::pull(site)
    plt <- lapply(sites, site_spec_dbh_plot,
                  observed_predicted = observed_predicted,
                  site_details = site_details) %>%
      wrap_plots(guides = "collect", ncol = 3) +
      guide_area()
    ggsave(
      file_out(!!path(figdir, "big-tree-sites.png")),
      plt,
      width = 12, height = 8, dpi = 300
    )
  },
  underpredict_sites = ggsave(
    file_out(!!path(figdir, "underpredict-sites.png")),
    lapply(
      c("IDS34", "OF01", "OF05", "OF04", "IDS40",
        "OF02", "AK06", "GR08", "AK60", "IDS36"),
      site_spec_dbh_plot,
      observed_predicted = observed_predicted,
      site_details = site_details
    ) %>%
      wrap_plots(guides = "collect", ncol = 3) +
      guide_area(),
    width = 12, height = 8, dpi = 300
  ),
  overpredict_sites = ggsave(
    file_out(!!path(figdir, "overpredict-sites.png")),
    lapply(
      c("NC22", "MN02", "NC17", "NC10", "MN04", "NC21"),
      site_spec_dbh_plot,
      observed_predicted = observed_predicted,
      site_details = site_details
    ) %>%
      wrap_plots(guides = "collect", ncol = 3) +
      guide_area(),
    width = 12, height = 8, dpi = 300
  ),
  all_sites_spec_dbh = ggsave(
    file_out(!!path(figdir, "all-sites-spec-dbh.png")),
    lapply(
      sort(unique(site_details$site)),
      site_spec_dbh_plot,
      observed_predicted = observed_predicted,
      site_details = site_details
    ) %>%
      wrap_plots(guides = "collect", ncol = 2),
    width = 8, height = 50, dpi = 300,
    limitsize = FALSE
  ),
  both_ndvi = calc_ndvi_bysite(observed_spectra, predicted_spectra,
                               site_structure_data),
  ndvi_dbh_png = ggsave(
    file_out(!!path(figdir, "ndvi-dbh.png")),
    ndvi_dbh_plot(both_ndvi),
    width = 6.4, height = 5.2, dpi = 300
  ),
  sail_predictions = tidy_sail_predictions(site_details, site_lai_total,
                                           tidy_posteriors)
)
