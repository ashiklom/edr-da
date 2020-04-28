plan <- drake_plan(
  pda_result_file = last_result_file(),
  tidy_posteriors = target({
    samples_bt <- readRDS(pda_result_file)
    posterior_matrix <- BayesianTools::getSample(
      samples_bt,
      thin = "auto",
      start = 15000
    )
    colnames(posterior_matrix) <- readLines(file_in("param_names.txt"))
    tidy_param_matrix(posterior_matrix, "posterior")
  }, format = "fst_tbl"),
  tidy_priors = target(tidy_prior(), format = "fst_tbl"),
  prior_posterior = ggsave(
    file_out(!!path(figdir, "posterior-pft.png")),
    pft_posterior_plot(tidy_priors, tidy_posteriors),
    width = 6, height = 7
  ),
  site_structure_data = file_in(!!site_structure_file) %>%
    read_csv(col_types = "cnnnnnnclc") %>%
    mutate(site_tag = paste0("sitesoil_", row_number())),
  soil_moisture_plt = ggsave(
    file_out(!!path(figdir, "posterior-soil.png")),
    soil_moisture_plot(tidy_posteriors, site_structure_data),
    width = 15, height = 2.9
  ),
  site_details = full_site_info(
    file_in(!!site_list_file),
    file_in("sites")
  ),
  site_lai_samples = target(
    predict_lai(site_details, tidy_posteriors),
    format = "fst_dt"
  ),
  site_lai_summary = summarize_lai_samples(site_lai_samples)
)
