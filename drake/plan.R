plan <- drake_plan(
  pda_result_file = last_result_file(),
  samples_bt = readRDS(pda_result_file),
  posterior_matrix = BayesianTools::getSample(
    samples_bt,
    thin = "auto",
    start = 15000
  ) %>%
    `colnames<-`(readLines(file_in("param_names.txt"))),
  tidy_posteriors = tidy_param_matrix(posterior_matrix, "posterior"),
  tidy_priors = tidy_prior(),
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
  )
)
