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
    file_out("text/figures/pda-summary-pft.png"),
    prior_posterior_plot(tidy_priors, tidy_posteriors),
    width = 6, height = 7
  )
)
