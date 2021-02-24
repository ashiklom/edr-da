# Interactive code for evaluating relative constraints

drake::loadd(tidy_posteriors)
drake::loadd(tidy_priors)

tidy_summary <- . %>%
  dplyr::group_by(biome, pft, variable, type) %>%
  dplyr::summarize(
    q025 = quantile(value, 0.025),
    q975 = quantile(value, 0.975),
    ci95 = q975 - q025
  ) %>%
  dplyr::ungroup()

priors <- tidy_priors %>% tidy_summary
posteriors <- tidy_posteriors %>% tidy_summary

both <- dplyr::bind_rows(priors, posteriors)

both %>%
  dplyr::arrange(biome, pft, variable, type)

both_wide <- both %>%
  dplyr::select(-q025, -q975) %>%
  tidyr::pivot_wider(names_from = type, values_from = ci95) %>%
  dplyr::mutate(rel_reduction = posterior / prior)

both_sub <- both_wide %>%
  dplyr::filter(!grepl("sitesoil|residual", variable))

# All parameters
both_sub %>%
  dplyr::arrange(rel_reduction) %>%
  print(n = Inf)

# Overall constraint
both_sub %>%
  dplyr::summarize(reduction = mean(rel_reduction))

# By variable
both_sub %>%
  dplyr::group_by(variable) %>%
  dplyr::summarize(reduction = mean(rel_reduction)) %>%
  dplyr::arrange(reduction)

# By PFT
both_sub %>%
  dplyr::group_by(pft) %>%
  dplyr::summarize(reduction = mean(rel_reduction)) %>%
  dplyr::arrange(reduction)
