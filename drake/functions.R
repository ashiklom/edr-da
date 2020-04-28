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
