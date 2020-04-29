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
    dplyr::ungroup()
  site_df
}

predict_lai <- function(site_details, tidy_posteriors, max_samples = 5000) {
  site_dt <- site_details %>%
    dplyr::mutate(pft = factor(ipft, 1:5, c(
      "Early_Hardwood", "North_Mid_Hardwood", "Late_Hardwood",
      "Northern_Pine", "Late_Conifer"
    )))
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
    dplyr::left_join(site_dt, "pft") %>%
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
}
