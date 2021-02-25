plan <- drake_plan(
  traceplots = {
    samples_bt <- readRDS(file_in(!!pda_result_file))
    posterior <- BayesianTools::getSample(samples_bt, thin = "auto", coda = TRUE,
                                          start = pda_start)
    param_names <- readLines(file_in(!!param_names_file))
    param_traces <- param_names %>%
      str_remove("^temperate\\..*\\.") %>%
      str_subset("sitesoil", negate = TRUE) %>%
      str_remove("_(slope|intercept)") %>%
      unique()
    trace_file <- file_out(!!path(figdir, "traceplots.pdf"))
    pdf(trace_file, width = 8, height = 10)
    for (param in param_traces) {
      iparam <- str_detect(param_names, param)
      chains <- posterior[, iparam]
      plot(chains, density = FALSE, main = param)
    }
    dev.off()
  },
  traceplots_full = {
    samples_bt <- readRDS(file_in(!!pda_result_file))
    posterior <- BayesianTools::getSample(samples_bt, thin = "auto", coda = TRUE,
                                          start = 50000)
    param_names <- readLines(file_in(!!param_names_file))
    param_traces <- param_names %>%
      str_remove("^temperate\\..*\\.") %>%
      str_subset("sitesoil", negate = TRUE) %>%
      str_remove("_(slope|intercept)") %>%
      unique()
    trace_file <- file_out(!!path(figdir, "traceplots-full.pdf"))
    pdf(trace_file, width = 8, height = 10)
    for (param in param_traces) {
      iparam <- str_detect(param_names, param)
      chains <- posterior[, iparam]
      plot(chains, density = FALSE, main = param)
    }
    dev.off()
  },
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
  prior_posterior = {
    dat_prior <- tidy_priors %>%
      dplyr::mutate(pft = pft_factor(pft))
    dat_posterior <- tidy_posteriors %>%
      dplyr::mutate(pft = pft_factor(pft))
    ggsave(
      file_out(!!path(figdir, "posterior-pft.png")),
      pft_posterior_plot(dat_prior, dat_posterior, ncol = 3),
      width = 7, height = 5, dpi = 300
    )
    },
  prior_posterior_present = ggsave(
    file_out(!!path(figdir, "posterior-pft-presentation.png")),
    pft_posterior_plot(tidy_priors, tidy_posteriors, ncol = 4),
    width = 12, height = 6, dpi = 300
  ),
  site_structure_data = file_in(!!site_structure_file) %>%
    read_csv(col_types = "cnnnnnnclc") %>%
    mutate(site_tag = paste0("sitesoil_", row_number())) %>%
    dplyr::filter(site_name %in% readLines(here("other_site_data", "site_list"))),
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
  lai_pred_obs_plot = {
    site_lai_mod <- site_lai_total %>%
      left_join(biggest_pft, "site") %>%
      mutate(pft = fct_relabel(pft, ~str_replace_all(.x, "_" , " ")))
    ggsave(
      file_out(!!path(figdir, "lai-pred-obs.png")),
      lai_predicted_observed_plot(site_lai_mod, lai_observed),
      width = 6, height = 5, dpi = 300
    )
  },
  lai_bias_plots = {
    dat <- bias_data %>%
      dplyr::filter(band == band[1])
    p1 <- ggplot(dat) +
      aes(x = mean_dbh, y = LAI_diff) +
      geom_smooth(method = "lm") +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Mean DBH (cm)", y = "LAI residual (predicted - observed)") +
      theme_bw()
    ggsave(file_out(!!path(figdir, "lai-bias-dbh.png")), p1,
           width = 4, height = 4, units = "in", dpi = 300)
    p2 <- p1 + aes(x = tot_dens) + labs(x = expression("Stand density" ~ (trees ~ ha^-1)))
    ggsave(file_out(!!path(figdir, "lai-bias-dens.png")), p2,
           width = 4, height = 4, units = "in", dpi = 300)
    p2bypft <- p2 %+% facet_wrap(vars(pft))
    ggsave(file_out(!!path(figdir, "lai-bias-dens-bypft.png")), p2bypft,
           width = 6, height = 4, units = "in", dpi = 300)
    p3 <- p1 + aes(x = frac_evergreen_wtd) + labs(x = "Evergreen fraction")
    ggsave(file_out(!!path(figdir, "lai-bias-evergreen.png")), p3,
           width = 4, height = 4, units = "in", dpi = 300)
    p4 <- ggplot(dplyr::filter(bias_data, pft != "All")) +
      aes(x = mean_diff, y = LAI_diff) +
      geom_smooth(method = "lm", color = "black") +
      geom_point(aes(color = pft)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      facet_grid(cols = vars(band), scales = "free_x") +
      labs(x = "Reflectance bias (predicted - observed)",
           y = "LAI residual (predicted - observed)") +
      scale_color_brewer(palette = "Set1") +
      theme_bw()
    ggsave(file_out(!!path(figdir, "lai-bias-refl-bias.png")), p4,
           width = 6, height = 4, units = "in", dpi = 300)
    p5 <- p1 + facet_wrap(vars(pft))
    ggsave(file_out(!!path(figdir, "lai-bias-dbh-bypft.png")), p5,
           width = 6, height = 4, units = "in", dpi = 300)
  },
  inversion_site_list = readLines(file_in(!!site_list_file)),
  edr_site_inputs = {
    aviris <- aviris_data() %>% select(-starts_with("band_"))
    left_join(site_details, aviris, c("site" = "iPLOT"))
  },
  predicted_spectra = target(
    predict_site_spectra(posterior_matrix, inversion_site_list, edr_site_inputs,
                         nsamp = 500, dedup = TRUE),
    dynamic = map(inversion_site_list)
  ),
  observed_spectra = {
    aviris_long <- aviris_data() %>%
      pivot_longer(starts_with("band_"), names_to = "band", values_to = "observed")
    aviris_waves <- read_csv(here("aviris", "NASA_FFT", "aviris_c_wavelength.csv"),
                             col_types = "n") %>%
      arrange(wavelength) %>%
      mutate(band = paste0("band_", row_number()))
    observed_spectra <- aviris_long %>%
      left_join(aviris_waves, "band") %>%
      select(aviris_id, site = iPLOT, czen, direct_sky_frac,
             wavelength, observed, everything())
    observed_spectra
  },
  observed_predicted = {
    pred2 <- predicted_spectra %>%
      dplyr::select(aviris_id, site, result) %>%
      tidyr::unnest(result)
    pred_waves <- unique(pred2$wavelength)
    resamp <- function(x, y, xout) approx(x, y, xout)$y
    obs_sampled <- observed_spectra %>%
      select(aviris_id, site, wavelength, observed) %>%
      group_by(aviris_id, site) %>%
      summarize(
        observed = list(approx(wavelength, observed, pred_waves, rule = 2)$y),
        wavelength = list(pred_waves)
      ) %>%
      ungroup() %>%
      unnest(c(wavelength, observed))
    obs_sampled %>%
      inner_join(pred2, c("site", "wavelength", "aviris_id")) %>%
      mutate(bias = albedo_mean - observed)
  },
  spec_error_all = {
    plot_dat <- observed_predicted %>%
      dplyr::group_by(site) %>%
      dplyr::mutate(site_mean_bias = mean(bias)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(site_f = forcats::fct_reorder(site, site_mean_bias))
    plt <- ggplot(plot_dat) +
      aes(x = wavelength, group = aviris_id) +
      geom_ribbon(aes(ymin = pmax(albedo_r_q025, 0),
                      ymax = pmin(albedo_r_q975, 1),
                      fill = "95% PI")) +
      geom_ribbon(aes(ymin = pmax(albedo_q025, 0),
                      ymax = pmin(albedo_q975, 1),
                      fill = "95% CI")) +
      geom_line(aes(y = albedo_mean, color = "EDR"), size = 1) +
      geom_line(aes(y = observed, color = "AVIRIS")) +
      facet_wrap(vars(site_f), scales = "fixed", ncol = 6) +
      labs(x = "Wavelength (nm)", y = "Reflectance (0 - 1)") +
      scale_fill_manual(
        name = "",
        values = c("95% PI" = "gray80",
                   "95% CI" = "green3")
      ) +
      scale_color_manual(
        name = "",
        values = c("EDR" = "green4",
                   "AVIRIS" = "black")
      ) +
      theme_bw()
    ggsave(
      file_out(!!path(figdir, "spec-error-all.png")), plt,
      width = 10, height = 14, units = "in", dpi = 300
    )
  },
  spec_error_aggregate = {
    biggest_pft2 <- biggest_pft %>%
      mutate(pft = fct_relabel(pft, ~str_replace_all(.x, "_", " ")))
    plot_data <- observed_predicted %>%
      left_join(biggest_pft2, "site") %>%
      arrange(pft) %>%
      mutate(pft = as.character(pft)) %>%
      bind_rows(mutate(observed_predicted, pft = "All sites"), .) %>%
      mutate(pft = fct_inorder(pft))
    nsites <- biggest_pft2 %>%
      count(pft) %>%
      mutate(pft = as.character(pft)) %>%
      add_row(pft = "All sites", n = nrow(biggest_pft2), .before = 1) %>%
      mutate(pft = fct_inorder(pft),
             label = paste("N[site] ==", n))
    avg <- plot_data %>%
      group_by(pft, wavelength) %>%
      summarize(
        bias_lo = quantile(bias, 0.25),
        bias_mid = median(bias),
        bias_hi = quantile(bias, 0.75)
      ) %>%
      ungroup()
    plt <- ggplot(plot_data) +
      aes(x = wavelength) +
      geom_ribbon(aes(ymin = bias_lo, ymax = bias_hi, y = NULL),
                  fill = "deepskyblue", color = "deepskyblue", data = avg) +
      geom_line(aes(y = bias, group = aviris_id), color = "black", size = 0.2, alpha = 0.3) +
      geom_line(aes(y = bias_mid), color = "black", data = avg) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(750, 1100), linetype = "dashed", size = 0.2) +
      geom_text(aes(x = -Inf, y = -Inf, label = label), data = nsites,
                parse = TRUE, hjust = -0.1, vjust = -0.1) +
      facet_wrap(vars(pft)) +
      labs(x = "Wavelength (nm)",
           y = "Predicted (mean) - observed reflectance") +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      scale_x_continuous(breaks = c(400, 600, 750, 900, 1100, 1300))
    ggsave(
      file_out(!!path(figdir, "spec-error-aggregate.png")),
      plt,
      width = 8, height = 5, dpi = 300
    )
  },
  spec_error_binned = {
    biggest_pft2 <- biggest_pft %>%
      mutate(pft = fct_relabel(pft, ~str_replace_all(.x, "_", " ")))
    bin <- 50
    plot_data <- observed_predicted %>%
      group_by(site, aviris_id, wave_band = round(wavelength / bin) * bin) %>%
      summarize(bias = mean(bias, na.rm = TRUE)) %>%
      ungroup() %>%
      left_join(biggest_pft2, "site")
    nsites <- biggest_pft2 %>%
      count(pft) %>%
      mutate(label = paste("N[site] ==", n))
    plt <- ggplot(plot_data) +
      aes(x = wave_band, y = bias, group = wave_band) +
      geom_jitter(color = "gray70", size = 0.4, alpha = 0.6) +
      geom_boxplot(outlier.shape = NA, position = "identity", fill = NA, weight = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_text(aes(x = -Inf, y = -Inf, label = label, group = NULL), data = nsites,
                parse = TRUE, hjust = -0.1, vjust = -0.1) +
      facet_wrap(vars(pft)) +
      labs(x = "Wavelength (nm)",
           y = "Predicted (mean) - observed reflectance") +
      theme_bw() +
      theme(panel.grid = element_blank())
    ggsave(
      file_out(!!path(figdir, "spec-error-binned.png")),
      plt,
      width = 8, height = 5, dpi = 300
    )
  },
  tree_sites_q = target({
    sites <- site_details %>%
      dplyr::group_by(site) %>%
      dplyr::summarize(max_dbh = max(dbh)) %>%
      dplyr::filter(max_dbh >= quantile(max_dbh, qlo),
                    max_dbh < quantile(max_dbh, qlo + 0.25)) %>%
      dplyr::arrange(max_dbh) %>%
      dplyr::pull(site)
    ymax <- observed_predicted %>%
      dplyr::filter(site %in% sites) %>%
      dplyr::select(observed, starts_with("albedo")) %>%
      as.matrix() %>%
      max(na.rm = TRUE)
    plt <- lapply(sites, site_spec_dbh_plot,
                  observed_predicted = observed_predicted,
                  site_details = site_details,
                  ymax = ymax) %>%
      wrap_plots(guides = "collect", ncol = 3)
    f <- path(figdir, paste0("tree-sites-q", qlo * 100, ".png"))
    ggsave(f, plt, width = 12, height = 8, dpi = 300)
  }, transform = map(qlo = c(0, 0.25, 0.50, 0.75))),
  tree_sites_bypft = target({
    pp <- stringr::str_remove(PFT, "temperate\\.")
    dat <- site_details %>%
      dplyr::group_by(site, pft) %>%
      dplyr::summarize(pft_dbh = sum(dbh * nplant)) %>%
      dplyr::filter(pft_dbh == max(pft_dbh)) %>%
      dplyr::ungroup()
    sites <- dat %>%
      dplyr::filter(pft == pp) %>%
      dplyr::pull(site)
    ymax <- observed_predicted %>%
      dplyr::filter(site %in% sites) %>%
      dplyr::select(observed, starts_with("albedo")) %>%
      as.matrix() %>%
      max(na.rm = TRUE)
    plt <- lapply(sites, site_spec_dbh_plot,
                  observed_predicted = observed_predicted,
                  site_details = site_details,
                  ymax = ymax) %>%
      wrap_plots(guides = "collect", ncol = 3)
    f <- path(figdir, paste0("tree-sites-", pp, ".png"))
    ggsave(f, plt, width = 14, height = 8, dpi = 300)
  }, transform = map(PFT = !!pfts)),
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
  site_map = {
    site_structure_sf <- sf::st_as_sf(site_structure_data,
                                      coords = c("longitude", "latitude"),
                                      crs = 4326)
    bbox <- sf::st_bbox(site_structure_sf) %>%
      sf::st_as_sfc() %>%
      sf::st_transform(5070) %>%
      sf::st_buffer(100 * 1000) %>%
      sf::st_bbox()
    usa <- rnaturalearth::ne_states(country = c("United States of America", "Canada"),
                                    returnclass = "sf")

    basemap_file <- here("data-raw/CONUS_Natural_Color_Relief_Hydro.tif")
    stopifnot(file.exists(basemap_file))
    basemap <- raster::brick(basemap_file) %>%
      raster::crop(bbox)

    site_structure_selected <- site_structure_sf %>%
      dplyr::filter(site_name %in% spec_summary_sites)

    plot_map <- ggplot() +
      ggspatial::layer_spatial(data = basemap) +
      ggspatial::layer_spatial(usa, fill = NA) +
      ggspatial::layer_spatial(site_structure_sf, size = 1.2) +
      ggrepel::geom_label_repel(
        data = site_structure_selected,
        aes(label = site_name, geometry = geometry),
        size = 2,
        stat = "sf_coordinates",
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      coord_sf(xlim = bbox[c(1,3)], ylim = bbox[c(2,4)], crs = 5070, expand = FALSE) +
      theme_bw() +
      ggspatial::annotation_north_arrow(location = "tr") +
      ggspatial::annotation_scale(location = "bl") +
      theme(plot.background = element_blank(),
            axis.title = element_blank())
    ggsave(path(figdir, "site-map.png"), plot_map,
           width = 6, height = 4, units = "in", dpi = 300)

  },
  site_structure = {
    site_structure_data_plt <- site_structure_data %>%
      dplyr::left_join(biggest_pft, c("site_name" = "site"))
    site_structure_sub <- site_structure_data_plt %>%
      dplyr::filter(site_name %in% spec_summary_sites)

    plot_travis <- ggplot(site_structure_data_plt) +
      aes(x = mean_dbh, y = tot_dens * 5000) +
      geom_point(aes(color = pft)) +
      # Self-thinning curve
      geom_function(
        aes(linetype = "y == 500 * bgroup('(', frac(x, 25), ')')^-1.4"),
        fun = ~500 * (.x/25)^-1.4,
        key_glyph = draw_key_abline
      ) +
      scale_linetype_manual(values = "dashed", labels = scales::parse_format()) +
      scale_color_brewer(palette = "Set1") +
      ggrepel::geom_text_repel(data = site_structure_sub, aes(label = site_name),
                                min.segment.length = 0, size = 2.5) +
      labs(x = "Mean diameter (cm)",
           y = expression("Stem density" ~ (trees ~ ha^-1))) +
      coord_cartesian(ylim = range(site_structure_data$tot_dens * 5000)) +
      theme_bw() +
      theme(legend.position = c(1, 1),
            legend.justification = c(1, 1),
            legend.title = element_blank(),
            legend.background = element_blank())
    if (interactive()) plot_travis

    ggsave(path(figdir, "site-structure.png"), plot_travis,
           width = 5, height = 5, units = "in", dpi = 300)
  },
  site_structure_map = {
    site_structure_sf <- sf::st_as_sf(site_structure_data,
                                      coords = c("longitude", "latitude"),
                                      crs = 4326)

    bbox <- sf::st_bbox(site_structure_sf) %>%
      sf::st_as_sfc() %>%
      sf::st_transform(5070) %>%
      sf::st_bbox()
    usa <- rnaturalearth::ne_states(country = c("United States of America", "Canada"),
                                    returnclass = "sf")

    plot_map <- ggplot() +
      geom_sf(data = usa, fill = NA) +
      geom_sf(data = site_structure_sf, size = 1.2) +
      coord_sf(xlim = bbox[c(1,3)], ylim = bbox[c(2,4)], crs = 5070) +
      theme_bw() +
      theme(plot.background = element_blank())
    plot_travis <- ggplot(site_structure_data) +
      aes(x = mean_dbh, y = tot_dens * 5000) +
      geom_point() +
      labs(x = "Mean diameter (cm)",
           y = expression("Stem density" ~ (trees ~ ha^-1))) +
      theme_bw()
    layout <- c(
      area(1, 1, 100, 100),
      area(2, 40, 50, 98)
    )
    finalplot <- plot_travis + plot_map + plot_layout(design = layout)
    ggsave(
      file_out(!!path(figdir, "site-map-structure.png")),
      finalplot,
      width = 6, height = 4, units = "in", dpi = 300
    )
  },
  spec_summary_sites = c("NC18", "NC07", "NC12",
                         "IDS10", "BI01", "BH10",
                         "SF04", "AK60", "BI05"),
  spec_summary_plot = {
    sites <- spec_summary_sites
    plt <- lapply(
      sites, site_spec_dbh_plot,
      observed_predicted = observed_predicted,
      site_details = site_details,
      spec_additions = list(
        coord_cartesian(ylim = c(0, 0.65)),
        theme(axis.title = element_blank())
      ),
      dbh_additions = list(
        coord_cartesian(xlim = c(0, 75)),
        theme(axis.title = element_blank())
      )
    ) %>%
      wrap_plots(guides = "collect", ncol = 3)
    gt <- patchwork::patchworkGrob(plt + theme_bw())
    plt2 <- gridExtra::arrangeGrob(
      gt,
      left = "Left: Reflectance [0, 1]    Right: Count",
      bottom = "Left: Wavelength (nm)    Right: Stem diameter (cm)"
    )
    ggsave(
      file_out(!!path(figdir, "mainfig-sites.png")), plt2,
      width = 12, height = 5, dpi = 300
    )
  },
  biggest_pft = site_details %>%
    dplyr::group_by(site, pft) %>%
    dplyr::summarize(pft_dbh = sum(dbh * nplant)) %>%
    dplyr::filter(pft_dbh == max(pft_dbh)) %>%
    dplyr::ungroup(),
  bias_data = {
    band_diffs <- observed_predicted %>%
      dplyr::mutate(
        band = cut(wavelength, c(400, 750, 1100, 1300),
                   include.lowest = TRUE, dig.lab = 4) %>%
          fct_relabel(~paste(.x, "nm"))
      ) %>%
      dplyr::group_by(band, site) %>%
      dplyr::summarize(mean_diff = mean(bias, na.rm = TRUE)) %>%
      dplyr::ungroup()
    lai_both <- site_lai_total %>%
      dplyr::left_join(lai_observed, "site") %>%
      dplyr::mutate(LAI_diff = lai_mean - obs_LAI) %>%
      dplyr::select(site, LAI_diff, lai_mean, obs_LAI, dplyr::everything()) %>%
      dplyr::arrange(dplyr::desc(LAI_diff))
    bias_data_0 <- band_diffs %>%
      dplyr::left_join(site_structure_data, c("site" = "site_name")) %>%
      dplyr::left_join(lai_both, "site") %>%
      dplyr::mutate(mostly_evergreen = (frac_evergreen > 0.5) %>%
                      factor(labels = paste("Mostly", c("deciduous", "evergreen"))),
                    tot_dens = tot_dens * 5000)
    bias_data_pft <- bias_data_0 %>%
      dplyr::left_join(biggest_pft, "site")
    bias_data <- bias_data_0 %>%
      dplyr::mutate(pft = "All") %>%
      dplyr::bind_rows(bias_data_pft) %>%
      dplyr::mutate(pft = factor(pft, labels = c("All", "EH", "MH", "LH", "NP", "LC")))
    bias_data
  },
  supplement_pft_plots = {
    pft_boxplot <- ggplot(bias_data) +
      aes(x = pft, y = mean_diff) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(color = "gray50", width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(vars(band), scales = "free_y") +
      labs(x = "Plant functional type",
           y = "Mean reflectance bias (predicted - observed)") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5))
    ggsave(file_out(!!path(figdir, "bias-boxplot-pft.png")), pft_boxplot,
           width = 7, height = 3.75, units = "in", dpi = 300)
    fits <- bias_data %>%
      dplyr::group_by(band, mostly_evergreen) %>%
      dplyr::summarize(fit = list(lm(mean_diff ~ tot_dens))) %>%
      dplyr::mutate(
        intercept = purrr::map_dbl(fit, ~.x$coefficients[1]),
        slope = purrr::map_dbl(fit, ~.x$coefficients[2]),
        r2 = purrr::map_dbl(fit, ~summary(.x)$adj.r.squared),
        p = purrr::map_dbl(fit, ~summary(.x)$coefficients[2, 4])
      )
    pft_density_bias_plot <- ggplot(bias_data) +
      aes(x = tot_dens, y = mean_diff) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(vars(band), vars(pft), scales = "free_y") +
      labs(x = expression("Stem density" ~ (trees ~ ha^-1)),
           y = "Mean reflectance bias (predicted - observed)") +
      theme_bw() +
      theme(panel.grid = element_blank())
    ggsave(file_out(!!path(figdir, "bias-density-pft.png")),
           pft_density_bias_plot,
           width = 7, height = 4.5, units = "in", dpi = 300)
  },
  lai_sens_edr = sensitivity_plot(c(seq(0.2, 1, 0.2), seq(1, 5, 0.5), seq(5, 10, 1)), "lai", "LAI"),
  clumping_sens_edr = sensitivity_plot(seq(0.1, 1.0, 0.1), "clumping_factor", "Clumping factor"),
  orient_sens_edr = sensitivity_plot(seq(-0.4, 0.6, 0.1), "orient_factor", "Leaf orient"),
  N_sens_edr = sensitivity_plot(seq(1.1, 3.0, 0.2), "N", "PROSPECT N"),
  czen_sens_edr = sensitivity_plot(seq(0.5, 1.0, 0.05), "czen", expression(cos(theta))),
  dsf_sens_edr = sensitivity_plot(seq(0.0, 1.0, 0.05), "direct_sky_frac", expression(f[direct])),
  lai_cohort_sens_edr = {
    # All of these sum to LAI=4, but distributed differently across cohorts.
    # Question is: Does EDR produce different results depending on the number of
    # cohorts, even if the total LAI is the same?
    #
    # Answer: No, it does not -- results are the same.
    values <- list(
      4,
      c(3, 1),
      c(2, 1, 0.5, 0.3, 0.2)
    )
    varname <- "lai"
    label <- "# cohorts"
    sens <- purrr::map(
      values, do_sens,
      fun = edr_r,
      variable = varname,
      .dots = edr_sensitivity_defaults
    ) %>%
      tidy_albedo(seq_along(values))
    plt <- ggplot(sens) +
      aes(x = wavelength, y = value, color = factor(var_value),
          group = var_value) +
      geom_line() +
      scale_color_viridis_d() +
      labs(x = "Wavelength (nm)", y = "Albedo [0,1]",
           color = label) +
      theme_bw() +
      theme(
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank()
      )
    ggsave(path(figdir, "edr-sensitivity-cohort_lai.png"), plt,
           width = 4, height = 4, dpi = 300)
  },
  edr_sail_comparison_lai = {
    lai <- c(
      0.1, 0.5,
      seq(1, 5, 1)
    )
    lai_dsf <- 0.9
    lai_sens_edr <- purrr::map(
      lai, do_sens,
      fun = edr_r,
      variable = "lai",
      .dots = modifyList(edr_sensitivity_defaults, list(direct_sky_frac = lai_dsf))
    ) %>%
      tidy_albedo(lai)
    lai_sens_sail <- purrr::map(
      lai, do_sens, fun = rrtm::pro4sail_5,
      variable = "LAI",
      .dots = sail_sensitivity_defaults
    ) %>%
      tidy_sail(lai)
    tidy_both <- dplyr::left_join(lai_sens_edr, lai_sens_sail) %>%
      dplyr::mutate(
        sail_dr = bdr * lai_dsf + hdr * (1 - lai_dsf),
        sail_hr = dhr * lai_dsf + bhr * (1 - lai_dsf)
      )
    tidy_both_long <- tidy_both %>%
      tidyr::pivot_longer(c(value, bhr:sail_hr))

    plot_dat <- tidy_both_long %>%
      dplyr::rename(LAI = var_value) %>%
      dplyr::mutate(Source = factor(name, c("value", "sail_hr", "sail_dr", "bhr", "dhr", "hdr", "bdr"),
                                    c("EDR", "SAIL: HR", "SAIL: DR", "SAIL: BHR", "SAIL: DHR", "SAIL: HDR", "SAIL: BDR"))) %>%
      dplyr::group_by(LAI, Source, band = band_cut(wavelength)) %>%
      dplyr::summarize(value = mean(value)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(band))
    plt <- ggplot(plot_dat) +
      aes(x = LAI, y = value, color = Source, linetype = Source) +
      geom_line() +
      scale_color_brewer(palette = "Dark2") +
      scale_linetype_manual(values = c(1, 1, 1, 3, 3, 3, 3)) +
      facet_wrap(vars(band), scales = "free_y") +
      labs(x = "LAI", y = "Reflectance [0,1]", color = "Source") +
      theme_bw() +
      theme(legend.position = "bottom")
    ggsave(path(figdir, "edr-sail-comparison-lai.png"), plt,
           width = 6, height = 3.8, units = "in", dpi = 300)
  },
  edr_sail_comparison_czen = {
    czen <- seq(0.5, 1.0, 0.1)
    dsf <- 0.9
    sens_edr <- purrr::map(
      czen, do_sens,
      fun = edr_r,
      variable = "czen",
      .dots = modifyList(edr_sensitivity_defaults, list(direct_sky_frac = dsf))
    ) %>%
      tidy_albedo(czen)
    zen <- acos(czen) * 180/pi
    sens_sail <- purrr::map(
      zen, do_sens, fun = rrtm::pro4sail_5,
      variable = "solar_zenith",
      .dots = sail_sensitivity_defaults
    ) %>%
      tidy_sail(czen)
    tidy_both <- dplyr::left_join(sens_edr, sens_sail) %>%
      dplyr::mutate(
        sail_dr = bdr * dsf + hdr * (1 - dsf),
        sail_hr = dhr * dsf + bhr * (1 - dsf)
      )
    tidy_both_long <- tidy_both %>%
      tidyr::pivot_longer(c(value, bhr:sail_hr))

    plot_dat <- tidy_both_long %>%
      dplyr::rename(czen = var_value) %>%
      dplyr::mutate(Source = factor(name, c("value", "sail_hr", "sail_dr", "bhr", "dhr", "hdr", "bdr"),
                                    c("EDR", "SAIL: HR", "SAIL: DR", "SAIL: BHR", "SAIL: DHR", "SAIL: HDR", "SAIL: BDR"))) %>%
      dplyr::group_by(czen, Source, band = band_cut(wavelength)) %>%
      dplyr::summarize(value = mean(value)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(band))
    plt <- ggplot(plot_dat) +
      aes(x = czen, y = value, color = Source, linetype = Source) +
      geom_line() +
      scale_color_brewer(palette = "Dark2") +
      scale_linetype_manual(values = c(1, 1, 1, 3, 3, 3, 3)) +
      facet_wrap(vars(band), scales = "free_y") +
      labs(x = expression(cos(theta[s])), y = "Reflectance [0,1]", color = "Source") +
      theme_bw() +
      theme(legend.position = "bottom")
    ggsave(path(figdir, "edr-sail-comparison-czen.png"), plt,
           width = 6, height = 3.8, units = "in", dpi = 300)
  },
  posterior_correlations = {
    cormat <- cor(posterior_matrix)
    wide_cor <- corrr::as_cordf(cormat)
    long_cor <- wide_cor %>%
      corrr::shave() %>%
      corrr::stretch(na.rm = TRUE) %>%
      dplyr::mutate(dplyr::across(c(x, y), ~gsub("temperate\\.", "", .x)))

    lc2 <- long_cor %>%
      tidyr::separate(x, c("xpft", "xvar"), sep = "\\.") %>%
      tidyr::separate(y, c("ypft", "yvar"), sep = "\\.") %>%
      dplyr::mutate(
        dplyr::across(c(xvar, yvar), ~gsub("prospect_", "", .x)),
        dplyr::across(c(xvar, yvar), ~gsub("_factor", "", .x))
      )

    ## lc2 %>%
    ##   dplyr::filter(dplyr::across(c(xpft, ypft), ~!grepl("sitesoil", .x))) %>%
    ##   dplyr::arrange(dplyr::desc(r)) %>%
    ##   print(n = 20)

    ## lc2 %>%
    ##   ## dplyr::filter(xvar == "clumping", xpft == "North_Mid_Hardwood") %>%
    ##   dplyr::filter(xvar == "orient", xpft == "Late_Hardwood") %>%
    ##   dplyr::arrange(dplyr::desc(abs(r)))

    lc3 <- lc2 %>%
      dplyr::filter(xpft == ypft) %>%
      dplyr::mutate(pft = xpft)

    plot_dat <- lc3 %>%
      dplyr::filter(!is.na(pft), !grepl("sitesoil", pft),
                    !is.na(xvar), !is.na(yvar)) %>%
      dplyr::mutate(dplyr::across(c(pft, xvar, yvar), forcats::fct_inorder))

    plt <- ggplot(plot_dat) +
      aes(x = 1, y = r, fill = pft) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_col(position = "dodge") +
      facet_grid(vars(xvar), vars(yvar), drop = TRUE) +
      scale_fill_brewer(palette = "Set1", name = "") +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = c(0, 0),
            legend.justification = c(0, 0))

    ggsave(file_out(!!path(figdir, "posterior-correlations.png")), plt,
           width = 7, height = 7, dpi = 300, units = "in")
  }
)
