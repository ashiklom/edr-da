library(tidyverse)
library(redr)
library(cowplot)
import::from(here, inhere = here)
import::from(fs, path_file, dir_ls)
import::from(glue, glue)
import::from(progress, progress_bar)
import::from(viridis, scale_color_viridis, scale_fill_viridis)
import::from(broom, tidy)
import::from(ggrepel, geom_text_repel, geom_label_repel)

data_dir <- inhere("sync_data", "msp20180402")

site_list <- readLines(inhere("other_site_data", "site_list"))
nsite <- length(site_list)

spsummary <- function(mat, ...) {
  if (PEcAnRTM::is_spectra(mat)) {
    waves <- PEcAnRTM::wavelengths(mat)
  } else {
    waves <- 400:2500
  }
  tibble(
    waves = waves,
    ...,
    mu = rowMeans(mat, na.rm = TRUE),
    lo = apply(mat, 1, quantile, 0.025, na.rm = TRUE),
    hi = apply(mat, 1, quantile, 0.975, na.rm = TRUE)
  )
}

get_site_data <- function(site) {
  spec_history_file <- dir_ls(data_dir, glob = glue("*spec_history.{site}*.rds"))
  prior_file <- dir_ls(data_dir, glob = glue("*prior_sim.{site}*.rds"))
  obs_raw <- load_observations(site)
  obs <- spsummary(obs_raw, type = "observed", site = site)
  post_dat <- readRDS(spec_history_file)$result
  prior_dat <- readRDS(prior_file)
  post <- spsummary(post_dat, type = "posterior", site = site)
  prior <- spsummary(prior_dat, type = "prior", site = site)
  bind_rows(post, prior, obs)
}

summary_spec_file <- "proc_results/site_spectra_summaries.rds"

if (file.exists(summary_spec_file)) {
  dat <- readRDS(summary_spec_file)
} else {
  dat <- map_dfr(
    site_list,
    with_prog(
      possibly(get_site_data, NULL),
      nsite,
      list(format = ":current/:total (:eta)")
    )
  ) %>%
    mutate(site = str_extract(site, ".*?(?=_)"))
  saveRDS(dat, summary_spec_file)
}

dat_structure <- read_csv("other_site_data/site_structure.csv") %>%
  rename(site = site_label)

dat_wide <- dat %>%
  select(waves, type, site, mu) %>%
  spread(type, mu)

dat_errors <- dat_wide %>%
  mutate(
    error = posterior - observed,
    rel_error = error / observed,
    is_vis = waves < 750
  ) %>%
  group_by(site, is_vis) %>%
  summarize(
    bias = mean(error, na.rm = TRUE),
    rel_bias = mean(rel_error, na.rm = TRUE),
    rmse = mean(sqrt(error^2), na.rm = TRUE),
    rel_rmse = mean(sqrt(rel_error^2), na.rm = TRUE)
  )

dat_plt <- dat_errors %>%
  left_join(dat_structure, by = "site")

example_sites <- c(
  "NC21", "IDS40", "OF02", "BI04"
)
selected_sites <- readLines("other_site_data/selected_sites") %>%
  grep("SF03", ., value = TRUE, invert = TRUE)
plot_dat <- dat_plt %>%
  mutate(
    is_vis = recode_factor(as.character(is_vis), `TRUE` = "VIS", `FALSE` = "NIR"),
    tot_dens = tot_dens * 1000
  ) %>%
  select(site, is_vis, rel_bias, mean_dbh, tot_dens, frac_evergreen) %>%
  gather(variable, value, mean_dbh, tot_dens) %>%
  mutate(
    variable = recode(
      variable,
      mean_dbh = "Mean DBH (cm)",
      tot_dens = "Stand density (trees/ha)"
    )
  )
plt <- ggplot(plot_dat) +
  aes(x = value, y = rel_bias, color = frac_evergreen) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_label_repel(
    aes(label = site),
    data = filter(plot_dat, site %in% c(example_sites)),
    fontface = "bold",
    min.segment.length = 0
  ) +
  facet_grid(is_vis ~ variable, scales = "free") +
  scale_color_gradient(low = "green4", high = "red") +
  labs(y = "Relative bias", color = "Evergreen fraction") +
  theme(
    strip.placement = "outside",
    axis.title.x = element_blank()
  ) +
  panel_border()
save_plot("figures/bias.pdf", plt, base_width = 7, base_height = 7)

site_spec_ci <- function(sites) {
  dat_sub <- dat %>%
    filter(site %in% sites, type != "observed") %>%
    mutate(type = factor(type, c("prior", "posterior")))
  obs_raw <- load_observations(sites)
  obs <- as.data.frame(obs_raw) %>%
    cbind(waves = PEcAnRTM::wavelengths(obs_raw), .) %>%
    set_tidy_names() %>%
    as_tibble() %>%
    gather(site2, value, -waves) %>%
    separate(site2, c("site", "site_rep"))
  ggplot(dat_sub) +
    aes(x = waves) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = 0.5) +
    geom_line(aes(y = value, group = site_rep), data = obs, color = "black") +
    facet_wrap(~site) +
    scale_fill_manual(values = c("prior" = "grey70", "posterior" = "red")) +
    labs(x = "Wavelength (nm)", y = "Reflectance (0-1)", fill = "") +
    panel_border()
}

spec_example <- site_spec_ci(example_sites)
spec_selected <- site_spec_ci(selected_sites)

both_example <- plot_grid(
  plt + guides(color = FALSE),
  get_legend(plt),
  spec_example + guides(fill = FALSE),
  get_legend(spec_example),
  nrow = 2,
  rel_widths = c(1, 0.3)
)
save_plot("figures/spec_validation.pdf", both_example, base_width = 8, base_height = 10)

lai_raw <- readRDS(file.path(data_dir, "lai_values.rds"))
lai_dat <- lai_raw %>%
  map(rowSums, na.rm = TRUE) %>%
  imap_dfr(~tibble(site = .y, LAI = .x)) %>%
  mutate(site = str_extract(site, "^.*?(?=_)"))
lai_obs_all <- read_csv("other_site_data/NASA_FFT_LAI_FPAR_Data.csv")
lai_sub <- lai_obs_all %>%
  group_by(site = Site_Plot) %>%
  summarize(
    obs_LAI = mean(LAI_CLX_10_60_m2_m2, na.rm = TRUE),
    obs_LAI_SD = sd(LAI_CLX_10_60_m2_m2, na.rm = TRUE),
    obs_LAI_lo = obs_LAI - obs_LAI_SD,
    obs_LAI_hi = obs_LAI + obs_LAI_SD
  )

lai_joined <- lai_dat %>%
  group_by(site) %>%
  summarize(
    LAI_mean = mean(LAI, na.rm = TRUE),
    LAI_SD = sd(LAI, na.rm = TRUE),
    LAI_lo = LAI_mean - LAI_SD,
    LAI_hi = LAI_mean + LAI_SD
  ) %>%
  left_join(lai_sub, by = "site") %>%
  mutate(
    error = LAI_mean - obs_LAI,
    abs_error = abs(error),
    site = fct_reorder(site, abs_error)
  )

lai_plot <- lai_joined %>%
  select(-abs_error, -error) %>%
  gather(variable, value, -site) %>%
  mutate(
    type = if_else(grepl("obs", variable), "observed", "model"),
    variable = recode(variable, obs_LAI = "LAI_mean"),
    variable = str_remove(variable, "obs_")
  ) %>%
  spread(variable, value)

lai_bysite <- ggplot(lai_plot) +
  aes(x = site, y = LAI_mean, ymin = LAI_lo, ymax = LAI_hi, color = type) +
  geom_errorbar(position = "dodge") +
  background_grid() +
  theme(axis.text.x = element_text(angle = 90))
save_plot("figures/lai_errorbars.pdf", lai_bysite, base_width = 18, base_height = 10)

lai_scatter <- lai_plot %>%
  gather(variable, value, -site, -type) %>%
  unite(var_type, variable, type) %>%
  spread(var_type, value) %>%
  left_join(dat_structure) %>%
  ggplot() +
  aes(
    x = LAI_mean_model, xmin = LAI_lo_model, xmax = LAI_hi_model,
    y = LAI_mean_observed, ymin = LAI_lo_observed, ymax = LAI_hi_observed,
    color = frac_evergreen
  ) +
  geom_point(size = 3) +
  geom_errorbarh() +
  geom_errorbar() +
  geom_abline(linetype = "dashed") +
  labs(x = "EDR predicted LAI", y = "Observed LAI", color = "Frac. evergreen") +
  #geom_label_repel(aes(label = site), size = 5) +
  scale_color_gradient(low = "green4", high = "red")
save_plot("figures/lai_sites.pdf", lai_scatter, base_width = 10, base_height = 10)

spec_lai_dat <- plot_dat %>%
  left_join(lai_joined) %>%
  filter(grepl("Mean DBH", variable))

slmod <- spec_lai_dat %>%
  ungroup() %>%
  group_by(is_vis) %>%
  nest() %>%
  mutate(
    fit = map(data, lm, formula = rel_bias ~ error),
    coefs = map(fit, tidy)
  ) %>%
  unnest(coefs)

#slmod_plot <- slmod %>%
  #group_by(is_vis) %>%
  #mutate(maxp = max(p.value)) %>%
  #ungroup() %>%
  #filter(maxp < 0.05) %>%
  #mutate(term = recode())


spec_lai_plot <- ggplot(spec_lai_dat) +
  aes(x = error, y = rel_bias, color = frac_evergreen) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(is_vis ~ ., scales = "free") +
  labs(
    x = "LAI error (model - observed)",
    y = "Relative spectral bias",
    color = "Frac. evergreen"
  ) +
  scale_color_gradient(low = "green4", high = "red") +
  panel_border()
save_plot("figures/spec_lai.pdf", spec_lai_plot, base_width = 7, base_height = 7)
