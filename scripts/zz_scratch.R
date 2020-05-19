samples <- BayesianTools::getSample(
  samples_bt,
  ## thin = "auto",
  ## start = 10000,
  coda = TRUE
)

# Gelman diagnostic
block_size <- 10000
block_step <- 5000
nsamples <- nrow(samples[[1]])
nblocks <- (samples)
gd <- list()
gw <- list()

iblock <- 0
done <- FALSE
while (!done) {
  iblock <- iblock + 1
  message("Block ", iblock)
  block_start <- 1 + (iblock - 1) * block_step
  block_end <- block_start + (block_size - 1)
  if (block_end > nsamples) {
    block_end <- nsamples
    done <- TRUE
    message("Last one!")
  }
  s_sub <- window(samples, block_start, block_end)
  tag <- paste(block_start, block_end, sep = ":")
  gd[[tag]] <- coda::gelman.diag(s_sub, autoburnin = FALSE, multivariate = FALSE)$psrf[,1]
  gw[[tag]] <- coda::geweke.diag(s_sub)
}

coda::gelman.diag(samples)

s_sub <- window(samples, start = 5000)
coda::gelman.diag(s_sub, autoburnin = FALSE)
plot(s_sub[,1])

plot(samples[,1])

plot(purrr::map_dbl(gd, "par 1"), type = 'l')

pnorm(gw[[17]][[1]]$z) > 0.05
plot(s_sub[,"par 13"])


zzz <- coda::as.mcmc(rnorm(5000))
coda::geweke.diag(zzz)
plot(zzz)


##################################################
zzz <- rnorm(10)

# Convert to CDF -- 0-1
zzz_p <- pnorm(zzz)

# Transform to target distribution
devtools::load_all()
mvtnorm::qmvnorm(0.1, mean = prospect_means[1,], sigma = prospect_covars[,,1])

##################################################
sites <- readLines(here::here("other_site_data", "site_list"))
site_structure <- readr::read_csv("other_site_data/site_structure.csv") %>%
  dplyr::mutate(site_tag = paste0("sitesoil_", dplyr::row_number()))


site_posterior_summary <- site_posterior %>%
  dplyr::group_by(variable) %>%
  dplyr::summarize(
    Mean = mean(value),
    SD = sd(value),
    lo = quantile(value, 0.025),
    hi = quantile(value, 0.975)
  ) %>%
  dplyr::ungroup()

site_data <- site_structure %>%
  dplyr::inner_join(site_posterior_summary, c("site_tag" = "variable")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

ggplot(site_data) +
  aes(x = tot_dens, y = Mean, ymin = lo, ymax = hi,
      color = frac_evergreen_wtd) +
  geom_pointrange() +
  scale_color_viridis_c() +
  labs(x = "Stand density", y = "Relative soil moisture (0 - 1)",
       color = "Weighted evergreen fraction") +
  guides(color = guide_colorbar(title.position = "top", direction = "horizontal")) +
  theme_bw() +
  theme(legend.position = c(0.95, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_blank())

basemap <- rnaturalearth::ne_states() %>%
  sf::st_as_sf()

basemap_sub <- basemap %>%
  sf::st_crop(sf::st_buffer(sf::st_as_sfc(sf::st_bbox(site_data)), 2))
ggplot(site_data) +
  geom_sf(data = basemap_sub) +
  geom_sf(aes(color = Mean), size = 4, pch = "x") +
  scale_color_viridis_c() +
  theme_bw()

clrs <- c("prior" = "gray70", "posterior" = "black")
ggplot(dplyr::filter(tidy_posteriors, grepl("sitesoil", variable))) +
  aes(x = variable, y = value,
      fill = type, color = type) +
  geom_violin(position = position_identity()) +
  scale_color_manual(values = clrs, aesthetics = c("color", "fill")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

##################################################
ggplot(site_lai_summary) +
  aes(x = elai_mean, y = hite, color = pft) +
  geom_segment(aes(yend = hite, xend = 0)) +
  geom_point() +
  facet_wrap(vars(site), scales = "free")

plt_sites <- c("OF04", "SF03", "OF01", "GR08")
plt_sites <- c("IDS35", "IDS34", "BH05")
site_lai_summary %>%
  ## filter(site %in% plt_sites) %>%
  ggplot() +
  aes(x = cum_elai, y = hite, color = pft) +
  geom_segment(aes(yend = hite, xend = 0)) +
  geom_point() +
  facet_wrap(vars(site), scales = "fixed")

site_lai_summary

##################################################
bad_sites <- site_lai_total %>%
  filter(elai_mean > 10) %>%
  pull(site)

inner_join(site_lai_total, lai_observed, "site") %>%
  ## filter(!site %in% bad_sites) %>%
  ggplot() +
  aes(x = lai_mean, y = obs_LAI,
      xmin = lai_lo, xmax = lai_hi) +
  geom_pointrange() +
  geom_abline(linetype = "dashed")
  ## coord_equal()

##################################################
nc18_spectra <- predict_site_spectra(
  posterior_matrix, "NC18",
  nsamp = 1000, dedup = TRUE, progress = TRUE
)

nc18_obs <- load_observations("NC18") %>%
  `colnames<-`(., as.character(seq_len(NCOL(.)))) %>%
  as.data.frame(., row.names = PEcAnRTM::wavelengths(.)) %>%
  tibble::as_tibble(rownames = "wavelength") %>%
  dplyr::mutate(wavelength = as.numeric(wavelength)) %>%
  tidyr::pivot_longer(-wavelength, names_to = "iobs", values_to = "observed")

ggplot(nc18_spectra) +
  aes(x = wavelength, y = albedo_mean,
      ymin = albedo_q025, ymax = albedo_q975) +
  ## geom_ribbon(aes(ymin = pmax(albedo_r_q025, 0), ymax = albedo_r_q975),
  ##             fill = "orange", alpha = 0.5) +
  geom_ribbon(fill = "skyblue") +
  geom_line(color = "blue4") +
  geom_line(aes(x = wavelength, y = observed, group = iobs),
            inherit.aes = FALSE,
            data = nc18_obs)

##################################################
mubar <- function(orient_factor) {
  phi1 <- 0.5 - orient_factor * (0.633 + 0.33 * orient_factor)
  phi2 <- 0.877 * (1 - 2 * phi1)
  mu_bar <- (1 - phi1 * log(1 + phi2 / phi1) / phi2) / phi2
  mu_bar
}

##################################################
Sys.setenv("RUN_CONFIG" = "hetero-pooled")
loadd("posterior_matrix", cache = dc)
loadd("inversion_site_list", cache = dc)
loadd("site_lai_total", cache = dc)
loadd("site_details", cache = dc)
loadd("lai_observed", cache = dc)

site_lai_total %>%
  arrange(desc(lai_mean))

##################################################
site <- "BI06"
site_df <- site_details %>% filter(site == !!site)
obs <- lai_observed %>% filter(site == !!site)
b1Bl <- purrr::map_dbl(allom_mu, "b1Bl") %>% exp()
names(b1Bl) <- gsub("temperate\\.", "", names(b1Bl))
b2Bl <- purrr::map_dbl(allom_mu, "b2Bl")
names(b2Bl) <- gsub("temperate\\.", "", names(b2Bl))
pm <- prospect_means
rownames(pm) <- gsub("temperate\\.", "", rownames(pm))
SLA <- pm[,"SLA"]

albl <- list()
lail <- list()
for (bbnp in c(0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.5, 1, 2, 3, 4)) {
  b1Bl["Northern_Pine"] <- bbnp
  site_lai <- site_df %>%
    mutate(bleaf = size2bl(dbh, b1Bl[pft], b2Bl[pft]),
           lai = nplant * SLA[pft] * bleaf,
           hite = dbh2h(dbh, ipft)) %>%
    arrange(hite)
  tag <- as.character(bbnp)
  lail[[tag]] <- site_lai
  albl[[tag]] <- edr_r(
    site_lai$pft,
    site_lai$lai,
    rep(0, nrow(site_lai)),
    rep(1, nrow(site_lai)),
    N = pm[,"prospect_N"],
    Cab = pm[, "prospect_Cab"],
    Car = pm[, "prospect_Car"],
    Cw = pm[, "prospect_Cw"],
    Cm = pm[, "prospect_Cm"],
    orient_factor = rep(0.5, 5),
    clumping_factor = rep(1, 5),
    soil_moisture = 0.8,
    direct_sky_frac = 0.9,
    czen = 1
  )[["albedo"]]
}

lai_sums <- purrr::map_dbl(lail, function(x) sum(x$lai))

r <- do.call(cbind, albl)
matplot(400:2500, r, type = "l")
legend(
  "topright",
  as.character(round(lai_sums, 1)),
  lty = seq_along(lai_sums),
  col = seq_along(lai_sums)
)

plot(out$albedo, type = 'l')

##################################################
drake::loadd("site_details", cache = dc)
drake::loadd("lai_observed", cache = dc)
drake::loadd("site_lai_total", cache = dc)
drake::loadd("tidy_posteriors", cache = dc)

ggplot(sail_tidy) +
  aes(x = wavelength, y = value, linetype = "stream") +
  geom_line() +
  facet_wrap(vars(site), scales = "fixed")

tidy_posteriors
