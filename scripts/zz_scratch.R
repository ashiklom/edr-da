old_plan_targets <- list(
  sail_predictions = tidy_sail_predictions(site_details, site_lai_total,
                                           tidy_posteriors),
  spec_error_all_presentation = ggsave(
    file_out(!!path(figdir, "spec-error-all-presentation.png")),
    spec_error_all_f(observed_predicted, sail_predictions, ncol = 8) +
      theme(axis.text.x = element_text(size = rel(0.75)),
            legend.position = c(1, 0),
            legend.justification = c(1, 0),
            legend.direction = "horizontal"),
    width = 10, height = 8, dpi = 300
  ),
  both_ndvi = calc_ndvi_bysite(observed_spectra, predicted_spectra,
                               site_structure_data),
  ndvi_dbh_png = ggsave(
    file_out(!!path(figdir, "ndvi-dbh.png")),
    ndvi_dbh_plot(both_ndvi),
    width = 6.4, height = 5.2, dpi = 300
  )
)

##################################################

drake::loadd(inversion_site_list)
params_matrix <- drake::readd(posterior_matrix)
site <- inversion_site_list[3]
dedup <- TRUE
progress <- FALSE

predicted_spectra <- results
predited_spectra <- ungroup(predicted_spectra)

z <- tidy_site_spec(site)

##################################################

drake::loadd(observed_predicted)

get_ar <- function(x) {
  fit <- ar(x, order.max = 1)
  rho <- fit$ar
  if (!length(rho) > 1) {
    rho
  } else {
    NA_real_
  }
}

bias_data <- observed_predicted %>%
  group_by(site) %>%
  arrange(wavelength, .by_group = TRUE) %>%
  summarize(rho = get_ar(bias)) %>%
  ungroup()

bias_data %>%
  pull(rho) %>%
  mean()

ggplot(bias_data) +
  aes(x = rho) +
  geom_histogram()

lapply(bias_data$fit, "[[", "ar")
bias_data$fit[[1]]$ar

##################################################

## Old process aviris code

# Observation has:
# - Reflectance
# - Site ID
# - Solar zenith angle (and azimuth?)
# - ...or just store all as EDR arguments as a closure?

## aviris_sf <- aviris_raw %>%
##   select(aviris_id, img_date, lon = ptX, lat = ptY) %>%
##   st_as_sf(coords = c("lon", "lat"), crs = 4326)



aviris_spectra <- aviris_raw %>%
  select(starts_with("band_")) %>%
  as.matrix() %>%
  t() %>%
  spectra(wavelengths = aviris_waves)

colnames(aviris_spectra) <- aviris_raw$iPLOT

# Remove water bands
aviris_wl <- wavelengths(aviris_spectra)
aviris_filter <- (aviris_wl > 1340 & aviris_wl < 1425) |
  (aviris_wl > 1800 & aviris_wl < 1950) |
  (aviris_wl < 400)
aviris_spectra[aviris_filter, ] <- NA

saveRDS(aviris_spectra, "aviris/aviris.rds")

## Store AVIRIS spectra in HDF5 file
#hdf_file <- h5file(here("aviris", "aviris.h5"), mode = "a")

#hdf_file[["spectral_data"]] <- aviris_spectra
#hdf_file[["wavelengths"]] <- aviris_wl
#hdf_file[["band_numbers"]] <- seq_along(aviris_wl)

#aviris_meta <- aviris_raw %>%
  #select(-starts_with("band_"))

#for (column in colnames(aviris_meta)) {
  #hdf_file[[column]] <- aviris_raw[[column]]
#}

#h5close(hdf_file)

##################################################
# Evaluate PROSPECT covariances
drake::loadd(posterior_matrix)

posterior_cov <- cov(posterior_matrix)

pcov_df <- corrr::as_cordf(posterior_cov)

pcov_tidy <- pcov_df %>%
  corrr::shave() %>%
  corrr::stretch(na.rm = TRUE) %>%
  dplyr::mutate(across(c(x, y), ~gsub("temperate\\.", "", .x)))

pcov_tidy %>%
  dplyr::filter(across(c(x, y), ~grepl("prospect", .x))) %>%
  arrange(desc(abs(r))) %>%
  print(n = Inf)

##################################################

drake::loadd(site_lai_total)
drake::loadd(lai_observed)
drake::loadd(observed_predicted)

lai_both <- site_lai_total %>%
  left_join(lai_observed, "site") %>%
  left_join(site_structure_data, c("site" = "site_name")) %>%
  mutate(
    LAI_diff = lai_mean - obs_LAI,
    eLAI_diff = elai_mean - obs_LAI
  ) %>%
  select(site, LAI_diff, lai_mean, obs_LAI, mean_dbh, tot_dens, everything()) %>%
  arrange(desc(LAI_diff))

lai_both %>%
  select(site, lai_mean, elai_mean, obs_LAI, mean_dbh, tot_dens) %>%
  tidyr::pivot_longer(c(lai_mean, elai_mean)) %>%
  ggplot() +
  aes(x = value, y = obs_LAI) +
  geom_abline(linetype = "dashed") +
  geom_smooth() +
  geom_point() +
  facet_grid(vars(name))

lai_both %>%
  select(site, LAI_diff, eLAI_diff, obs_LAI, mean_dbh, tot_dens) %>%
  tidyr::pivot_longer(c(LAI_diff, eLAI_diff)) %>%
  ggplot() +
  aes(x = mean_dbh, y = value) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(cols = vars(name))


xmx <- 40
ymx <- 60
plot(x = 0, y = 0, type = "n", xlim = c(0, xmx), ylim = c(0, ymx), xlab = "DBH", ylab = "BL")
for (b2 in seq(1.5, 2.2, 0.05)) {
  curve(size2bl(x, 0.1, b2), 0, xmx, add = TRUE, col = alpha("blue4"))
}
## plot(x = 0, y = 0, type = "n", xlim = c(0, 60), ylim = c(0, 100), xlab = "DBH", ylab = "BL")
for (b1 in seq(0.05, 0.2, 0.01)) {
  curve(size2bl(x, b1, 1.8), 0, ymx, add = TRUE, col = alpha("red4"))
}

p1 <- ggplot(lai_both) +
  aes(x = elai_mean, y = obs_LAI) +
  geom_abline(linetype = "dashed") +
  geom_smooth(method = "lm") +
  geom_point()

p2 <- ggplot(lai_both) +
  aes(x = lai_mean, y = obs_LAI) +
  geom_abline(linetype = "dashed") +
  geom_smooth(method = "lm") +
  geom_point()

p1 + p2

ggplot(lai_both) +
  aes(x = mean_dbh, y = eLAI_diff) +
  ## aes(x = tot_dens, y = LAI_diff) +
  ## aes(x = frac_evergreen_wtd, y = LAI_diff) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed")

band_diffs <- observed_predicted %>%
  dplyr::mutate(band = if_else(wavelength < 750, "VIS", "NIR")) %>%
  dplyr::group_by(band, site) %>%
  dplyr::summarize(mean_diff = mean(bias, na.rm = TRUE)) %>%
  dplyr::ungroup()

resid2 <- band_diffs %>%
  left_join(lai_both, "site")

ggplot(resid2) +
  aes(x = LAI_diff, y = mean_diff) +
  geom_smooth(method = "lm") +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(vars(band), scales = "free_y")


##################################################
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

both_wide %>%
  dplyr::filter(!grepl("sitesoil", variable)) %>%
  dplyr::arrange(variable, rel_reduction) %>%
  print(n = Inf)

##################################################

sites <- readLines("other_site_data/site_list")

obs <- load_observations()

site_structure_data %>%
  dplyr::arrange(dplyr::desc(tot_dens))

site_structure_data %>%
  dplyr::filter(mean_dbh > 35)

library(dplyr)

count(tidy_posteriors, variable)

resid <- tidy_posteriors %>%
  filter(variable %in% c("residual_intercept", "residual_slope"),
         type == "posterior") %>%
  group_by(variable) %>%
  summarize(Mean = mean(value))
rm <- resid %>% filter(variable == "residual_slope") %>% pull(Mean)
rb <- resid %>% filter(variable == "residual_intercept") %>% pull(Mean)

obs_sub <- observed_predicted %>%
  filter(iobs == iobs[1], site == unique(site)[8]) %>%
  arrange(wavelength)
wl <- obs_sub[["wavelength"]]
obs <- obs_sub[["observed"]]
pred <- obs_sub[["albedo_mean"]]

rsd <- rb + rm * pred

plot(wl, obs, type = "l")
lines(wl, pred, col = "red")

ar_fit <- ar(obs, order.max = 1)

times <- seq_along(wl)
rho <- ar_fit$ar
H <- abs(outer(times, times, "-"))
covar <- diag(rsd) %*% (rho^H) %*% diag(rsd)

image(wl, wl, covar)

logl_diag <- mvtnorm::dmvnorm(pred, obs, diag(rsd^2), log = TRUE)
logl_ar1 <- mvtnorm::dmvnorm(pred, obs, covar, log = TRUE)

exp(logl_ar1 - logl_diag)

observed_predicted

bias_bysite <- observed_predicted %>%
  dplyr::select(iobs, site, wavelength, bias) %>%
  dplyr::group_by(iobs, site) %>%
  dplyr::arrange(wavelength, .by_group = TRUE)

results <- bias_bysite %>%
  dplyr::summarize(
    autocorr = list(acf(bias, lag.max = 50)),
    my_neff = neff(bias),
    coda_neff = coda::effectiveSize(bias)
  )

summary(results$coda_neff)
summary(results$my_neff)
hist(results$coda_neff)
hist(results$my_neff)
mean(results$coda_neff)

plot(results$autocorr[[3]])

observed_predicted %>%
  dplyr::count(site, iobs)

##################################################

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

##################################################
curve(dlnorm(x, 1, 0.5), 0, 12, n = 1000)
curve(dlnorm(x, 0, 1), 0, 5, n = 1000)
## curve(plnorm(x, 1, 1), 0, 12, n = 1000)

##################################################

drake::loadd(posterior_matrix)
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

lc2 %>%
  dplyr::filter(dplyr::across(c(xpft, ypft), ~!grepl("sitesoil", .x))) %>%
  dplyr::arrange(dplyr::desc(r)) %>%
  print(n = 20)

lc2 %>%
  ## dplyr::filter(xvar == "clumping", xpft == "North_Mid_Hardwood") %>%
  dplyr::filter(xvar == "orient", xpft == "Late_Hardwood") %>%
  dplyr::arrange(dplyr::desc(abs(r)))

lc3 <- lc2 %>%
  dplyr::filter(xpft == ypft) %>%
  dplyr::mutate(pft = xpft)

plot_dat <- lc3 %>%
  dplyr::filter(!is.na(pft), !grepl("sitesoil", pft),
                !is.na(xvar), !is.na(yvar)) %>%
  dplyr::mutate(dplyr::across(c(pft, xvar, yvar), forcats::fct_inorder))

ggplot(plot_dat) +
  aes(x = 1, y = r, fill = pft) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(position = "dodge") +
  facet_grid(vars(xvar), vars(yvar), drop = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0))

ggsave("text/figures/posterior-correlations.png", width = 7, height = 7, dpi = 300, units = "in")

wide_cor %>%
  corrr::rearrange() %>%
  dplyr::mutate(term = gsub("temperate\\.", "", term)) %>%
  dplyr::slice(1:10) %>%
  dplyr::select(1:5)
## corrr::network_plot(wide_cor)

## wide_cor %>%
##   corrr::focus(matches("Early_Hardwood"), mirror = TRUE) %>%
##   corrr::rplot()

long_cor %>%
  dplyr::filter(dplyr::across(c(x, y), ~grepl("Early_Hardwood", .x))) %>%
  corrr::rplot()

curve(dlnorm(x, 1, 0.5), 10, 20, add = TRUE, col = "red")


plnorm(c(2, 5), 1, 0.5)
qlnorm(c(0.025, 0.5, 0.975), 1, 0.5)

curve(dlnorm(x, 1, 0.5), 0.1, 10)
