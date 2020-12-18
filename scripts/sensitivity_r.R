pkgload::load_all(".")

library(ggplot2)
library(magrittr, include.only = "%>%")

# Single cohort
edr_defaults <- list(
  N = 1.4,
  Cab = 40,
  Car = 10,
  Cw = 0.01,
  Cm = 0.01,
  lai = 1,
  cai = 1,
  clumping_factor = 1,
  orient_factor = 0,
  direct_sky_frac = 0.9,
  pft = 1,
  czen = 1,
  wai = 0,
  soil_moisture = 0.5
)

sail_defaults <- list(
  N = edr_defaults[["N"]],
  Cab = edr_defaults[["Cab"]],
  Car = edr_defaults[["Car"]],
  Cw = edr_defaults[["Cw"]],
  Cm = edr_defaults[["Cm"]],
  Cbrown = 0,
  LAI = edr_defaults[["lai"]] * edr_defaults[["clumping_factor"]],
  soil_moisture = edr_defaults[["soil_moisture"]],
  hot_spot = 0,
  solar_zenith = acos(edr_defaults[["czen"]]),
  LIDFa = edr_defaults[["orient_factor"]],
  LIDFb = 0
)

do_sens <- function(value, variable, fun, .dots = list()) {
  stopifnot(is.list(.dots))
  varlist <- list()
  varlist[[variable]] <- value
  arglist <- modifyList(.dots, varlist)
  do.call(fun, arglist)
}

tidy_albedo <- function(result_list, values) {
  stopifnot(length(result_list) == length(values))
  values_df <- tibble::tibble(
    variable = paste0("V", seq_along(values)),
    var_value = values
  )
  names(result_list) <- values_df[["variable"]]
  albedo_dfw <- purrr::map_dfc(result_list, "albedo")
  albedo_dfw[["wavelength"]] <- seq(400, 2500)
  albedo_long <- tidyr::gather(albedo_dfw, variable, value, -wavelength)
  dplyr::left_join(albedo_long, values_df, by = "variable")
}

## result_list <- lai_sens_sail
## values <- lai
tidy_sail <- function(result_list, values) {
  stopifnot(length(result_list) == length(values))
  results_dfl <- purrr::map(result_list, tibble::as_tibble) %>%
    purrr::map(function(x) {x$wavelength <- seq(400, 2500); x})
  values_df <- tibble::tibble(
    variable = paste0("V", seq_along(values)),
    var_value = values,
    saildata = results_dfl
  )
  tidyr::unnest(values_df, saildata)
}

sens_plot <- function(tidy_sens) {
  ggplot(tidy_sens) +
    aes(x = wavelength, y = value, color = var_value,
        group = factor(var_value)) +
    geom_line() +
    labs(x = "Wavelength (nm)", y = "Albedo")
}

lai <- c(
  seq(0.2, 1, 0.2),
  seq(1, 2, 0.4),
  seq(2, 5, 0.5)
)
lai_sens <- purrr::map(
  lai, do_sens,
  fun = edr_r,
  variable = "lai",
  .dots = modifyList(edr_defaults, list(direct_sky_frac = 0.8))
) %>%
  tidy_albedo(lai)
## sens_plot(lai_sens) + ggtitle("LAI sensitivity") + scale_color_viridis_c()
lai_sens_sail <- purrr::map(lai, do_sens, fun = rrtm::pro4sail_5,
                            variable = "LAI", .dots = sail_defaults) %>%
  tidy_sail(lai)
tidy_both <- dplyr::left_join(lai_sens, lai_sens_sail)
tidy_both_long <- tidy_both %>%
  tidyr::pivot_longer(c(value, bhr:bdr))
## ggplot(tidy_both_long) +
##   aes(x = wavelength, y = value, group = variable, color = var_value) +
##   geom_line() +
##   facet_grid(cols = vars(name)) +
##   scale_color_viridis_c()
ggplot(tidy_both_long) +
  aes(x = wavelength, y = value, color = name) +
  geom_line() +
  facet_wrap(vars(var_value))

ggplot(tidy_both) +
  aes(x = wavelength, group = variable, color = value) +
  geom_line(aes(color = ))
# Compare sensitivity of EDR and SAIL

# How does sensitivity to soil moisture change with different LAI values?
do_soil_sens <- function(lai, soil_vals = seq(0, 1, 0.1)) {
  purrr::map(
    soil_vals,
    edr_sens,
    variable = "soil_moisture",
    .dots = modifyList(edr_defaults, list(lai = lai))
  ) %>%
    tidy_albedo(soil_vals) %>%
    dplyr::mutate(lai = lai)
}

soil_lai <- purrr::map_dfr(c(0.5, seq(1, 5)), do_soil_sens)

ggplot(soil_lai) +
  aes(x = wavelength, y = value, color = var_value, group = var_value) +
  geom_line() +
  facet_wrap(~lai) +
  scale_color_viridis_c() +
  labs(x = "Wavelength", y = "Albedo", color = "Soil moisture frac") +
  ggtitle("Albedo sensitivity to soil moisture for various LAI")

ggsave("figures/lai_soil_sensitivity.png", width = 8, height = 4)

# Do different LAI combinations of identical PFTs produce identical
# results?
edr_defaults_sub <- edr_defaults[!names(edr_defaults) %in% c("lai", "pft", "cai", "wai")]
laicomb <- tibble::tibble(
  !!!edr_defaults_sub,
  lai = list(
    # Combinations of 0.5
    c(0.25, 0.25), 0.5,
    # Combinations of 1
    c(0.5, 0.5), 1,
    # Combinations of 2
    c(1, 1), c(1.5, 0.5), c(0.5, 1.5), 2,
    # Combinations of 3
    c(1.5, 1.5), c(2, 1), c(1, 2), 3
  ),
  pft = list(
    c(1, 1), 1, c(1, 1), 1, c(1, 1), c(1, 1), c(1, 1), 1,
    c(1, 1), c(1, 1), c(1, 1), 1
  ),
  cai = pft,
  wai = list(
    c(0, 0), 0, c(0, 0), 0, c(0, 0), c(0, 0), c(0, 0), 0,
    c(0, 0), c(0, 0), c(0, 0), 0
  ),
  ) %>%
  dplyr::mutate(
    lailab = purrr::map_chr(lai, paste, collapse = ", "),
    total_lai = purrr::map_dbl(lai, sum),
    raw_out = purrr::pmap(., edr_r),
    albedo = purrr::map(raw_out, "albedo") %>%
      purrr::map(~dplyr::mutate(
        tibble::as_tibble(.),
        wavelength = seq(400, 2500)
      ))
  )
laicomb_long <- tidyr::unnest(laicomb, albedo)

ggplot(laicomb_long) +
  aes(x = wavelength, y = value, color = lailab) +
  geom_line() +
  scale_color_viridis_d() +
  facet_wrap(~total_lai)

# Answer: Yes they do.
