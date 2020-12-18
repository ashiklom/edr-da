devtools::load_all(".")

library(ggplot2)
import::from(magrittr, "%>%")

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

edr_sens <- function(value, variable, .dots = list()) {
  stopifnot(is.list(.dots))
  varlist <- list()
  varlist[[variable]] <- value
  arglist <- modifyList(.dots, varlist)
  do.call(edr_r, arglist)
}

tidy_albedo <- function(result_list, values) {
  stopifnot(length(result_list) == length(values))
  values_df <- tibble::tibble(
    variable = paste0("V", seq_along(values)),
    var_value = values
  )
  albedo_dfw <- purrr::map_dfc(result_list, "albedo")
  albedo_dfw[["wavelength"]] <- seq(400, 2500)
  albedo_long <- tidyr::gather(albedo_dfw, variable, value, -wavelength)
  dplyr::left_join(albedo_long, values_df, by = "variable")
}

sens_plot <- function(tidy_sens) {
  ggplot(tidy_sens) +
    aes(x = wavelength, y = value, color = factor(var_value)) +
    geom_line() +
    labs(x = "Wavelength (nm)", y = "Albedo")
}

lai <- c(
  seq(0.2, 1, 0.2),
  seq(1, 2, 0.4),
  seq(2, 5, 0.5)
)
lai_sens <- purrr::map(lai, edr_sens, variable = "lai", .dots = edr_defaults) %>%
  tidy_albedo(lai)
sens_plot(lai_sens) + ggtitle("LAI sensitivity") + scale_color_viridis_d()
albedo_list <- purrr::map(sensitivity_raw, "albedo")
abedo_mat <- Reduce(cbind, albedo_list)
matplot(albedo_mat, type = "l", main = "LAI sensitivity")

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

ggsave("figures/lai_soil_sensitivity.PDP", width = 8, height = 4)

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
