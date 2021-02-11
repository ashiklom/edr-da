library(conflicted)
conflict_prefer("filter", "dplyr")

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(lubridate)
library(here)
library(sf)

pkgload::load_all()

aviris_file <- here("aviris", "NASA_FFT",
                    "allPlotsAggregated_NIT_A_SpecJoined.csv")

aviris_wl <- read_csv("aviris/NASA_FFT/aviris_c_wavelength.csv", col_types = "d")$wavelength
aviris_filter <- (aviris_wl > 1340 & aviris_wl < 1425) |
  (aviris_wl > 1800 & aviris_wl < 1950) |
  (aviris_wl < 400)
aviris_na_bands <- paste0("band_", seq_along(aviris_wl))[aviris_filter]

aviris_flines <- read_csv("aviris/avirisc-flightlines.csv") %>%
  distinct(Name, .keep_all = TRUE)

aviris_raw <- read_csv(aviris_file) %>%
  mutate(
    aviris_id = row_number(),
    flightline = str_extract(IMG, ".*(?=rdn_)"),
    img_date = str_extract(flightline, "(?<=f)[[:digit:]]+") %>%
      lubridate::ymd(),
    img_doy = lubridate::yday(img_date),
  ) %>%
  left_join(aviris_flines, c("flightline" = "Name")) %>%
  mutate(
    # Assume solar time is at the mean longitude
    across(c(`UTC Hour`, `UTC Minute`, `Solar Elevation`), ~na_if(.x, -999)),
    lon_mean = (Lon1 + Lon2 + Lon3 + Lon4)/4,
    solar_time = (`UTC Hour` + `UTC Minute`/60) + lon_mean * 24/360,
    # If not provided, assume the solar time is 10:30am
    solar_time = if_else(is.na(solar_time), 10.5, solar_time),
    solar_zenith = 90 - `Solar Elevation`,
    czen = if_else(
      is.na(solar_zenith),
      cos_solar_zenith_angle(img_doy, ptY, ptX, solar_time),
      cos(solar_zenith * pi/180)
    )
  )

merra_data <- read_csv(here("other_site_data", "merra-radiation.csv")) %>%
  mutate(
    direct_sky_frac = nir_direct / (nir_direct + nir_diffuse),
    # HACK: If it's too early, try the next hour's DSF.
    direct_sky_frac = if_else(is.finite(direct_sky_frac), direct_sky_frac, dplyr::lead(direct_sky_frac))
  ) %>%
  filter(is.finite(direct_sky_frac))

merra_sub <- merra_data %>%
  select(latitude, longitude, date, time, direct_sky_frac)

aviris_complete <- aviris_raw %>%
  mutate(solar_date_hour = img_date + lubridate::dhours(round(solar_time))) %>%
  left_join(merra_sub, c("ptX" = "longitude", "ptY" = "latitude", "solar_date_hour" = "time"))

aviris_final <- aviris_complete %>%
  select(iPLOT, aviris_id, flightline, img_date, img_doy,
         czen, direct_sky_frac, starts_with("band_")) %>%
  mutate(across(!!aviris_na_bands, ~NA_real_))

write_csv(aviris_final, here("other_site_data", "aviris-processed.csv"))
