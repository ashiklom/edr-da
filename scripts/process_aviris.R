library(conflicted)
conflict_prefer("filter", "dplyr")

library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(here)
library(sf)

pkgload::load_all()

aviris_file <- here("aviris", "NASA_FFT",
                    "allPlotsAggregated_NIT_A_SpecJoined.csv")

aviris_waves <- read_csv("aviris/NASA_FFT/aviris_c_wavelength.csv", col_types = "d")$wavelength

aviris_raw <- read_csv(aviris_file) %>%
  mutate(
    aviris_id = row_number(),
    flightline = str_extract(IMG, ".*(?=rdn_)"),
    img_date = str_extract(flightline, "(?<=f)[[:digit:]]+") %>%
      lubridate::ymd(),
    img_doy = lubridate::yday(img_date),
    # HACK: Assume observations were collected at 10:30am
    czen = cos_solar_zenith_angle(img_doy, ptY, ptX, 10.5)
  )

flines <- read_csv("aviris/avirisc-flightlines.csv")

dat <- flines %>%
  semi_join(aviris_raw, c("Name" = "flightline"))

dat %>%
  pull(Comments)

# TODO: Grab sun angle information from flight lines table
# TODO: Parse cloud cover info?

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
