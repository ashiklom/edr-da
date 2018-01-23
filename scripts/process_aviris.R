library(tidyverse)
library(hdf5r)
library(here)
library(PEcAnRTM)

data(raw.sensor.data, package = "PEcAnRTM")

aviris_file <- here("aviris", "NASA_FFT",
                    "allPlotsAggregated_NIT_A_SpecJoined.csv")

aviris_raw <- read_csv(aviris_file) %>%
  mutate(
    aviris_id = row_number()
  )

aviris_spectra <- aviris_raw %>%
  select(starts_with("band_")) %>%
  as.matrix() %>%
  t() %>%
  spectra(wavelengths = fwhm.aviris.classic$avg)

# Remove water bands
aviris_wl <- wavelengths(aviris_spectra)
aviris_filter <- (aviris_wl > 1340 & aviris_wl < 1425) |
  (aviris_wl > 1800 & aviris_wl < 1950) |
  (aviris_wl < 400)
aviris_spectra[aviris_filter, ] <- NA

# Store AVIRIS spectra in HDF5 file
hdf_file <- h5file(here("aviris", "aviris.h5"), mode = "a")

hdf_file[["spectral_data"]] <- aviris_spectra
hdf_file[["wavelengths"]] <- aviris_wl
hdf_file[["band_numbers"]] <- seq_along(aviris_wl)

aviris_meta <- aviris_raw %>%
  select(-starts_with("band_"))

for (column in colnames(aviris_meta)) {
  hdf_file[[column]] <- aviris_raw[[column]]
}

h5close(hdf_file)
