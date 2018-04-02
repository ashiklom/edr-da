library(tidyverse)
library(here)
library(PEcAnRTM)

aviris_file <- here("aviris", "NASA_FFT",
                    "allPlotsAggregated_NIT_A_SpecJoined.csv")

aviris_waves <- read_csv("aviris/NASA_FFT/aviris_c_wavelength.csv", col_types = "d")$wavelength

aviris_raw <- read_csv(aviris_file) %>%
  mutate(
    aviris_id = row_number()
  )

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
