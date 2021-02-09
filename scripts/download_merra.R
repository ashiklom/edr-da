library(dplyr)
library(readr)
library(tidyr)
library(here)
library(fs)
library(ncdf4)

pkgload::load_all(".")

aviris_file <- here("aviris", "NASA_FFT",
                    "allPlotsAggregated_NIT_A_SpecJoined.csv")
aviris_raw <- read_csv(aviris_file) %>%
  mutate(
    aviris_id = row_number(),
    flightline = str_extract(IMG, ".*(?=rdn_)"),
    img_date = str_extract(flightline, "(?<=f)[[:digit:]]+") %>%
      lubridate::ymd(),
    img_doy = lubridate::yday(img_date),
  )

merra_outdir <- here("other_site_data", "merra")
merra_inputs <- aviris_raw %>%
  distinct(img_date, ptX, ptY) %>%
  rename(date = img_date, longitude = ptX, latitude = ptY)

# Download raw MERRA files
purrr:::pwalk(merra_inputs, get_merra_date, outdir = merra_outdir)

merra_files <- merra_inputs %>%
  mutate(
    ftag = sprintf("lat%.2f-lon%.2f-%s.nc", latitude, longitude, date),
    lfo_name = path(merra_outdir, paste0("merra-lfo-", ftag)),
    most_name = path(merra_outdir, paste0("merra-most-", ftag))
  )

read_lfo <- function(lfo_name, date) {
  nc <- nc_open(lfo_name)
  on.exit(nc_close(nc), add = TRUE)
  tibble::tibble(
    par_time = date + lubridate::dhours(0:23),
    par_direct = ncvar_get(nc, "PARDR"),
    par_diffuse = ncvar_get(nc, "PARDF")
  )
}

read_most <- function(most_name, date) {
  nc <- nc_open(most_name)
  on.exit(nc_close(nc), add = TRUE)
  tibble::tibble(
    nir_time = date + lubridate::dhours(0:23),
    nir_direct = ncvar_get(nc, "NIRDR"),
    nir_diffuse = ncvar_get(nc, "NIRDF")
  )
}

merra_data <- merra_files %>%
  mutate(
    par = purrr::map2(lfo_name, date, read_lfo),
    nir = purrr::map2(most_name, date, read_most)
  ) %>%
  unnest(c(par, nir)) %>%
  mutate(time = par_time) %>%
  select(latitude, longitude, date, time,
         par_direct, par_diffuse, nir_direct, nir_diffuse,
         lfo_name, most_name)

write_csv(merra_data, here("other_site_data", "merra-radiation.csv"))

if (interactive()) {
  # Summary plot
  library(ggplot2)
  merra_summary <- merra_data %>%
    filter(lubridate::hour(time) > 9, lubridate::hour(time) < 15) %>%
    mutate(par_dsf = par_direct / (par_direct + par_diffuse),
           nir_dsf = nir_direct / (nir_direct + nir_diffuse)) %>%
    group_by(longitude, latitude, date) %>%
    summarize(across(matches("^(par|nir)"), mean)) %>%
    ungroup()

  merra_summary %>%
    pivot_longer(c(par_direct:nir_dsf)) %>%
    ggplot() +
    aes(x = value) +
    geom_histogram() +
    facet_wrap(vars(name), scales = "free_x")
}
