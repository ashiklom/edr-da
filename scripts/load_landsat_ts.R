library(tidyverse)
library(lubridate)
import::from(here, inhere = here)

landsat_dir <- inhere("other_site_data", "landsat_time_series")
site_codes <- list.files(landsat_dir)

read_landsat_site <- function(site_code) {
  site_path <- file.path(landsat_dir, site_code)
  band_files <- list.files(site_path, full.names = TRUE)
  data_list <- map(band_files, read_csv, col_types = "cn")
  data_df <- reduce(data_list, full_join, by = "system:time_start")
  data_df %>%
    mutate(date = mdy(`system:time_start`), site = site_code) %>%
    select(site, date, BLUE:SWIR2)
}

landsat_data <- map_dfr(site_codes, read_landsat_site) %>%
  filter_at(vars(-site, -date), any_vars(!is.na(.)))

write_csv(landsat_data, inhere("other_site_data", "landsat_data.csv"))
