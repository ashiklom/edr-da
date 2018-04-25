library(tidyverse)
library(PEcAn.ED2)
import::from(here, inhere = here)

site_data <- tibble(
  site_bety = readLines(inhere("other_site_data", "site_list"))
) %>%
  mutate(
    site = str_extract(site_bety, "^.*?_"),
    site_dir = inhere("sites", site_bety),
    css_file = map(site_dir, list.files, "\\.css$", full.names = TRUE) %>% map_chr(1),
    prefix = str_match(css_file, "(.*)lat")[, 2],
    year = str_match(prefix, "FFT\\.([[:digit:]]{4})\\.")[, 2] %>% as.numeric(),
    veg = map(prefix, read_ed_veg),
    latitude = map_dbl(veg, "latitude"),
    longitude = map_dbl(veg, "longitude")
  ) %>%
  select(-veg, -css_file)

write_csv(site_data, "other_site_data/site_data.csv")
