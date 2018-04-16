library(FedData)
library(redr)
library(tidyverse)
import::from(here, here)
import::from(PEcAn.ED2, read_ed_veg)

site_df <- tibble(sites = list.files("sites")) %>%
  mutate(
    full_paths = here("sites", sites),
    prefixes = map_chr(full_paths, get_site_prefix),
    veg = map(prefixes, read_ed_veg),
    latitude = map_dbl(veg, "latitude"),
    longitude = map_dbl(veg, "longitude")
  )

latitude <- site_df$latitude[1]
longitude <- site_df$longitude[1]
label <- site_df$sites[1]

get_ssurgo_point <- function(latitude, longitude, label) {
  wgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  coord <- sp::SpatialPoints(list(x = longitude, y = latitude), proj4string = wgs84)
  out <- get_ssurgo(coord, label = label)
  return(out)
}

site_df2 <- site_df %>%
  mutate(
    ssurgo = pmap(list(latitude, longitude, sites), possibly(get_ssurgo_point, NULL))
  )

# Tables:
# - chorizon: Horizon information
# - cosoilmoist: Component soil moisture
# - cosoiltemp: Component soil temperature
#
# Other notes:
#   - RV (suffix `.r`) is the "representative value" -- generally, the best value to use
#   - Horizons are grouped by component, whose metadata is in the `component` table
get_ssurgo_texture <- function(ssurgo) {
  stopifnot("tabular" %in% names(ssurgo))
  tables <- ssurgo$tabular
  tex_group <- tables$chtexturegrp %>%
    as_tibble() %>%
    filter(grepl("yes", rvindicator, ignore.case = TRUE))
  tex_group %>%
    left_join(tables$chorizon, by = "chkey") %>%
    select(
      component = cokey,
      depth_top_cm = hzdept.r,
      depth_bottom_cm = hzdepb.r,
      texture_code = texture,
      texture_description = texdesc,
      horizon_name = hzname
    ) %>%
    arrange(component) %>%
    filter(component == unique(component)[1]) %>% # HACK: Grab only the first component
    select(-component) %>%
    arrange(depth_bottom_cm)
}

site_soil <- site_df2 %>%
  mutate(
    maptex = map(ssurgo, safely(get_ssurgo_texture)),
    texture = map(maptex, "result")
  ) %>%
  filter(!map_lgl(texture, is.null))

selected <- readLines("other_site_data/selected_sites")
selected_rxp <- paste(selected, collapse = "|")
site_soil %>%
  filter(grepl(selected_rxp, sites)) %>%
  distinct(sites)

site_soil %>%
  select(sites, latitude, longitude, texture) %>%
  unnest(texture) %>%
  write_csv("other_site_data/soil_texture.csv")
