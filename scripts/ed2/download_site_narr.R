library(PEcAn.ED2)
library(purrr)
import::from(PEcAn.data.atmosphere, download.NARR_site)
import::from(magrittr, "%>%")
import::from(here, here)
import::from(stringr, str_match, str_extract)
import::from(lubridate, as_date)

sites <- list.files(here("sites"))
select_sites <- readLines("other_site_data/selected_sites")
sites_sub <- grep(paste(select_sites, collapse = "|"), sites, value = TRUE)
sites_sub_dirs <- here("sites", sites_sub)

for (site_dir in sites_sub_dirs) {
  message("Processing: ", site_dir)
  stopifnot(file.exists(site_dir))
  if (file.exists(file.path(site_dir, "NARR", "NARR.2017.nc"))) {
    message("NARR file exists. Moving to next site.")
    next
  }
  site_files <- list.files(site_dir, "\\.css$")[1]
  stopifnot(length(site_files) == 1)
  prefix <- str_match(site_files, "(^FFT\\.[[:digit:]]{4}\\.)lat.*")[, 2]
  full_prefix <- file.path(site_dir, prefix)
  veg <- read_ed_veg(full_prefix)
  latitude <- veg$latitude
  longitude <- veg$longitude
  year <- as.numeric(str_match(prefix, "^FFT\\.([[:digit:]]{4})")[, 2])
  if (grepl("OF05", site_dir)) {
    year <- 1983
  }
  start_date <- as_date(paste0(year, "-01-01"))
  end_date <- as_date("2017-12-31")
  outfolder <- file.path(site_dir, "NARR")
  metresult <- tryCatch({
    download.NARR_site(outfolder, start_date, end_date, latitude, longitude,
                                  progress = TRUE, parallel = TRUE, ncores = parallel::detectCores())
  }, error = function(e) {message(e); NULL}, 
  interrupt = function(e) {message(e); NULL})
  if (!is.null(metresult)) {
    saveRDS(metresult, file.path(site_dir, "NARR_info.rds"))
  }
}
