library(PEcAnRTM)
library(tidyverse)
library(redr)
import::from(here, inhere = here)
import::from(progress, progress_bar)

data("sensor.rsr", package = "PEcAnRTM")

selected_sites <- readLines(inhere("other_site_data", "selected_sites"))
sync_dir <- inhere("sync_data")

site_spec_data <- list.files(sync_dir, full.names = TRUE)

site2landsat <- function(sitefile, landsat = "landsat7", pb = NULL) {
  on.exit(if (!is.null(pb)) pb$tick())
  d1 <- readRDS(sitefile)
  dgrouped <- d1 %>%
    group_by(date, site, ens_code) %>%
    summarize(landsat = list(spec2landsat(value)))
  dlandsat <- dgrouped %>% unnest()
  dl <- dlandsat %>%
    filter(landsat == !!landsat) %>%
    select(-landsat, -wavelength) %>%
    spread(band, value)
  dlsummary <- dl %>%
    group_by(date, site) %>%
    summarize_at(
      vars(-ens_code),
      funs(Mean = mean, SD = sd)
    )
  dlsummary
}

pb <- progress_bar$new(total = length(site_spec_data), format = "(:current/:total) :eta")
results_raw <- map(site_spec_data, safely(site2landsat), pb = pb)

results_flipped <- transpose(results_raw)
results <- results_flipped$result %>% bind_rows()

results_dir <- inhere("proc_results")
dir.create(results_dir, showWarnings = FALSE)
saveRDS(results, file.path(results_dir, "landsat7_spectra.rds"))
