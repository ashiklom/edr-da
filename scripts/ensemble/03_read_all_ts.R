library(tidyverse)
library(optparse)
import::from(here, inhere = here)

# [1] "OF05"  "IDS36" "SF03"  "BH07"  "AK60"  "OF02"  "BH05" 

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c("--site=BH05")
  #rmote::start_rmote(port = 4433)
}
argl <- OptionParser() %>%
  add_option("--site", action = "store", type = "character", default = "AK60") %>%
  parse_args(args)
print(argl)

site <- argl$site
root_dir <- inhere("ensemble_outputs", "msp_hf20180402")
site_dir <- list.files(root_dir, site, full.names = TRUE)
stopifnot(length(site_dir) == 1)

message("Finding CSV files")
csv_files_raw <- system2("find", c(site_dir, "-name", "spectra_sims.csv"),
                         stdout = TRUE)
#csv_files <- list.files(site_dir, "spectra_sims.csv",
                        #full.names = TRUE, recursive = TRUE)

message("Reading CSV files")
colspec <- cols(date = col_date(), .default = col_double())
csv_list <- map(csv_files, read_csv, col_types = colspec)

message("Row-binding CSV files")
allspec <- bind_rows(csv_list, .id = "id")

ens_ids <- str_extract(csv_files, "ens_[[:digit:]]{3}")
ens_df <- tibble(ens_code = ens_ids, id = seq_along(ens_ids))

message("Gathering CSV files")
longspec <- allspec %>%
  mutate(id = as.integer(id), site = site) %>%
  left_join(ens_df) %>%
  gather(wavelength, value, -id, -date, -ens_code, -site) %>%
  mutate(
    wavelength = str_remove(wavelength, "wave_") %>% as.numeric
  )

syncdir <- inhere("sync_data")
dir.create(syncdir, showWarnings = FALSE)
message("Saving to: ", syncdir)
saveRDS(longspec, file.path(syncdir, paste("ts_spectra", site, "rds", sep = ".")))
