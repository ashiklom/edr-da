library(tidyverse)
library(optparse)
import::from(here, inhere = here)
import::from(fs, dir_ls, path_file)
import::from(progress, progress_bar)

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
root_dir <- inhere("ensemble_outputs", "msp_hf20180402.2018-04-18")
site_dir <- list.files(root_dir, site, full.names = TRUE)
stopifnot(length(site_dir) == 1)

#rds_files <- dir_ls(site_dir, glob = "*ensemble_*.rds")
#names(rds_files) <- path_file(rds_files) %>%
  #str_match("ensemble_([[:digit:]]{3}).rds") %>%
  #.[,2]
#raw_in <- map(rds_files, with_prog(readRDS, length(rds_files)))

#lai_vals <- modify_depth(raw_in, 2, attr_getter("LAI")) %>%
  #map(~do.call(cbind, .))

message("Finding CSV files")
csv_files <- dir_ls(site_dir, recursive = TRUE, glob = "*spectra_sims.csv")

message("Reading CSV files")
colspec <- cols(date = col_date(), .default = col_double())
csv_list <- map(csv_files, with_prog(read_csv, length(csv_files)), col_types = colspec)

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
