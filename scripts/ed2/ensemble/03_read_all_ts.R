library(tidyverse)
library(optparse)
library(redr)
library(fst)
import::from(here, inhere = here)
import::from(fs, dir_ls, path_file, dir_create)
import::from(progress, progress_bar)

# [1] "OF05"  "IDS36" "SF03"  "BH07"  "AK60"  "OF02"  "BH05" 

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c("--site=BH05")
  #rmote::start_rmote(port = 4433)
}
argl <- OptionParser() %>%
  add_option("--site", action = "store", type = "character", default = "AK60") %>%
  add_option("--prefix", action = "store", type = "character", default = "msp20180402") %>%
  parse_args(args)
print(argl)

site <- argl$site
prefix <- argl$prefix

#site_list <- readLines("run_sites")
#for (site in site_list) {
  #message("site: ", site)

root_dir <- inhere("ensemble_outputs", argl$prefix)
site_dir <- list.files(root_dir, site, full.names = TRUE)
stopifnot(length(site_dir) == 1)

message("Finding CSV files")
csv_files <- dir_ls(site_dir, recursive = TRUE, regexp = "spectra_sims\\.csv")

message("Reading CSV files")
colspec <- cols(date = col_date(), .default = col_double())
csv_list <- map(
  csv_files,
  with_prog(read_csv, length(csv_files)),
  col_types = colspec,
  progress = FALSE
)

message("Row-binding CSV files")
allspec <- bind_rows(csv_list, .id = "id") %>%
  mutate(
    ens = str_extract(id, "ens_[[:digit:]]{3}") %>%
      str_remove("ens_") %>%
      as.integer(),
    site = site
  ) %>%
  select(ens, site, everything(), -id)

syncdir <- inhere("sync_data", prefix)
dir_create(syncdir)
message("Saving to: ", syncdir)
write_fst(allspec, file.path(syncdir, paste("ts_spectra", site, "fst", sep = ".")))
#}
