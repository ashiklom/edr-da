library(redr)
library(tidyverse)
library(optparse)
import::from(here, inhere = here)
import::from(progress, progress_bar)
import::from(fs, dir_ls, path_file)
import::from(hdf5r, H5File)
import::from(lubridate, as_date, month, mday)
import::from(udunits2, ud.convert)

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c(
    "--site=OF02",
    "--ens=1"
  )
}
argl <- OptionParser() %>%
  add_option("--site", action = "store", type = "character", default = "OF02") %>%
  add_option("--ens", action = "store", type = "integer", default = 1) %>%
  parse_args(args)
print(argl)

if (FALSE) {
  var_table <- read_csv("inst/ed_state_vars.csv")
  var_table %>%
    filter(type == "PFT") %>%
    filter(grepl("resp", variable, ignore.case = TRUE)) %>%
    select(variable, group_type) %>%
    print(n = Inf)
}

site_code <- argl$site
ens <- argl$ens

root_dir <- inhere("ensemble_outputs", "msp_hf20180402")
site_dir <- list.files(root_dir, site_code, full.names = TRUE)
ens_dir <- file.path(site_dir, sprintf("ens_%03d", ens))
run_dir <- file.path(ens_dir, "out")
stopifnot(file.exists(site_dir), file.exists(ens_dir), file.exists(run_dir))

histfiles <- dir_ls(run_dir, glob = "*history-S*.h5")
dates <- path_file(histfiles) %>%
  str_extract("[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}") %>%
  as_date()
names(histfiles) <- dates
stopifnot(length(histfiles) > 0)

readhist <- function(hfile, pb = NULL) {
  on.exit(if (!is.null(pb)) pb$tick())
  safely(read_ed_history)(hfile)
}
pb <- progress_bar$new(total = length(histfiles), format = ":current/:total (:eta)")
hist_raw <- map(histfiles, readhist, pb = pb) %>% transpose()
hist_result <- map_if(hist_raw$result, ~!is.null(.), ~.)

add_ens <- . %>% mutate(ens = ens, date = as_date(date)) %>% select(ens, date, everything())

scalars <- map_dfr(hist_result, "scalar", .id = "date") %>%
  spread(variable, value) %>%
  add_ens
cohorts <- map_dfr(hist_result, "cohort", .id = "date") %>% add_ens
pfts <- map_dfr(hist_result, "pft", .id = "date") %>% add_ens
soils <- map_dfr(hist_result, "soil", .id = "date") %>% add_ens
other <- map(hist_result, "other")
all_results <- list(
  ens = ens,
  scalars = scalars,
  cohorts = cohorts,
  pfts = pfts,
  soils = soils,
  other = other
)
saveRDS(all_results, file.path(ens_dir, "history_results.rds"))
