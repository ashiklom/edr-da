library(redr)
library(PEcAnRTM)
library(tidyverse)
library(optparse)
library(PEcAn.ED2)
import::from(here, inhere = here)
import::from(lubridate, as_date, year, month, mday)
import::from(progress, progress_bar)

import::from(foreach, foreach, "%dopar%")

import::from(snow, makeCluster)
import::from(doSNOW, registerDoSNOW)

# [1] "OF05"  "IDS36" "SF03"  "BH07"  "AK60"  "OF02"  "BH05" 

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c("--site=OF05", "--ncores=2")
}

argl <- OptionParser() %>%
  add_option("--site", action = "store", type = "character", default = "AK60") %>%
  add_option("--ncores", action = "store", type = "integer", default = 1) %>%
  parse_args(args)
print(argl)

site <- argl$site
ensemble_root <- inhere("ensemble_outputs", "msp_hf20180402")
site_dir <- list.files(ensemble_root, argl$site, full.names = TRUE)
stopifnot(length(site_dir) == 1)

simulate_spec_ensemble <- function(ens, site_dir, by = 10) {
  ens_dir <- file.path(site_dir, sprintf("ens_%03d", ens))
  stopifnot(file.exists(ens_dir))
  message("Running ensemble: ", ens)
  trait_values <- readRDS(file.path(ens_dir, "trait_values.rds"))
  history_files <- list.files(file.path(ens_dir, "out"), "history-S.*\\.h5")
  has_dates <- stringr::str_extract(history_files, "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}") %>%
    lubridate::as_date()
  has_dates <- has_dates[has_dates > lubridate::as_date("2000-01-01")]
  dates <- seq(min(has_dates), max(has_dates), by = by)
  message("Running from ", min(dates) " to ", max(dates))
  message("In total, ", length(dates), " runs.")
  ed2in <- PEcAn.ED2::read_ed2in(file.path(ens_dir, "ED2IN"))
  ed2in$RK4_TOLERANCE <- 1e-5
  run_edr_date2 <- function(date, ...) {
    message("Running date: ", date)
    purrr::safely(run_edr_date)(date, ...)
  }
  PEcAn.logger::logger.setLevel("INFO")
  purrr::map(
    dates,
    run_edr_date2,
    ed2in = ed2in, trait_values = trait_values
  )
}

############################################################
PEcAn.logger::logger.setLevel("INFO")
ens_vec <- seq(1, 50)
cl <- makeCluster(argl$ncores, outfile = "")
registerDoSNOW(cl)
opts <- list(progress = function(n) message("Running ensemble: ", n))
result <- foreach(
  ens = ens_vec,
  .packages = c("magrittr", "redr"),
  .options.snow = opts
) %dopar%
  purrr::safely(simulate_spec_ensemble)(ens, site_dir, by = 5)
saveRDS(result, file.path(site_dir, "spec_sims.rds"))
