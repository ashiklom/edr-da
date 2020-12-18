library(redr)
library(PEcAnRTM)
library(tidyverse)
library(optparse)
library(PEcAn.ED2)
import::from(here, inhere = here)
import::from(lubridate, as_date, year, month, mday)
import::from(progress, progress_bar)

import::from(foreach, foreach, "%do%")

# [1] "OF05"  "IDS36" "SF03"  "BH07"  "AK60"  "OF02"  "BH05" 

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c(
    #"--site=AK60"
    glue("--site={argl$site}"),
    glue("--ens=1"),
    glue("--prefix={argl$outdir}"),
    glue("--by=1")
  )
}

argl <- OptionParser() %>%
  add_option("--site", action = "store", type = "character", default = "AK60") %>%
  add_option("--ens", action = "store", type = "integer", default = 1) %>%
  add_option("--by", action = "store", type = "integer", default = 10) %>%
  add_option("--prefix", action = "store", type = "character", default = "msp20180402") %>%
  parse_args(args)
print(argl)

site <- argl$site
ensemble_root <- inhere("ensemble_outputs", argl$prefix)
site_dir <- list.files(ensemble_root, argl$site, full.names = TRUE)
stopifnot(length(site_dir) == 1)

write_line <- function(a, b, file, append = TRUE, ...) {
  btxt <- paste(b, collapse = ",")
  txt <- paste(a, btxt, sep = ",")
  cat(txt, "\n", file = file, append = append, ...)
  txt
}

simulate_spec_ensemble <- function(ens, site_dir, by = 10) {
  ens_dir <- file.path(site_dir, sprintf("ens_%03d", ens))
  stopifnot(file.exists(ens_dir))
  message("Running ensemble: ", ens)
  trait_values <- readRDS(file.path(ens_dir, "trait_values.rds"))
  history_files <- list.files(file.path(ens_dir, "out"), "history-S.*\\.h5")
  has_dates <- stringr::str_extract(history_files, "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}") %>%
    lubridate::as_date()
  #has_dates <- has_dates[has_dates > lubridate::as_date("2000-01-01")]
  dates <- seq(min(has_dates), max(has_dates), by = by)
  message("Running from ", min(dates), " to ", max(dates))
  message("In total, ", length(dates), " runs.")
  ed2in <- PEcAn.ED2::read_ed2in(file.path(ens_dir, "ED2IN"))
  #ed2in$RK4_TOLERANCE <- 1e-5
  PEcAn.logger::logger.setLevel("INFO")
  outfile <- file.path(ens_dir, "spectra_sims.csv")
  file.remove(outfile)
  file.create(outfile)
  wavecols <- paste0("wave_", 400:2500)
  write_line("date", wavecols, outfile)
  foreach(date = dates, .errorhandling = "pass") %do% {
    message("Running date: ", date)
    spec <- run_edr_date(
      date = date,
      ed2in = ed2in,
      trait_values = trait_values
    )
    write_line(date, spec, outfile)
    spec
  }
}

############################################################
PEcAn.logger::logger.setLevel("INFO")
out <- simulate_spec_ensemble(argl$ens, site_dir, by = argl$by)
saveRDS(out, file.path(site_dir, sprintf("ensemble_%03d.rds", argl$ens)))
