library(optparse)
library(stringr)
import::from(magrittr, "%>%")
import::from(here, inhere = here)
import::from(lubridate, as_date, year, month, mday)
import::from(PEcAn.ED2, read_ed2in, write_ed2in, modify_ed2in)

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c()
}

argl <- OptionParser() %>%
  add_option("--prefix", action = "store", type = "character", default = "msp_hf20180402") %>%
  add_option("--outdir", action = "store", type = "character", default = NULL) %>%
  add_option("--site", action = "store", type = "character", default = "AK60_site_1-25674") %>%
  add_option("--ens", action = "store", type = "integer", default = 1) %>%
  parse_args(args)
str(argl)

site <- argl$site
ens <- argl$ens

ens_root <- inhere("ensemble_outputs")
outdir <- if (is.null(argl$outdir)) argl$prefix else argl$outdir
ens_dir <- file.path(ens_root, outdir)
ens_site_dir <- file.path(ens_dir, site)
stopifnot(file.exists(ens_site_dir))
run_dir <- file.path(ens_site_dir, sprintf("ens_%03d", ens))
stopifnot(file.exists(run_dir))
out_dir <- file.path(run_dir, "out")
stopifnot(file.exists(out_dir))

all_histfiles <- list.files(out_dir, "history-S")
hist_dates <- str_extract(all_histfiles, "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}") %>%
  as_date()
last_date <- max(hist_dates)

ed2in_old <- read_ed2in(file.path(run_dir, "ED2IN"))
ed2in_new <- modify_ed2in(
  ed2in_old,
  runtype = "HISTORY",
  start_date = last_date,
  pecan_defaults = FALSE,
  SFILIN = file.path(out_dir, "history-S")
)
write_ed2in(ed2in_new, file.path(run_dir, "ED2IN_restart"))
