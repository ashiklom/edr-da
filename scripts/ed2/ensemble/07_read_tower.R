library(redr)
library(tidyverse)
library(optparse)
library(PEcAn.ED2)
import::from(here, inhere = here)
import::from(progress, progress_bar)
import::from(fs, dir_ls, path_file, dir_create)
import::from(hdf5r, H5File)
import::from(fst, write_fst)

options(prog_default = list(format = ":current/:total (:eta)"))

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c(
    "--site=OF05",
    "--prefix=msp20180402"
  )
}

argl <- OptionParser() %>%
  add_option("--prefix", action = "store", type = "character", default = "msp20180402-nodrought") %>%
  add_option("--site", action = "store", type = "character", default = "AK60") %>%
  parse_args(args)

root_dir <- inhere("ensemble_outputs", argl$prefix)
stopifnot(file.exists(root_dir))
site_dir <- dir_ls(root_dir, regexp = argl$site)
site_full <- basename(site_dir)
stopifnot(length(site_dir) == 1)

read_tower_ensemble <- function(site_dir, ens, lai = FALSE, soil = FALSE) {
  ens_dir <- file.path(site_dir, sprintf("ens_%03d", ens))
  stopifnot(file.exists(ens_dir))

  tfiles <- fs::dir_ls(ens_dir, recursive = TRUE, regexp = "-T-")
  stopifnot(length(tfiles) > 0)
  tdates <- file_date(tfiles)
  tstart <- as.POSIXct(tdates[1], "UTC")
  tend <- as.POSIXct(tail(tdates, 1) + lubridate::years(1) - lubridate::seconds(1), "UTC")
  times <- seq(tstart, tend, "30 mins")

  soil_names <- c(
    "FMEAN_SENSIBLE_GG_PY",
    "FMEAN_SMOIST_GG_PY",
    "FMEAN_SOIL_ENERGY_PY",
    "FMEAN_SOIL_FLIQ_PY",
    "FMEAN_SOIL_MSTPOT_PY",
    "FMEAN_SOIL_TEMP_PY",
    "FMEAN_SOIL_WATER_PY",
    "FMEAN_TRANSLOSS_PY"
  )

  except <- character()
  if (!lai) except <- c(except, "LAI_PY")
  if (!soil) except <- c(except, soil_names)
  tdata_raw <- purrr::map(tfiles, get_h5vars, except = except)

  vector_data <- tdata_raw %>%
    purrr::map(`[`, setdiff(names(tdata_raw[[1]]), c(soil_names, "LAI_PY", "date", "file"))) %>%
    purrr::map_dfr(dplyr::bind_cols) %>%
    dplyr::mutate(times = times, ens = 1)
  out <- vector_data
  if (soil){
    soil_data <- purrr::map(tdata_raw, `[`, soil_names) %>%
      purrr::modify_depth(2, ~split(., col(.))) %>%
      dplyr::bind_rows()
    out <- dplyr::bind_cols(out, soil_data)
  }
  if (lai) {
    lai_data <- purrr::map(tdata_raw, "LAI_PY") %>%
      purrr::map(aperm, c(3, 1, 2)) %>%
      purrr::map(~lapply(seq(dim(.)[1]), function(x) .[x, , ])) %>%
      purrr::flatten()
    out <- tibble::add_column(out, LAI_CO = lai_data)
  }
  out %>%
    dplyr::select(times, ens, dplyr::everything()) %>%
    dplyr::rename_all(~{stringr::str_remove(., "^FMEAN_") %>% stringr::str_remove("_PY$")})
}

nens <- 5
enss <- seq(nens)
all_ens_raw <- map(enss, with_prog(read_tower_ensemble, nens), site_dir = site_dir)
all_ens <- bind_rows(all_ens_raw)

all_ens %>%
  group_by(times)

ggplot(all_ens) +
  aes(x = times, y = LAI) +
  geom_line()

all_nonlist <- select_if(all_ens, ~!is.list(.))
all_list <- select_if(all_ens, is.list)

save_dir <- dir_create(inhere("sync_data", argl$prefix))
write_fst(all_nonlist, file.path(save_dir, paste("analysis", argl$site, "fst", sep = ".")))
saveRDS(all_list, file.path(save_dir, paste("analysis_list", argl$site, "rds", sep = ".")))
