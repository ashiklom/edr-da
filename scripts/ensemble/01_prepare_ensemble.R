library(redr)
library(tidyverse)
library(optparse)
library(PEcAn.ED2)
import::from(here, here)
import::from(BayesianTools, getSample)
import::from(PEcAnRTM, params2edr)
import::from(progress, progress_bar)

# [1] "OF05"  "IDS36" "SF03"  "BH07"  "AK60"  "OF02"  "BH05" 

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c(
    "--burnin=40000",
    "--nens=50",
    "--site=OF05_site_1-25710"
    #"--site=IDS36_site_1-25686"
    #"--site=SF03_site_1-25721"
    #"--site=BH07_site_1-25669"
    #"--site=AK60_site_1-25674"
    #"--site=OF02_site_1-25708"
    #"--site=BH05_site_1-25667"
  )
}

argl <- OptionParser() %>%
  add_option("--prefix", action = "store", type = "character", default = "msp_hf20180402") %>%
  add_option("--outdir", action = "store", type = "character", default = NULL) %>%
  add_option("--site", action = "store", type = "character", default = "BH02_site_1-25665") %>%
  add_option("--final", action = "store_true", default = FALSE) %>%
  add_option("--initial", action = "store_true", default = FALSE) %>%
  add_option("--burnin", action = "store", type = "integer", default = 8000) %>%
  add_option("--nens", action = "store", type = "integer", default = 50) %>%
  parse_args(args)

print(argl)

############################################################
# Set up paths
############################################################
site <- argl$site
is_initial <- argl$initial

ens_root <- here("ensemble_outputs")
dir.create(ens_root, showWarnings = FALSE)

pda_dir <- here("ed-outputs", argl$prefix)
stopifnot(file.exists(pda_dir))

outdir <- if (is.null(argl$outdir)) argl$prefix else argl$outdir
ens_dir <- file.path(ens_root, outdir)
dir.create(ens_dir, showWarnings = FALSE)

ens_site_dir <- file.path(ens_dir, site)
dir.create(ens_site_dir, showWarnings = FALSE)

site_dir <- file.path(pda_dir, site)
stopifnot(file.exists(site_dir))

# Parse prefix for other flags
hetero <- grepl("_.?h", argl$prefix)
fix_allom2 <- grepl("_.?f", argl$prefix)

############################################################
# Prepare ensemble trait samples
############################################################
samp_base <- if (argl$final) "results.rds" else "progress.rds"
samp_fname <- file.path(pda_dir, samp_base)
stopifnot(file.exists(samp_fname))
edr_samples_all <- readRDS(samp_fname)
edr_samples <- getSample(edr_samples_all, start = argl$burnin)
other_posteriors <- readRDS(here("ed-inputs", "istem-posteriors", "processed.rds"))

edr_sub <- edr_samples[sample.int(nrow(edr_samples), argl$nens), ]
other_sub <- other_posteriors[sample.int(nrow(other_posteriors), argl$nens), ]
combined <- cbind(edr_sub, other_sub)

sum4ed <- function(params, vis = 400:700, nir = 701:1300) {
  pedr <- params2edr(params)
  pspec <- pedr$spectra_list
  traits <- pedr$trait.values
  pvis <- purrr::map(pspec, ~colMeans(.[[vis]])) %>%
    purrr::map(setNames, c("leaf_reflect_vis", "leaf_trans_vis")) %>%
    .[names(traits)]
  pnir <- purrr::map(pspec, ~colMeans(.[[nir]])) %>%
    purrr::map(setNames, c("leaf_reflect_nir", "leaf_trans_nir")) %>%
    .[names(traits)]
  stopifnot(
    all(names(traits) == names(pvis)),
    all(names(traits) == names(pnir))
  )
  traits %>%
    purrr::map2(pvis, c) %>%
    purrr::map2(pnir, c)
}

ensemble_trait_list <- apply(combined, 1, sum4ed)

############################################################
# Prepare meteorology
############################################################
site_input_dir <- here("sites", site)
stopifnot(file.exists(site_input_dir))

site_info <- tibble(site_files = list.files(site_input_dir, "\\.css$")[1]) %>%
  mutate(
    prefix = str_match(site_files, "(^.*)lat.*")[, 2],
    full_prefix = file.path(site_input_dir, prefix),
    year = str_extract(prefix, "[[:digit:]]{4}") %>% as.numeric(),
    year = if_else(is_initial, 1983, year),
    veg = map(full_prefix, read_ed_veg),
    latitude = map_dbl(veg, "latitude"),
    longitude = map_dbl(veg, "longitude"),
    start_date = paste0(year, "-07-01"),
    end_date = "2017-12-31"
  )

ed_met_dir <- file.path(site_input_dir, "ED_NARR")

if (!file.exists(file.path(ed_met_dir, "2017DEC.h5"))) {
  met_info <- met2model.ED2(
    in.path = file.path(site_input_dir, "NARR"),
    in.prefix = "NARR",
    outfolder = ed_met_dir,
    start_date = site_info$start_date,
    end_date = site_info$end_date,
    lat = site_info$latitude,
    lon = site_info$longitude
  )
} else {
  message("Found existing ED met. No need to do met2model.")
  met_info <- tibble(
    file = file.path(ed_met_dir, "ED_MET_DRIVER_HEADER"),
    host = "",
    mimetype = "text/plain",
    formatname = "ed.met_driver_header files format",
    startdate = as.POSIXct(site_info$start_date, tz = "UTC"),
    enddate = as.POSIXct(site_info$end_date, tz = "UTC"),
    dbfile.name = "ED_MET_DRIVER_HEADER"
  )
}

############################################################
# Add soil information to vegetation site file
############################################################
soil_all <- read_csv("other_site_data/soil_texture.csv")
soil_site <- soil_all %>%
  filter(sites == site)

split_soil <- function(soil_df) {
  toprow <- botrow <- soil_df
  toprow$depth_bottom_cm <- toprow$depth_bottom_cm / 2
  botrow$depth_top_cm <- toprow$depth_bottom_cm
  dplyr::bind_rows(toprow, botrow)
}

nsoil <- nrow(soil_site)
# HACK: ED complains if there's only one soil layer.
# So, if there's only one, split it into two identical layers.
if (nsoil == 1) {
  soil_site <- split_soil(soil_site)
  nsoil <- nrow(soil_site)
}

add_site_soil <- function(veg, soil_df) {
  site <- veg$site
  nsite <- attr(site, "nsite")
  site$soil <- NULL
  texture_vec <- soil_df %>%
    left_join(texture_codes, by = "texture_code") %>%
    select(ed_texture_id) %>%
    purrr::transpose() %>%
    purrr::flatten() %>%
    rev() %>%
    setNames(paste0("soil", seq(length(.))))
  site <- tibble::add_column(site, !!!texture_vec)
  attr(site, "file_format") <- 3
  attr(site, "nsite") <- nsite
  veg$site <- site
  veg
}

veg_mod <- add_site_soil(site_info$veg[[1]], soil_site)
veg_prefix <- file.path(ens_site_dir, "FFT.ensemble")
write_ed_veg(veg_mod, veg_prefix)

############################################################
# Prepare ED2IN
############################################################
site_ed2in <- file.path(site_dir, "ED2IN")
ed2in_temp <- read_ed2in(site_ed2in) %>%
  modify_ed2in(
    veg_prefix = veg_prefix,
    latitude = site_info$latitude,
    longitude = site_info$longitude,
    met_driver = met_info$file,
    start_date = site_info$start_date,
    end_date = site_info$end_date,
    EDI_path = here("ed-inputs", "EDI"),
    output_types = c("monthly", "restart"),
    runtype = "INITIAL",
    IED_INIT_MODE = if_else(is_initial, 0, 3),
    NZG = nsoil,  # number of soil layers
    SLZ = rev(-soil_site$depth_bottom_cm / 100),  # soil layer depths
    SLMSTR = rep(0.65, nsoil),
    STGOFF = rep(0, nsoil),
    UNITSTATE = 1,
    FRQSTATE = 1,
    DTLSM = 600,  # Lower these values to reduce the likelihood of integration failures
    RADFRQ = 600,
    RK4_TOLERANCE = 1e-4,
    INTEGRATION_SCHEME = 1,
    INCLUDE_THESE_PFT = pft_lookup$pft_num  # Limit to PDA PFTs
  )

pb <- progress_bar$new(total = argl$nens)
for (i in seq_len(argl$nens)) {
  pb$tick()
  run_dir <- file.path(ens_site_dir, sprintf("ens_%03d", i))
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  out_dir <- file.path(run_dir, "out")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ed2in_ens <- ed2in_temp %>%
    modify_ed2in(
      output_dir = out_dir,
      run_dir = run_dir,
      pecan_defaults = FALSE
    )
  ed2in_path <- file.path(run_dir, "ED2IN")
  write_ed2in(ed2in_ens, ed2in_path)
  trait_values <- ensemble_trait_list[[i]]
  saveRDS(trait_values, file.path(run_dir, "trait_values.rds"))
  xml_path <- file.path(run_dir, "config.xml")
  PEcAn.logger::logger.setLevel("INFO")
  xml <- PEcAn.ED2::write.config.xml.ED2(
    defaults = list(),
    settings = list(revision = "git", config.header = NULL),
    trait.values = trait_values
  )
  PREFIX_XML <- "<?xml version=\"1.0\"?>\n<!DOCTYPE config SYSTEM \"ed.dtd\">\n"
  XML::saveXML(xml, file = xml_path, indent = TRUE, prefix = PREFIX_XML)
}

############################################################
# Write submission script
############################################################
ed_exe_path <- "/projectnb/dietzelab/ashiklom/ED2/ED/build/ed_2.1-dbg"
log_path <- "logs_ed"
dir.create(log_path, showWarnings = FALSE)
submit_script <- c(
  "#!/bin/bash",
  "#$ -q \"geo*\"",
  "#$ -pe omp 4",
  "#$ -j y",
  paste("#$ -o", log_path),
  paste0("#$ -t 1-", argl$nens),
  paste("#$ -N", site),
  "",
  "printf -v ENSDIR \"ens_%03d\" $SGE_TASK_ID",
  "",
  paste(
    "/usr3/graduate/ashiklom/.singularity/sexec",
    ed_exe_path,
    "-f",
    file.path(ens_site_dir, "$ENSDIR", "ED2IN")
  )
)
writeLines(submit_script, paste0("qsub_ed_", site, ".sh"))

if (interactive()) {
  # Test one run
  system2(ed_exe_path, c("-s", "-f", ed2in_path))
}

if (FALSE) {
  devtools::install("~/dietzelab/pecan/models/ed")
}
