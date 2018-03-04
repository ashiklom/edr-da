library(PEcAn.ED2)
library(PEcAn.data.atmosphere)
library(here)
library(tibble)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

PEcAn.logger::logger.setLevel("DEBUG")

############################################################
# User configuration
root <- here()
site_dir <- normalizePath("~/Dropbox/NASA_TE_PEcAn-RTM_Project/Data/NASA_FFT_Project/site_pss_css/")
#root <- normalizePath("~/Projects/nasa-rtm/edr-da")

############################################################
# Additional paths (should be relative, so no user config required)
inputs_dir <- file.path(root, "ed-inputs")
edi_dir <- file.path(inputs_dir, "EDI")
pda_dir <- file.path(root, "ed-outputs", "multi_site_pda")
dir.create(pda_dir, showWarnings = FALSE)

############################################################
# Load site data
sites <- tibble(files = list.files(site_dir)) %>%
  mutate(
    years = map(file.path(site_dir, files), list.files) %>%
      map(~str_match(., "^FFT.([[:digit:]]+)")[,2]) %>%
      map_dbl(~unique(as.numeric(.))),
    prefixes = file.path(site_dir, files, paste0("FFT.", years, ".")),
    veg = map(prefixes, read_ed_veg),
    latitude = map_dbl(veg, "latitude"),
    longitude = map_dbl(veg, "longitude"),
    GFDL_lat = floor(floor(latitude) / 2) + 45 + 1,
    GFDL_lon = floor(floor(longitude) / 2.5) + 1,
    GFDL_dir = file.path(pda_dir, paste("met", GFDL_lat, GFDL_lon, "GFDL", sep = "_"))
  )

############################################################
# Download meteorology data
start_run <- "2006-07-01"
end_run <- "2006-07-14"

gfdl_dat <- distinct(sites, GFDL_lat, GFDL_lon, GFDL_dir) %>%
  mutate(
    GFDL_dat = pmap(
      list(
        lat.in = GFDL_lat,
        lon.in = GFDL_lon,
        outfolder = GFDL_dir
      ),
      download.GFDL,
      start_date = start_run,
      end_date = end_run,
      site_id = NULL
    )
  )

############################################################
# Convert met to ED format
ed_met_dat <- gfdl_dat %>%
  full_join(sites) %>%
  mutate(
    rundir = file.path(pda_dir, files),
    full_name = map_chr(GFDL_dat, "file"),
    ed_met = pmap(
      list(
        outfolder = rundir,
        in.path = dirname(full_name),
        in.prefix = gsub("\\.[[:digit:]]{4}\\.nc$", "", basename(full_name)),
        start_date = start_run,
        end_date = end_run,
        lat = latitude,
        lon = longitude
      ),
      met2model.ED2,
      overwrite = TRUE
    )
  )

############################################################
# Set global ED2IN settings
ed2in <- read_ed2in(system.file("ED2IN.rgit", package = "PEcAn.ED2")) %>%
  modify_ed2in(
    run_name = "EDR multi-site PDA",
    EDI_path = edi_dir,
    output_types = "all",
    runtype = "INITIAL"
  )

############################################################
# Set site-specific ED2IN settings
site_ed2in <- ed_met_dat %>%
  mutate(
    ed2in = pmap(
      list(
        latitude = latitude,
        longitude = longitude,
        start_date = map(ed_met, "startdate"),
        end_date = map(ed_met, "enddate"),
        met_driver = map(ed_met, "file"),
        run_dir = rundir,
        output_dir = file.path(rundir, "output"),
        veg_prefix = prefixes
      ),
      modify_ed2in,
      ed2in = ed2in
    ),
    ed2in_path = file.path(rundir, "ED2IN")
  )

############################################################
# Check ED2IN validity and write to file
checks <- map_lgl(site_ed2in$ed2in, possibly(check_ed2in, FALSE, quiet = FALSE))
walk2(site_ed2in$ed2in, site_ed2in$ed2in_path, write_ed2in)

############################################################
# Run ED at each site
img_path <- "~/Projects/ED2/ed2.simg"
walk(site_ed2in$ed2in_path, ~run_ed_singularity(img_path, ed2in_path = .))
