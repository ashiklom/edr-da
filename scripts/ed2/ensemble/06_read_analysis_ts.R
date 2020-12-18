library(redr)
library(tidyverse)
library(optparse)
import::from(here, inhere = here)
import::from(progress, progress_bar)
import::from(fs, dir_ls, path_file)
import::from(hdf5r, H5File)
import::from(lubridate, as_date, month, mday)
import::from(udunits2, ud.convert)
import::from(glue, glue)

options(prog_default = list(format = ":current/:total (:eta)"))

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c(
    #"--site=BH07",
    "--site=OF05",
    #"--ens"
    "--ens=1",
    #"--prefix=msp20180402-nodrought"
    "--prefix=msp20180402"
    #glue("--site={argl$site}"),
    #glue("--prefix={argl$outdir}")
  )
}
argl <- OptionParser() %>%
  add_option("--site", action = "store", type = "character", default = "OF02") %>%
  add_option("--ens", action = "store", type = "integer", default = 1) %>%
  add_option("--prefix", action = "store", type = "character", default = "msp20180402") %>%
  parse_args(args)
print(argl)

site <- argl$site
ens <- argl$ens

root_dir <- inhere("ensemble_outputs", argl$prefix)
site_dir <- list.files(root_dir, site, full.names = TRUE)

ens_dir <- file.path(site_dir, sprintf("ens_%03d", ens))
run_dir <- file.path(ens_dir, "out")
stopifnot(file.exists(site_dir), file.exists(ens_dir), file.exists(run_dir))

anfiles <- dir_ls(run_dir, regexp = "analysis-E.*\\.h5")
dates <- file_date(anfiles)
if (interactive()) {
  anfiles <- anfiles[dates < "2012-01-01"]
  dates <- dates[dates < "2012-01-01"]
}
nfile <- length(anfiles)
stopifnot(nfile > 0)

out <- map(anfiles, with_prog(safely(get_h5vars), nfile), vars = NULL)

if (interactive()) {

  lai_co <- map(out, c("result", "LAI_CO"))
  lai_df <- tibble(
    date = dates,
    lai = map_dbl(lai_co, sum)
  )

  ggplot(lai_df) +
    aes(x = date, y = lai) +
    geom_line() +
    ggtitle(site, ens)

  slz <- map(out, c("result", "SLZ"))[[1]]
  soil_water_l <- map(out, c("result", "MMEAN_SOIL_WATER_PY"))
  soil_water <- soil_water_l %>%
    bind_cols() %>%
    mutate(slz = slz) %>%
    gather(file, soil_water, -slz) %>%
    mutate(date = file_date(file))

  ggplot(soil_water) +
    aes(x = date, y = soil_water) +
    geom_line() +
    facet_grid(. ~ slz) +
    ggtitle(site, ens)

  prod_dat <- tibble(
    date = dates,
    gpp = map_dbl(out, c("result", "MMEAN_GPP_PY")),
    nep = map_dbl(out, c("result", "MMEAN_NEP_PY")),
    npp = map_dbl(out, c("result", "MMEAN_NPP_PY")),
    tot_resp = map_dbl(out, c("result", "MMEAN_PLRESP_PY")),
    atm_temp = map_dbl(out, c("result", "MMEAN_ATM_TEMP_PY")),
    atm_vpd = map_dbl(out, c("result", "MMEAN_ATM_VPDEF_PY")),
    avail_water = map_dbl(out, c("result", "MMEAN_AVAILABLE_WATER_PY"))
  ) %>%
    left_join(lai_df) %>%
    left_join(soil_water %>% group_by(date) %>% summarize(soil_water = mean(soil_water)))

  prod_dat %>%
    gather(variable, value, -date) %>%
    ggplot() +
    aes(x = date, y = value) +
    geom_line() +
    facet_grid(variable ~ . , scales = "free") +
    xlim(as_date(c("2000-01-01", "2008-01-01")))
    #coord_cartesian(xlim = as_date(c("2000-01-01", "2010-01-01")))
    #coord_cartesian(xlim = as_date(c("2000-01-01", "2005-12-31")))

}

if (FALSE) {
  var_table <- read_csv("inst/ed_state_vars.csv")
  var_table %>%
    filter(type == "PFT") %>%
    filter(grepl("resp", variable, ignore.case = TRUE)) %>%
    select(variable, group_type) %>%
    print(n = Inf)
}

saveRDS(out, file.path(ens_dir, "monthly_output.rds"))
