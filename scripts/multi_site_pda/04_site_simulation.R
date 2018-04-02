library(PEcAnRTM)
library(PEcAn.ED2)
library(here)
library(purrr)

img_path <- NULL
edr_exe_path <- "/projectnb/dietzelab/ashiklom/ED2/EDR/build/ed_2.1-dbg"
prefix <- "multi_site_pda"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args <- c("20", 1:5)
}

n_sim <- as.numeric(args[1])
site_index <- as.numeric(args[-1])

if (system2("hostname", stdout = TRUE) == "ashiklom") {
  img_path <- "~/Projects/ED2/ed2.simg"
  edr_exe_path <- NULL
}

source("scripts/multi_site_pda/site_sensitivity.R")

pda_dir <- here("ed-outputs", "multi_site_pda")

site_dirs <- list.files(pda_dir, "_site_")
ed2in_paths <- list.files(file.path(pda_dir, site_dirs, "edr"), "ED2IN", full.names = TRUE)[site_index]
names(ed2in_paths) <- strsplit(ed2in_paths, "/") %>% map_chr(8)

samples <- readRDS(file.path(pda_dir, "progress.rds"))
burnin <- 1000
waves <- seq(400, 1300, by = 10)

message("Loading observations")
obs_list_full <- map(names(ed2in_paths), possibly(get_obs, NULL), use_wl = waves)
names(obs_list_full) <- names(ed2in_paths)
has_aviris <- !map_lgl(obs_list_full, is.null)
obs_list <- obs_list_full[has_aviris]
sim_sites <- ed2in_paths[names(obs_list)]

message("Running simulations for sites: \n", paste(names(sim_sites), collapse = "\n"))

PEcAn.logger::logger.setLevel("INFO")
message("Running site simulations")
result <- site_sensitivity(samples, sim_sites, burnin, n = n_sim)

suffix <- paste(site_index, collapse = "_")

pdf(paste(prefix, "refl_valid", suffix, "pdf", sep = "."))
pwalk(
  list(
    result = result,
    obs = obs_list,
    main = names(obs_list)
  ),
  possibly(plot_sens, NULL)
)
dev.off()
