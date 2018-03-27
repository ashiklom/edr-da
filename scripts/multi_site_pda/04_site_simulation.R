library(PEcAnRTM)
library(PEcAn.ED2)
library(here)
library(purrr)
library(foreach)
library(doParallel)

img_path <- NULL
edr_exe_path <- "/projectnb/dietzelab/ashiklom/ED2/EDR/build/ed_2.1-dbg"
source("scripts/multi_site_pda/site_sensitivity.R")

pda_dir <- here("ed-outputs", "multi_site_pda")

site_dirs <- list.files(pda_dir, "_site_")
ed2in_paths <- list.files(file.path(pda_dir, site_dirs, "edr"), "ED2IN", full.names = TRUE)
names(ed2in_paths) <- strsplit(ed2in_paths, "/") %>% map_chr(8)

samples <- readRDS(file.path(pda_dir, "progress.rds"))
burnin <- 1000
ssites <- ed2in_paths[1]
waves <- seq(400, 1300, by = 10)
print(ssites)

message("Loading observations")
obs_list <- map(names(ssites), get_obs, use_wl = waves)

PEcAn.logger::logger.setLevel("INFO")
result <- site_sensitivity(samples, ssites, burnin, n = 50)

r1 <- result[[1]]
rmu <- rowMeans(r1)
rlo <- apply(r1, 1, quantile, 0.025)
rhi <- apply(r1, 1, quantile, 0.975)
