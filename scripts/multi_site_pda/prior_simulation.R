library(PEcAnRTM)
library(PEcAn.ED2)
library(here)
library(purrr)
library(redr)
library(optparse)

parser <- OptionParser() %>%
  add_option("--nsite", action = "store", type = "integer", default = 3L) %>%
  add_option("--nsim", action = "store", type = "integer", default = 20L) %>%
  add_option("--prefix", action = "store", type = "character", default = "multi_site_pda") %>%
  add_option("--geo", action = "store_true", default = FALSE) %>%
  add_option("--hetero", action = "store_true", default = FALSE) %>%
  add_option("--fix_allom2", action = "store_true", default = FALSE) %>%
  add_option("--nprior", action = "store", type = "integer", default = 5000L)

argl <- parse_args(parser)
print(argl)

if (argl$geo) {
  options(
    redr.img_path = NULL,
    redr.edr_exe_path = "/projectnb/dietzelab/ashiklom/ED2/EDR/build/ed_2.1-opt"
  )
} else {
  options(redr.img_path = "~/Projects/ED2/ed2.simg")
}

prefix <- argl$prefix
n_sim <- argl$nsim
site_index <- argl$nsite

site_list <- readLines(here("other_site_data", "site_list"))
stopifnot(site_index <= length(site_list))
site <- site_list[site_index]

pda_dir <- here("ed-outputs", prefix)
stopifnot(file.exists(pda_dir))

prior <- create_prior(fix_allom2 = argl$fix_allom2, heteroskedastic = argl$hetero)
message("Drawing priors")
sims <- map(seq_len(argl$nprior), ~prior$sampler()) %>% do.call(rbind, .)

message("Running simulations for site: ", site)
PEcAn.logger::logger.setLevel("INFO")
result <- site_sensitivity(sims, site, pda_dir, n = n_sim)

saveRDS(result, file.path(pda_dir, paste("prior_sim", site, "rds", sep = ".")))
