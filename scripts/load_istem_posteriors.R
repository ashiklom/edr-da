library(purrr)
import::from(here, here)
import::from(stringr, str_remove)

posterior_dir <- here("ed-inputs", "istem-posteriors")
pft_dirs <- list.files(posterior_dir, "temperate.")
pfts <- str_remove(pft_dirs, "\\.IF$")

get_pft_post <- function(posterior_dir, pft) {
  pd <- list.files(posterior_dir, pft)
  stopifnot(length(pd) == 1)
  mcmc_file <- list.files(file.path(posterior_dir, pd), "trait\\.mcmc\\.pda\\..*\\.Rdata") %>%
    tail(1)
  stopifnot(length(mcmc_file) == 1)
  mcmc_fname <- file.path(posterior_dir, pd, mcmc_file)
  load(mcmc_fname)
  trait_samps <- trait.mcmc
  burn <- purrr::map(trait_samps, coda::niter) %>% purrr::map(~floor(. / 2))
  burned <- purrr::map2(trait_samps, burn, ~window(.x, start = .y))
  samps_mat <- purrr::map(burned, as.matrix) %>% purrr::reduce(cbind)
  colnames(samps_mat) <- paste(pft, names(trait_samps), sep = ".")
  samps_mat
}

pft_list <- map(pfts, get_pft_post, posterior_dir = posterior_dir)
pft_mat <- do.call(cbind, pft_list)

saveRDS(pft_mat, file.path(posterior_dir, "processed.rds"))
