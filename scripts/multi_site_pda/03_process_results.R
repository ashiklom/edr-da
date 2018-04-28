library(redr)
library(BayesianTools)
library(optparse)
library(tidyverse)
library(PEcAnRTM)
import::from(here, inhere = here)
import::from(progress, progress_bar)

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  args <- c("--prefix=msp_hf20180402",
            "--burnin=80000")
}

argl <- OptionParser() %>%
  add_option("--prefix", action = "store", type = "character", default = "msp20180402") %>%
  add_option("--hetero", action = "store_true", default = FALSE) %>%
  add_option("--fix_allom2", action = "store_true", default = FALSE) %>%
  add_option("--final", action = "store_true", default = FALSE) %>%
  add_option("--burnin", action = "store", type = "integer", default = 8000L) %>%
  parse_args(args) 

if (grepl("_h?f", argl$prefix)) argl$fix_allom2 <- TRUE
if (grepl("_h", argl$prefix)) argl$hetero <- TRUE
print(argl)

pda_dir <- inhere("ed-outputs", argl$prefix)
stopifnot(file.exists(pda_dir))

resultfile <- file.path(pda_dir, ifelse(argl$final, "results.rds", "progress.rds"))
stopifnot(file.exists(resultfile))

result <- readRDS(resultfile)
figdir <- inhere("msp_figures", argl$prefix)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE) 

message("Generating trace plots")
pdf(file.path(figdir, "traces.pdf"))
tracePlot(result)
dev.off()
message("Done!")

samples_burned <- getSample(result, start = argl$burnin, coda = TRUE)
params <- coda::varnames(samples_burned)

subparams <- params[!grepl("residual", params)] %>% str_split("\\.")
pft <- map_chr(subparams, 2)
upft <- unique(pft)
variable <- map_chr(subparams, 3)
uvariable <- unique(variable)

samples_mat_full <- as.matrix(samples_burned)
samples_mat <- samples_mat_full[, !grepl("residual", params)]
samples_bypft <- map(upft, ~samples_mat[, pft == .]) %>%
  map(`colnames<-`, uvariable)
samples_byvar <- map(uvariable, ~samples_mat[, variable == .]) %>%
  map(`colnames<-`, upft)

message("Generating pairs plots grouped by PFT")
pdf(file.path(figdir, "pairs.bypft.pdf"))
walk2(upft, samples_bypft, ~pairs(.y, main = .x, pch = "."))
dev.off()

message("Generating pairs plots grouped by variable")
pdf(file.path(figdir, "pairs.byvar.pdf"))
walk2(uvariable, samples_byvar, ~pairs(.y, main = .x, pch = "."))
dev.off()

message("Summarizing samples")
results_summary <- summarize_samples(samples_burned) %>%
  mutate(type = "posterior")

message("Drawing priors")
prior <- create_prior(fix_allom2 = argl$fix_allom2, heteroskedastic = argl$hetero)
prior_draws <- map(seq_len(2000), ~prior$sampler()) %>% invoke(rbind, .)
message("Done drawing priors")

prior_summary <- tibble(
  params = colnames(prior_draws),
  `Mean` = colMeans(prior_draws),
  `2.5%` = apply(prior_draws, 2, quantile, 0.025),
  `97.5%` = apply(prior_draws, 2, quantile, 0.975),
  type = "prior"
) %>%
  split_params("params")

all_summary <- bind_rows(results_summary, prior_summary) %>%
  filter(!grepl("residual", variable)) %>%
  mutate(
    pft = factor(pft, levels = unique(pft)),
    type = factor(type, c("prior", "posterior")),
    variable = factor(variable, levels = unique(variable)),
    pft_type = interaction(type, pft)
  )

summary_dir <- inhere("sync_data", argl$prefix)
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(all_summary, file.path(summary_dir, "pda_summary.rds"))

message("Generating summary plot")
summary_plot <- ggplot(all_summary) +
  aes(x = pft_type, y = Mean, ymin = `2.5%`, ymax = `97.5%`,
      color = pft, linetype = type) +
  geom_pointrange() +
  facet_wrap(~ variable, scales = "free") +
  theme(axis.text.x = element_blank())

ggsave(file.path(figdir, "summary.pdf"), summary_plot)

message("Generating density plots")
samples_matrix <- as.matrix(samples_burned)
post_densities <- col_densities(samples_matrix)
prior_densities <- col_densities(prior_draws[, params])

pdf(file.path(figdir, "densplots.pdf"))
pwalk(
  list(post_densities, prior_densities, params),
  prior_posterior_density
)
dev.off()

sitelist <- readLines(inhere("other_site_data", "site_list"))

message("LAI histograms for sites")
get_lai <- function(site, pda_dir, burnin, pb = NULL, ...) {
  if (!is.null(pb)) pb$tick()
  lai_file <- file.path(pda_dir, site, "edr", "lai_store")
  stopifnot(file.exists(lai_file))
  nlines <- wcl(lai_file)
  skip <- ifelse(nlines < burnin, 0, burnin)
  data.table::fread(lai_file, header = FALSE, blank.lines.skip = TRUE, skip = skip) %>%
    as.matrix()
}

pb <- progress::progress_bar$new(total = length(sitelist))
lai_list <- map(
  sitelist,
  get_lai,
  pda_dir = pda_dir,
  burnin = argl$burnin,
  pb = pb
)
names(lai_list) <- sitelist
saveRDS(lai_list, file.path(summary_dir, "lai_values.rds"))

lai_sums <- map(lai_list, rowSums)
pdf(file.path(figdir, "lai_hist.pdf"))
iwalk(lai_sums, ~hist(.x, xlab = "Total LAI", main = .y))
dev.off()

pb <- progress_bar$new(total = length(sitelist), format = ":current/:total (ETA: :eta)")
spec_history <- map(
  sitelist,
  ~{pb$tick(); safely(read_spectra_history)(., pda_dir, argl$burnin)}
)
names(spec_history) <- sitelist
for (s in sitelist) {
  print(s)
  saveRDS(spec_history[[s]], file.path(summary_dir, paste0("spec_history.", s, ".rds")))
} 

message("Spectral confidence intervals for sites")
pb <- progress::progress_bar$new(
  total = length(sitelist),
  format = ":current/:total (ETA: :eta)"
)
pdf(file.path(figdir, "refl_validation.pdf"))
walk(
  sitelist,
  possibly(plot_site_spectra, NULL, quiet = FALSE),
  pda_dir = pda_dir,
  burnin = argl$burnin,
  pb = pb
)
dev.off()
