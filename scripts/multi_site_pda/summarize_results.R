summarize_results <- function(samples, prefix, burnin = 30000, sites = character(), site_iter = 50) {

  pdf(paste(prefix, "traces", "pdf", sep = "."))
  BayesianTools::tracePlot(samples)
  dev.off()

  samples_burned <- getSample(samples, start = burnin, coda = TRUE)
  params <- coda::varnames(samples_burned)

  subparams <- params[params != "residual"] %>%
    str_split("\\.")
  pft <- map_chr(subparams, 2)
  upft <- unique(pft)
  variable <- map_chr(subparams, 3)
  uvariable <- unique(variable)

  samples_mat_full <- as.matrix(samples_burned)
  samples_mat <- samples_mat_full[, params != "residual"]
  samples_bypft <- map(upft, ~samples_mat[, pft == .]) %>%
    map(`colnames<-`, uvariable)
  samples_byvar <- map(uvariable, ~samples_mat[, variable == .]) %>%
    map(`colnames<-`, upft)

  message("Generating pairs plots grouped by PFT")
  pdf(paste(prefix, "pairs", "bypft", "pdf", sep = "."))
  walk2(upft, samples_bypft, ~pairs(.y, main = .x, pch = "."))
  dev.off()

  message("Generating pairs plots grouped by variable")
  pdf(paste(prefix, "pairs", "byvariable", "pdf", sep = "."))
  walk2(uvariable, samples_byvar, ~pairs(.y, main = .x, pch = "."))
  dev.off()

  message("Summarizing samples")
  samples_summary_raw <- summary(samples_burned)

  samples_summary <- samples_summary_raw %>%
    .[c("quantiles", "statistics")] %>%
    map(~as_tibble(., rownames = "parameter")) %>%
    reduce(left_join) %>%
    split_params("parameter") %>%
    mutate(type = "posterior")

  all_summary <- bind_rows(samples_summary, prior_summary) %>%
    filter(variable != "residual") %>%
    mutate(
      pft = factor(pft, levels = unique(pft)),
      variable = factor(variable, levels = unique(variable)),
      pft_type = interaction(type, pft)
    )

  # Summary of parameter estimates
  message("Generating summary plot")
  summary_plot <- ggplot(all_summary) +
    aes(x = pft_type, y = Mean, ymin = `2.5%`, ymax = `97.5%`,
        color = pft, linetype = type) +
    geom_pointrange() +
    facet_wrap(~ variable, scales = "free") +
    theme(axis.text.x = element_blank())

  ggsave(paste(prefix, "summary", "pdf", sep = "."), summary_plot)

  # Density plot
  message("Generating density plots")
  samples_matrix <- as.matrix(samples_burned)
  densities <- apply(samples_matrix, 2, density)
  prior_densities <- apply(prior_draws[, params], 2, density)

  pdf(paste(prefix, "densities", "pdf", sep = "."))
  pwalk(
    list(densities, prior_densities, params),
    prior_posterior_density
  )
  dev.off()
}

prior_posterior_density <- function(post, pri, param, ...) {
  xrange <- range(c(post$x, pri$x))
  yrange <- range(c(post$y, pri$y))
  plot(0, 0, type = "n", xlim = xrange, ylim = yrange,
       main = param, xlab = "value", ylab = "density")
  lines(post$x, post$y, col = "black")
  lines(pri$x, pri$y, col = "red")
  legend("topright", legend = c("posterior", "prior"),
         lty = "solid", col = c("black", "red"))
}

split_params <- function(.data, param_col, ...) {
  separate(.data, !!param_col, c("biome", "pft", "variable"),
           sep = "\\.", fill = "left", ...)
}
