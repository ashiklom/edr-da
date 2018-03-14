summarize_results <- function(samples, prefix, burnin = 30000) {

  samples_burned <- getSample(samples, start = burnin, coda = TRUE)

  params <- coda::varnames(samples_burned)

  samples_summary_raw <- summary(samples_burned)

  samples_summary <- samples_summary_raw %>%
    .[c("quantiles", "statistics")] %>%
    map(~as_tibble(., rownames = "parameter")) %>%
    reduce(left_join) %>%
    separate(parameter, c("biome", "pft", "param"), sep = "\\.", fill = "left")

  # Summary of parameter estimates
  summary_plot <- samples_summary %>%
    filter(
      param != "residual"
    ) %>%
    mutate(
      pft = factor(pft, levels = unique(pft)),
      param = factor(param, levels = unique(param))
    ) %>%
    ggplot() +
    aes(x = pft, y = Mean, ymin = `2.5%`, ymax = `97.5%`, color = pft) +
    geom_pointrange() +
    facet_wrap(~ param, scales = "free") +
    theme(axis.text.x = element_blank())
  ggsave(paste(prefix, "summary", "pdf", sep = "."), summary_plot)

  # Density plot
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
