import::from("BayesianTools", "getSample")
import::from("here", "here")
import::from("coda", "densplot", "traceplot")
import::from("purrr", "walk")

load("priors/mvtraits_priors.RData")
pfts <- rownames(means)

results_path <- here("pda_progress")

plots_path <- here("pda_plots")
dir.create(plots_path, showWarnings = FALSE)

flist <- dir(results_path)

plot_param <- function(param, out_coda, ...) {
  if (param == "residual") {
    message("Param is residual. Returning NULL.")
    return(NULL)
  }
  pft_regex <- paste0("(", paste(pfts, collapse = "|"), ")")
  out_pft <- gsub(paste0(pft_regex, "\\..*"), "\\1", param)
  out_param <- gsub(paste0(pft_regex, "\\."), "", param)
  par(mfrow = c(1, 2))
  traceplot(out_coda[, param], main = param)
  out_d <- density(as.matrix(out_coda[, param]))
  if (out_param %in% colnames(means)) {
    prior_mu <- means[out_pft, out_param]
    prior_sd <- sqrt(covars[out_param, out_param, out_pft])
    prior_samps <- rnorm(5000, prior_mu, prior_sd)
  } else if (out_param == "clumping_factor") {
    prior_samps <- runif(5000, 0, 1)
  } else if (out_param == "orient_factor") {
    prior_samps <- runif(5000, -0.5, 0.5)
  }
  prior_d <- density(prior_samps)
  xlim <- range(out_d$x, prior_d$x)
  ylim <- range(out_d$y, prior_d$y)
  plot(0, 0, type = 'n', xlim = xlim, ylim = ylim,
       xlab = out_param, ylab = "density", ...)
  lines(out_d$x, out_d$y, col = "black")
  lines(prior_d$x, prior_d$y, col = "red")
  legend("topright", c("PDA", "Prior"), col = c("black", "red"), lty = "solid")
}

analyze_file <- function(f, saveplot = TRUE) {
  out <- readRDS(file.path(results_path, f))
  out_coda <- getSample(out, coda = TRUE)
  param_names <- colnames(out_coda[[1]])

  pdf_fname <- gsub("rds$", "pdf", f)
  if (saveplot) pdf(file.path(plots_path, pdf_fname))
  walk(param_names, plot_param, out_coda = out_coda, main = f)
  if (saveplot) dev.off()
}

walk(flist, analyze_file, saveplot = TRUE)
