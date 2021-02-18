library(conflicted)

library(fs)
library(here)
library(BayesianTools)
library(future)
library(doFuture)

stopifnot(
  requireNamespace("mvtnorm", quietly = TRUE)
)

options(future.rng.onMisuse = "ignore")

future::plan(future::multicore())
future::plan()
registerDoFuture()

argv <- commandArgs(trailingOnly = TRUE)
stopifnot(all(argv %in% "resume"))
resume <- "resume" %in% argv

pkgload::load_all()

# Load inversion data
load(here("data-raw", "inversion-data", "inversion-data.RData"))

# Inversion settings
niter <- 5000
max_iter <- 1e7
max_attempts <- floor(max_iter / niter)
attempt <- 0
threshold <- 1.2
target_neff <- 500
nchains <- 3

start_value <- prior$sampler(nchains)
start_value[, grep("b1Bl", colnames(start_value))] <- runif(nchains * 5, 0, 0.001)
start_value[, grep("b1Bw", colnames(start_value))] <- runif(nchains * 5, 0, 0.001)

message("Running on ", availableCores(), " cores.")
for (i in seq(1, nchains)) {
  message("Testing ", i)
  params <- start_value[i,]
  print(overall_likelihood(params))
}

param_names <- names(params)

# Create directory for storage
sampdir <- strftime(Sys.time(), "%Y-%m-%d-%H%M")
outtag <- "revision"
base_outdir <- file.path("multi_site_pda_results", outtag)
outdir <- file.path(base_outdir, sampdir)
dir_create(outdir)
message("Storing results in: ", outdir)
writeLines(param_names, file.path(outdir, "param_names.txt"))

# Create BayesianTools setup
message("Creating setup")
newsamples <- BayesianTools::createBayesianSetup(overall_likelihood, prior)
settings <- list(
  iterations = niter,
  consoleUpdates = 10,
  startValue = start_value
)

if (resume) {
  message("Resuming from previous sampling")
  last_samplefile <- tail(list.files(base_outdir, ".rds",
                          recursive = TRUE, full.names = TRUE), 1)
  message("Resuming from sample file: ", last_samplefile)
  stopifnot(length(last_samplefile) > 0, file.exists(last_samplefile))
  samples <- readRDS(last_samplefile)
} else {
  message("Starting fresh inversion")
  samples <- newsamples
}

repeat {
  attempt <- attempt + 1
  if (attempt > (max_attempts + 1)) {
    message("Failed to converge after attempt: ", attempt)
    break
  }
  message("Sampling attempt: ", attempt)
  samples <- BayesianTools::runMCMC(samples, settings = settings, sampler = "DEzs")
  saveRDS(samples, file.path(outdir, "current_samples.rds"))
  nsamp <- nrow(samples$chain[[1]])
  maxsamp <- 20000
  coda_samples <- BayesianTools::getSample(
    samples,
    start = if (nsamp > maxsamp) maxsamp else floor(nsamp / 2),
    thin = "auto",
    coda = TRUE
  )
  gd <- tryCatch(
    coda::gelman.diag(
      coda_samples,
      multivariate = TRUE,
      autoburnin = FALSE
    ),
    error = function(e) NULL)

  # Check PSRF can be calculated. If not, probably all the values are the same,
  # so no convergence.
  if (is.null(gd)) {
    message("GD calc failed. No convergence.")
    next
  }

  # Check PSRF values. If any parameters haven't converged, continue. Don't
  # worry about the multivariate PSRF -- only use the univariate ones.
  ## psrf <- c(gd[["psrf"]][, 1], gd[["mpsrf"]])
  ## names(psrf) <- c(param_names, "mpsrf")
  psrf <- gd[["psrf"]][, 1]
  names(psrf) <- param_names
  exceeds <- psrf > threshold
  if (any(exceeds)) {
    message(sum(exceeds), " parameters have not converged, including:")
    print(head(psrf[exceeds], 20))
    next
  }

  # Check for trend in the posterior samples. If there's a slope, no
  # convergence.
  has_slope <- function(y, threshold = 0.01) {
    x <- seq_along(y)
    fit <- lm(y ~ x)
    sfit <- summary(fit)
    pval <- sfit$coefficients["x", "Pr(>|t|)"]
    pval < threshold
  }
  is_sloped <- rowSums(do.call(cbind, lapply(coda_samples, function(x) apply(x, 2, has_slope)))) > 0
  if (any(is_sloped)) {
    message("Converged, but slope detected for ", sum(is_sloped), " parameters, including: ")
    print(head(param_names[is_sloped], 10))
  }

  # Finally, check effective sample size.
  too_few <- neff < target_neff
  if (any(too_few)) {
    message("Converged, but the following parameters have too few samples:")
    print(neff[too_few])
    next
  }

  message("Converged!")
}
