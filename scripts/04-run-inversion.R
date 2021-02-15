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
  # Need this to reset parallelism
  samples$setup$likelihood <- newsamples$likelihood
  samples$setup$posterior <- newsamples$posterior
} else {
  message("Starting fresh inversion")
  samples <- newsamples
}

repeat {
  attempt <- attempt + 1
  message("Sampling attempt: ", attempt)
  samples <- BayesianTools::runMCMC(samples, settings = settings, sampler = "DEzs")
  saveRDS(samples, file.path(outdir, "current_samples.rds"))
  nsamp <- nrow(samples$chain[[1]])
  coda_samples <- BayesianTools::getSample(
    samples,
    start = if (nsamp > 5000) 5000 else floor(nsamp / 2),
    thin = "auto",
    coda = TRUE
  )
  gd <- tryCatch(
    coda::gelman.diag(
      coda_samples,
      multivariate = TRUE,
      autoburnin = FALSE
    ),
    error = function (e) NULL)
  if (is.null(gd)) {
    message("GD calc failed. No convergence.")
  } else {
    psrf <- c(gd[["psrf"]][, 1], gd[["mpsrf"]])
    names(psrf) <- c(param_names, "mpsrf")
    exceeds <- psrf > threshold
    if (any(exceeds)) {
      message("The following parameters have not converged:")
      print(psrf[exceeds])
    } else {
      neff <- coda::effectiveSize(coda_samples)
      too_few <- neff < target_neff
      if (!any(too_few)) {
        message("Converged!")
        break
      }
      message("Converged, but the following parameters have too few samples:")
      print(neff[too_few])
    }
  }

  message("Resuming sampling...")
  if (attempt > max_attempts) {
    message("Failed to converge after attempt: ", attempt)
    break
  }
}
