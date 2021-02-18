library(drake)

run_config <- Sys.getenv("RUN_CONFIG")
if (run_config == "") run_config <- "revision-fixed"

message("Processing ", run_config)

outpath <- here::here("multi_site_pda_results", run_config)
stopifnot(file.exists(outpath))

cachepath <- file.path(outpath, ".drake")
if (!file.exists(cachepath)) {
  dc <- new_cache(cachepath)
} else {
  dc <- drake_cache(cachepath)
}

source("drake/packages.R")
source("drake/functions.R")
source("drake/globals.R")
source("drake/plan.R")

requireNamespace("future", quietly = TRUE)

future::plan("multicore")

drake_config(
  plan,
  parallelism = "future",
  jobs = 6,
  cache = dc
)
