library(drake)

source("drake/packages.R")
source("drake/functions.R")
source("drake/globals.R")
source("drake/plan.R")

requireNamespace("future", quietly = TRUE)

future::plan("multiprocess")

drake_config(
  plan,
  parallelism = "future",
  jobs = 6
)
