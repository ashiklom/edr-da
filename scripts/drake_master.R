library(drake)
import::from(magrittr, "%>%")
import::from(tibble, as_tibble)
import::from(dplyr, mutate)
pkgconfig::set_config("drake::strings_in_dots" = "literals")

## source(here::here("scripts/drake_functions"))

nens <- 50

pre_plan <- drake_plan(
  other_posteriors = readRDS(file_in("ed-inputs/istem-posteriors/processed.rds")),
  samplefile = file_in("multi_site_pda_results/testsamples/current_samples.rds"),
  param_names = readLines(file_in("param_names.txt")),
  ensemble_trait_list = preprocess_samples(samplefile, param_names, other_posteriors, nens)
)

ed_plan_template <- drake_plan(
  edresult = run_ed_site_ens(
    site = "site__",
    trait_values = ensemble_trait_list[[ens__]],
    outdir_prefix = file.path("testsamples", ens__),
    end_date = "2011-12-31"
  )
)

selected_sites <- c(
    ## "--site=OF05_site_1-25710"   # Initial
    #"SF03_site_1-25721" # Skipped for now
    "IDS36_site_1-25686",
    "BH07_site_1-25669",
    "AK60_site_1-25674",
    "OF02_site_1-25708",
    "BH05_site_1-25667"
)

ed_plan <- evaluate_plan(
  ed_plan_template,
  rules = list(
    site__ = selected_sites,
    ens__ = seq_len(nens)
  )
)

plan <- bind_plans(pre_plan, ed_plan)
dconf <- drake_config(plan)
make(plan)
## make(plan, parallelism = "parLapply", jobs = 2)
