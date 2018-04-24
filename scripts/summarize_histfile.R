library(redr)
library(tidyverse)
import::from(here, inhere = here)

run_dir <- inhere("ed-outputs", "EDR_sim_output_BH02_2008", "outputs")
histfiles <- list.files(run_dir, "history-S", full.names = TRUE)
histfile <- histfiles[1]

l <- read_ed_history(histfile)

scalar_vars <- tibble(
  variable = l$scalar$variable,
  type = "scalar",
  size_per_group = 1
)

cohort_vars <- tibble(
  variable = colnames(l$cohort),
  type = "cohort",
  size_per_group = map(l$cohort, 1) %>% map_int(length)
) %>%
  mutate(
    group_type = recode(
      size_per_group,
      `1` = "single cohort",
      `5` = "disturbance type1",
      `10` = "radiation profile",
      `13` = "monthly+"
    )
  ) %>%
  select(-size_per_group) %>%
  filter(variable != "cohort_id")

pft_vars <- tibble(
  variable = colnames(l$pft),
  type = "PFT",
  size_per_group = map(l$pft, 1) %>% map_int(length)
) %>%
  mutate(
    group_type = recode(
      size_per_group,
      `1` = "single PFT",
      `8` = "height class",
      `11` = "DBH class"
    )
  ) %>%
  select(-size_per_group) %>%
  filter(variable != "pft_id")

soil_vars <- tibble(
  variable = colnames(l$soil),
  type = "soil",
  group_type = "soil depth"
)

other_vars <- tribble(
  ~variable, ~type, ~group_type,
  "AVG_MONTHLY_PCPG", "other", "monthly",
  "DISTURBANCE_MEMORY", "other", "disturbance type2",
  "DISTURBANCE_RATES", "other", "disturbance type2",
  "HGT_CLASS", "other", "height class",
  "LAMBDA_FIRE", "other", "monthly",
  "WORKLOAD", "other", "monthly+"
)

allvars <- bind_rows(
  scalar_vars,
  cohort_vars,
  pft_vars,
  soil_vars,
  other_vars
)

write_csv(allvars, "inst/ed_state_vars.csv")
