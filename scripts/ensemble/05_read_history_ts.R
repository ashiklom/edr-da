library(redr)
library(tidyverse)
import::from(here, inhere = here)
import::from(progress, progress_bar)
import::from(fs, dir_ls)
import::from(hdf5r, H5File)
import::from(lubridate, as_date, month)

var_table <- read_csv("inst/ed_state_vars.csv")

var_table %>%
  filter(type == "cohort") %>%
  select(variable, group_type) %>%
  print(n = Inf)

run_dir <- inhere("ed-outputs", "EDR_sim_output_BH02_2008", "outputs")
histfiles <- list.files(run_dir, "history-S", full.names = TRUE)

read_hist_variable <- function(histfile, variable) {
  hf <- hdf5r::H5File$new(histfile, "r")
  hf[[variable]]$read()
}

base_tibble <- tibble(files = histfiles) %>%
  mutate(
    datetime = str_extract(files, "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}") %>%
      as_date()
  )

bstore <- base_tibble %>%
  mutate(
    pft = map(files, read_hist_variable, "PFT"),
    storage = map(files, read_hist_variable, "BSTORAGE"),
    height = map(files, read_hist_variable, "HITE"),
    data = pmap(list(PFT = pft, storage = storage, height = height), tibble)
  ) %>%
  unnest(data) %>%
  select(-files) %>%
  group_by(datetime) %>%
  arrange(desc(height)) %>%
  mutate(cohort_id = row_number()) %>%
  ungroup() %>%
  arrange(datetime, cohort_id)

ggplot(bstore) +
  aes(x = datetime, y = storage, color = factor(PFT), group = cohort_id) +
  geom_line()

#cb_list <- base_tibble %>%
  #mutate(
    #month = month(datetime),
    #cb_full = map(files, read_hist_variable, "CB"),
    #rn = map(cb_full, ~seq_len(nrow(.))),
    #cb_full = map2(cb_full, rn, `rownames<-`),
    #cb_full = map(cb_full, as_tibble, rownames = "month_id")
  #) %>%
  #unnest(cb_full)

#cb_list %>%
  #gather(cohort, value, -(files:month_id)) %>%
  #filter(month_id == month) %>%
  #ggplot() +
  #aes(x = datetime, y = value, color = cohort) +
  #geom_line()

#cb_list %>% unnest(cb_full)

  #map(histfiles, read_hist_variable, "CB") %>%
  #map(~.[13, ]) %>%
  #reduce(rbind) %>%
  #as_tibble() %>%
  #gather("cohort", "CB", -datetime)

#ggplot(cb_list) +
  #aes(x = datetime, y = CB, color = cohort) +
  #geom_line()


#hist_raw <- map(histfiles[1:10], safely(read_ed_history)) %>% transpose()
#hist_results <- hist_raw$result

#cohort_df <- tibble(
  #datetime = map(hist_raw$result, "datetime") %>% reduce(c),
  #cohort = map(hist_raw$result, "cohort")
#) %>%
  #unnest()

## Pull out carbon balance

#ggplot(cohort_df) +
  #aes(x = )
