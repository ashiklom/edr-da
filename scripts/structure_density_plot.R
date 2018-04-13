library(tidyverse)
library(redr)
library(PEcAn.ED2)
import::from(here, here)
import::from(ggrepel, geom_text_repel)

site_dir <- here("sites")
site_list <- list.files(site_dir)

site_df <- tibble(
  site_name = site_list,
  site_path = file.path(site_dir, site_name)
)

site_df2 <- site_df %>%
  mutate(
    prefix = map_chr(site_path, get_site_prefix),
    veg = map(prefix, read_ed_veg)
  )

site_df3 <- site_df2 %>%
  mutate(
    css = map(veg, "css"),
    latitude = map_dbl(veg, "latitude"),
    longitude = map_dbl(veg, "longitude")
  ) %>%
  unnest(css) %>%
  left_join(ed_pft_mapping, by = c("pft" = "ed_pft_number"))

avgs <- site_df3 %>%
  filter(dbh >= 12.5) %>%
  group_by(site_name, latitude, longitude) %>%
  summarize(
    mean_dbh = mean(dbh, na.rm = TRUE),
    tot_dens = sum(n, na.rm = TRUE),
    frac_evergreen = weighted.mean(
      grepl("conifer|pine", pft_name, ignore.case = TRUE),
      dbh
    )
  ) %>%
  mutate(site_label = str_extract(site_name, "^.*?_") %>% str_remove("_")) %>%
  ungroup()

select_labels <- c("NC17", "PB11", "OF05", "SF01", "BI08", "BI02")

avgs_select <- avgs %>%
  mutate(selected = site_label %in% select_labels)

ggplot(avgs_select) +
  aes(x = tot_dens, y = mean_dbh, color = frac_evergreen, label = site_label) +
  geom_point(aes(size = selected)) +
  geom_text_repel(max.iter = 2000, min.segment.length = 0.2) +
  scale_color_continuous(high = "brown", low = "green4") +
  scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1))
