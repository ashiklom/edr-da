library(tidyverse)
library(redr)
library(PEcAn.ED2)
library(cowplot)
import::from(here, here)
import::from(ggrepel, geom_text_repel, geom_label_repel)
import::from(maps, geogmap = map)
import::from(broom, tidy)

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
    frac_evergreen = mean(
      grepl("conifer|pine", pft_name, ignore.case = TRUE)
    ),
    frac_evergreen_wtd = weighted.mean(
      grepl("conifer|pine", pft_name, ignore.case = TRUE),
      dbh
    )
  ) %>%
  mutate(site_label = str_extract(site_name, "^.*?_") %>% str_remove("_")) %>%
  ungroup()

select_labels <- readLines("other_site_data/selected_sites")
avgs_select <- avgs %>%
  mutate(
    selected = site_label %in% select_labels,
    fontface = if_else(selected, "bold", "plain")
  )
write_csv(avgs_select, "site_structure.csv")
plt <- ggplot(avgs_select) +
  aes(x = tot_dens * (1000), y = mean_dbh, color = frac_evergreen, label = site_label) +
  geom_point(aes(size = selected)) +
  geom_text_repel(
    aes(fontface = fontface),
    max.iter = 2000, min.segment.length = 1,
    point.padding = 0.5,
    data = filter(avgs_select, selected)
  ) +
  scale_color_continuous(high = "brown", low = "green4") +
  scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.2)) +
  labs(
    color = "Frac. evergreen",
    size = "Forward simulation",
    x = expression("Stem density" ~ ("trees"~"ha"^-1)),
    y = "Mean DBH (cm)"
  ) +
  theme(
    legend.position = c(0.98, 0.5),
    legend.justification = c(1, 1)
  )
ggsave("figures/selected_sites.pdf", plt)

mapdat <- geogmap("state", plot = FALSE, fill = TRUE) %>%
  tidy() %>%
  as_tibble()

avgs_select %>%
  arrange(selected) %>%
  ggplot() +
  aes(x = longitude, y = latitude, color = frac_evergreen) +
  geom_polygon(aes(x = long, y = lat, group = group), data = mapdat,
               color = "grey50", fill = NA) +
  geom_point() +
  geom_label_repel(
    aes(label = site_label),
    max.iter = 2000, min.segment.length = 0.2,
    point.padding = 0.5,
    label.padding = 0.1,
    data = filter(avgs_select, selected)
  ) +
  coord_map("sinusoidal", xlim = c(-100, -70), ylim = c(35, 50)) +
  scale_color_continuous(high = "brown", low = "green4") +
  scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.2)) +
  guides(color = FALSE) +
  theme(
    axis.title = element_blank()
  ) -> plt2
ggsave("figures/sitemap.pdf", plt2, width = 7, height = 5)

# Both plots together
plt_both <- ggdraw() +
  draw_plot(
    plt + theme(
      legend.position = c(1, 1),
      legend.text = element_text(size = rel(0.6)),
      legend.title = element_text(size = rel(0.8))
    ),
    0, 0, 1, 1
  ) +
  draw_plot(
    plt2 + theme(axis.text = element_text(size = rel(0.6))),
    -0.05, 0.3,
    scale = 0.55
  )
save_plot("figures/sites_both.pdf", plt_both, base_width = 7, base_height = 7)
