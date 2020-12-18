library(tidyverse)

l7dat <- readRDS("proc_results/landsat7_spectra.rds")

ndvi <- function(b3, b4) (b4 - b3) / (b4 + b3)

plot_dat <- l7dat %>%
  mutate(
    NDVI_mean = ndvi(`3_Mean`, `4_Mean`),
    NDVI_lo = ndvi(`3_Mean` - `3_SD`, `4_Mean` - `4_SD`),
    NDVI_hi = ndvi(`3_Mean` + `3_SD`, `4_Mean` + `4_SD`)
  )

l7long <- l7dat %>%
  gather(variable, value, -date, -site) %>%
  separate(variable, c("band", "stat")) %>%
  filter(band != "NDVI") %>%
  mutate(band = as.numeric(band)) %>%
  spread(stat, value)

l7long %>%
  filter(band %in% 2:4) %>%
  ggplot() +
  aes(
    x = date, y = Mean, ymin = Mean - SD, ymax = Mean + SD,
    color = site, fill = site
  ) +
  geom_ribbon() +
  geom_line() +
  facet_grid(band ~ site, scales = "free")

ggplot(l7long %>% filter(band %in% 2:4))

plt <- plot_dat %>%
  filter(date > "2008-06-15", site != "IDS36") %>%
  mutate_at(
    vars(matches("NDVI")),
    ~if_else(. > 1, NA_real_, .)
  ) %>%
  ggplot() +
  aes(x = date, y = NDVI_mean, ymin = NDVI_lo, ymax = NDVI_hi) +
  geom_ribbon(fill = "blue", alpha = 0.5) +
  geom_line() +
  facet_wrap(~ site) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(
    x = "Date", y = expression("NDVI" ~ ("Mean" %+-% "SD"))
  ) +
  theme_bw()
ggsave("figures/landsat_ts.pdf", plt, width = 7, height = 7)
