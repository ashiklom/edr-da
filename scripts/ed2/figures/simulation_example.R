library(ggplot2)
library(purrr)
library(dplyr)
library(forcats)
library(redr)
library(PEcAnRTM)
library(PEcAn.ED2)
import::from(lubridate, as_date, year, month, mday)
import::from(progress, progress_bar)
import::from(imguR, imgur, imgur_off)
import::from(tidyr, unnest, spread)


ens_dir <- "ensemble_outputs/msp_hf20180402"
date <- "2009-07-02"
ens <- 1

run_edr_site <- function(date, site, ens_dir, ens = 1, pb = NULL) {
  on.exit(if (!is.null(pb)) pb$tick())
  site_dir <- list.files(ens_dir, site, full.names = TRUE)
  stopifnot(length(site_dir) == 1)
  ens_dir <- file.path(site_dir, sprintf("ens_%03d", ens))
  stopifnot(file.exists(ens_dir))
  ens_out_dir <- file.path(ens_dir, "out")
  stopifnot(file.exists(ens_out_dir))
  history_file <- list.files(ens_out_dir, date, full.names = TRUE)
  stopifnot(length(history_file) == 1)
  ed2in <- read_ed2in(file.path(ens_dir, "ED2IN"))
  ed2in$RK4_TOLERANCE <- 1e-5
  trait_values <- readRDS(file.path(ens_dir, "trait_values.rds"))
  run_edr_date(date, ed2in, trait_values)
}

sites <- readLines("other_site_data/selected_sites")
pb <- progress_bar$new(total = length(sites))
site_out <- map(
  sites, safely(run_edr_site),
  date = date, ens = ens, ens_dir = ens_dir, pb = pb
)

s2 <- transpose(site_out)
spec <- s2$result %>% setNames(sites) %>% discard(is.null)
tidyspec <- spec %>%
  imap(~tibble(site = .y, waves = 400:2500, refl = .x)) %>%
  bind_rows()

i1 <- imgur("png", width = 5, height = 5, units = "in", res = 300)
ggplot(tidyspec) +
  aes(x = waves, y = refl, color = fct_reorder(site, refl, max, .desc = TRUE)) +
  geom_line() +
  labs(color = "Site code", x = "Wavelength", y = "Reflectance") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95), legend.justification = c(1, 1))
i2 <- imgur_off(i1)
i2$link

data("sensor.rsr")
lsat <- tidyspec %>%
  group_by(site) %>%
  summarize(data = list(spec2landsat(refl))) %>%
  unnest()

l7 <- lsat %>%
  filter(landsat == "landsat7")

i1 <- imguR::imgur("png", width = 5, height = 5, units = "in", res = 300)
ggplot(l7) +
  aes(x = wavelength, y = value, color = fct_reorder(site, value, max, .desc = TRUE)) +
  geom_line() +
  labs(color = "Site code", x = "Wavelength", y = "Reflectance") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95), legend.justification = c(1, 1))
i2 <- imguR::imgur_off(i1)
i2$link
