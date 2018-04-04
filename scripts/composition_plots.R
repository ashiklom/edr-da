library(redr)
library(here)
library(tidyverse)
library(PEcAn.ED2)
library(ggsci)
import::from("RColorBrewer", "brewer.pal")

sitelist <- readLines(here("other_site_data", "site_list"))

site <- sitelist[1]

pft_names <- levels(ed_pft_mapping$pft_name)
pft_colors <- brewer.pal(length(pft_names), "Set1")
names(pft_colors) <- pft_names

get_css <- function(site) {
  site_path <- file.path("sites", site)
  cssfile <- list.files(site_path, "\\.css$", full.names = TRUE)[1]
  css <- read_delim(cssfile, " ") %>%
    dplyr::left_join(ed_pft_mapping, by = c("pft" = "ed_pft_number")) %>%
    dplyr::mutate(site = !!site)
  css
}

all_css <- map_dfr(sitelist, get_css)
ylims <- c(0, max(all_css$dbh))

composition_plot <- function(site, all_css) {
  css <- dplyr::filter(all_css, site == !!site)
  plt <- ggplot(css) +
    aes(x = pft_name, y = dbh, color = pft_name) +
    geom_jitter() +
    scale_color_manual(values = pft_colors) +
    ggtitle(site) +
    coord_cartesian(ylim = ylims)
  print(plt)
}

pdf("other_site_data/composition_plots.pdf")
walk(sitelist, composition_plot, all_css = all_css)
dev.off()
