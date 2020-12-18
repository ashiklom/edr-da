library(dplyr)
library(ggplot2)
library(ggmap)
library(sf)

usa <- rnaturalearth::ne_states(
  country = c("United States of America", "Canada"),
  returnclass = "sf"
)

drake::loadd(site_structure_data)
site_structure_sf <- st_as_sf(
  site_structure_data,
  coords = c("longitude", "latitude"),
  crs = 4326
)

## bbox <- st_bbox(site_structure_sf)
## names(bbox) <- c("left", "bottom", "right", "top")

bbox <- c(left = -95, bottom = 42, right = -84, top = 48)
basemap <- get_stamenmap(bbox, maptype = "terrain-background", zoom = 8,
                         where = "~/.local/share/stamenmaps") %>%
  ggmap()

best_sites <- c("BH07", "PB09")

sitemap <- basemap +
  geom_sf(data = site_structure_sf, inherit.aes = FALSE,
          aes(color = site_name %in% best_sites)) +
  guides(color = FALSE) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  theme(axis.title = element_blank())
if (interactive()) sitemap
ggsave("scripts/figures/hl-sitemap.png", sitemap,
       width = 5, height = 4.5, units = "in", dpi = 300)


# Spectra
drake::loadd(observed_spectra)
drake::loadd(predicted_spectra)

for (i in seq_along(best_sites)) {
  s1 <- best_sites[i]
  obs1 <- observed_spectra %>%
    filter(site == s1)
  pred1 <- predicted_spectra %>%
    filter(site == s1)
  p1 <- ggplot() +
    aes(x = wavelength) +
    geom_ribbon(aes(ymin = albedo_q025, ymax = albedo_q975, fill = "ed"),
                alpha = 0.5, data = pred1) +
    geom_line(aes(y = albedo_mean, color = "ed"), data = pred1,
              size = 2) +
    geom_line(aes(y = observed, group = iobs, color = "aviris"), data = obs1) +
    scale_color_manual(values = c(aviris = "black", ed = "green4")) +
    scale_fill_manual(values = c(ed = "green1")) +
    labs(x = "Wavelength (nm)", y = "Reflectance") +
    ylim(0, 0.55) +
    guides(color = FALSE, fill = FALSE) +
    theme_bw()
  if (interactive()) p1
  ggsave(glue::glue("scripts/figures/hl-p{i}.png"), p1,
         width = 3.6, height = 2.9, units = "in", dpi = 300)
}

# Rayrend
library(rayrender)
library(lhs)

set.seed(999)


f1 <- "sites/BH07/FFT.2008.lat43.4199lon-89.8011.css"
d1 <- read.table(f1, header = TRUE) %>%
  as_tibble() %>%
  mutate(
    hite = redr::dbh2h(dbh, redr::get_ipft(pft)),
    lai = n * 15 * redr::size2bl(dbh, exp(-3.7), 1.8)
  )
xz1 <- (randomLHS(nrow(d1), 2) - 0.5) * 4

f2 <- "sites/PB09/FFT.2008.lat46.4325lon-91.5199.css"
d2 <- read.table(f2, header = TRUE) %>%
  as_tibble() %>%
  mutate(
    hite = redr::dbh2h(dbh, redr::get_ipft(pft)),
    lai = n * 15 * redr::size2bl(dbh, exp(-3.7), 1.8)
  )
xz2 <- (randomLHS(nrow(d2), 2) - 0.5) * 4

tree <- function(hite, rad, x = 0, z = 0, pft = 6) {
  color <- as.character(factor(pft, c(6, 8:11),
                               c("orange", "purple", "forestgreen", "blue4", "red3")))
  center_hite <- hite - rad
  trunkwid <- 0.02 * d2$hite
  trunk <- cube(x = x, y = center_hite / 2, z = z,
                xwidth = trunkwid,
                ywidth = center_hite,
                zwidth = trunkwid,
                material = diffuse(color = "brown"))
  if (pft %in% c(6, 8)) {
    # Conifer -- use a cone
    canopy <- cone(
      start = c(x, 0.2 * hite, z),
      end = c(x, hite, z),
      radius = rad,
      material = diffuse(color = color)
    )
  } else {
    # Deciduous
    canopy <- sphere(x = x, y = center_hite, z = z,
                     radius = rad, scale = c(0.7, 1, 0.7),
                     material = diffuse(color = color))
  }
  rbind(trunk, canopy)
}
trees1 <- Map(tree, hite = d1$hite/10, rad = sqrt(d1$lai),
              x = xz1[,1], z = xz1[,2], pft = d1$pft)
trees2 <- Map(tree, hite = d2$hite/10, rad = sqrt(d2$lai),
              x = xz2[,1], z = xz2[,2], pft = d2$pft)

png("scripts/figures/hl-sim1.png", width = 2.8, height = 2.8,
    units = "in", res = 300)
generate_ground(depth = 0, material = diffuse(color = "gray99")) %>%
  add_object(do.call(rbind, trees1)) %>%
  render_scene(lookfrom = c(0, 10, 15), samples = 400)
dev.off()

png("scripts/figures/hl-sim2.png", width = 2.8, height = 2.8,
    units = "in", res = 300)
generate_ground(depth = 0, material = diffuse(color = "gray99")) %>%
  add_object(do.call(rbind, trees2)) %>%
  render_scene(lookfrom = c(0, 10, 15), samples = 400)
dev.off()
