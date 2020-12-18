library(tidyverse)
import::from(here, inhere = here)

summary_dat <- readRDS(inhere("sync_data", "msp20180402", "pda_summary.rds"))

vars_rename <- tribble(
  ~variable, ~level,
  "prospect_N", expression("PROSPECT" ~ "No. meso"),
  "prospect_Cab", expression("PROSPECT Chl." ~ (mu * g ~ cm^-2)),
  "prospect_Car", expression("PROSPECT Car." ~ (mu * g ~ cm^-2)),
  "prospect_Cw", expression("PROSPECT Water" ~ (g ~ cm^-2)),
  "prospect_Cm", expression("PROSPECT LMA" ~ (g ~ cm^2)),
  "SLA", expression("Specific leaf area" ~ (mm^2 ~ g^-1)),
  "b1Bl", expression("Leaf biomass allom." ~ "base"),
  "b2Bl", expression("Leaf biomass allom." ~ "exp."),
  "clumping_factor", expression("Clumping factor" ~ (0 - 1)),
  "orient_factor", expression("Orientation factor" ~ (-1 - 1))
)

vars_rename_vec <- vars_rename %>%
  select(2, 1) %>%
  deframe()

max_b1bl <- 5
summary_plot <- summary_dat %>%
  mutate(
    Mean = if_else(variable == "b1Bl" & Mean > max_b1bl, NA_real_, Mean),
    `2.5%` = if_else(variable == "b1Bl" & Mean > max_b1bl, NA_real_, `2.5%`),
    `97.5%` = if_else(variable == "b1Bl" & Mean > max_b1bl, NA_real_, `97.5%`),
    variable = factor(variable) %>% fct_recode(!!!vars_rename_vec),
  )
dodge <- position_dodge(width = 0.5)
plt <- ggplot(summary_plot) +
  aes(
    x = pft, y = Mean, ymin = `2.5%`, ymax = `97.5%`,
    linetype = type, color = pft, group = pft_type, shape = type
  ) +
  geom_linerange(position = dodge, size = 0.5) +
  geom_point(position = dodge, size = 2) +
  facet_wrap(~ variable, scales = "free", labeller = "label_parsed") +
  scale_linetype_manual(values = c("prior" = "dashed", "posterior" = "solid")) +
  labs(
    x = "Plant functional type",
    y = "Mean and 95% CI",
    color = "",
    linetype = "",
    shape = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
  )
ggsave("figures/pda_summary.pdf", plt, width = 9, height = 6)
