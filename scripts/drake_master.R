#!/usr/bin/env Rscript
library(drake)
library(redr)
library(magrittr)
library(ggplot2)
pkgconfig::set_config("drake::strings_in_dots" = "literals")

nens <- 50

pre_plan <- drake_plan(
  other_posteriors = readRDS(file_in("ed-inputs/istem-posteriors/processed.rds")),
  samplefile = tail(
    list.files("multi_site_pda_results", "\\.rds$",
               recursive = TRUE, full.names = TRUE),
    1
  ),
  samples_bt = readRDS(samplefile),
  param_names = readLines(file_in("param_names.txt")),
  ensemble_trait_list = preprocess_samples(samples_bt, param_names,
                                                 other_posteriors, 50,
                                           fix_allom2 = TRUE),
  params_matrix = BayesianTools::getSample(samples_bt, start = 6000) %>%
    `colnames<-`(param_names)
)
make(pre_plan)

spec_validation_template <- drake_plan(
  sitespec_predicted = predict_site_spectra(params_matrix, "site__",
                                            nsamp = 1000, dedup = TRUE, progress = FALSE),
  sitespec_observed = load_observations("site__") %>%
    `colnames<-`(., as.character(seq_len(NCOL(.)))) %>%
    as.data.frame(., row.names = PEcAnRTM::wavelengths(.)) %>%
    tibble::as_tibble() %>%
    tibble::rownames_to_column("wavelength") %>%
    dplyr::mutate(wavelength = as.numeric(wavelength)) %>%
    tidyr::gather(iobs, observed, -wavelength)
)

do_valid_plot <- function(pred, obs) {
  suppressWarnings({
    ggplot() +
      aes(x = wavelength) +
      ## geom_ribbon(aes(ymin = albedo_q025, ymax = albedo_q975,
      ##                 fill = "Predicted", color = "Predicted"),
      ##             data = pred) +
      geom_ribbon(aes(ymin = albedo_r_mean - albedo_r_sd, ymax = albedo_r_mean + albedo_r_sd,
                      fill = "Predicted", color = "Predicted"),
                  data = pred) +
      ## geom_line(aes(y = albedo_r_mean, fill = "Predicted", color = "Predicted"),
      ##           data = pred) +
      geom_line(aes(y = observed, group = iobs,
                    fill = "Observed", color = "Observed"),
                data = obs) +
      scale_fill_manual(name = "",
                        values = c(Observed = "white", Predicted = "green1")) +
      scale_color_manual(name = "",
                         values = c(Observed = "grey25", Predicted = "green4")) +
      labs(x = "Wavelength (nm)", y = "Surface reflectance") +
      theme_bw()
  })
}

spec_validation_plot_template <- drake_plan(
  sitespec_valid_plot = do_valid_plot(
    sitespec_predicted_site__,
    sitespec_observed_site__
  )
)

selected_sites <- readLines("other_site_data/site_list")

spec_validation_plan <- evaluate_plan(
  spec_validation_template,
  rules = list(site__ = selected_sites)
)

spec_validation_plot_plan <- evaluate_plan(
  spec_validation_plot_template,
  rules = list(site__ = make.names(selected_sites))
)

current_plan <- bind_plans(pre_plan, spec_validation_plan, spec_validation_plot_plan)
current_config <- drake_config(current_plan)

lai_plan <- drake_plan(
  site_lai = purrr::map_dfr(selected_sites, calc_site_lai, param_matrix = param_matrix)
)

current_plan <- bind_plans(pre_plan, spec_validation_plan, spec_validation_plot_plan,
                           lai_plan)
current_config <- drake_config(current_plan)

make(current_plan, parallelism = "parLapply", jobs = 8)

ed_plan_template <- drake_plan(
  edresult = run_ed_site_ens(
    site = "site__",
    trait_values = ensemble_trait_list[[ens__]],
    outdir_prefix = file.path("testsamples", ens__),
    end_date = "2011-12-31"
  )
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
make(plan, parallelism = "parLapply", jobs = parallel::detectCores())
