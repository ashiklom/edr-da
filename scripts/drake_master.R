#!/usr/bin/env Rscript
library(drake)
library(redr)
library(magrittr)
library(ggplot2)
library(PEcAnRTM)
pkgconfig::set_config("drake::strings_in_dots" = "literals")

expose_imports(redr)

nens <- 50
pda_dir <- here::here("multi_site_pda_results")
selected_sites <- readLines(here::here("other_site_data", "site_list"))

pre_plan <- drake_plan(
  other_posteriors = readRDS(
    file_in("ed-inputs/istem-posteriors/processed.rds")
  ),
  samplefile = target(
    command = tail(list.files(pda_dir, "\\.rds$",
                              recursive = TRUE, full.names = TRUE), 1),
    trigger = trigger(
      change = list.files(pda_dir, "\\.rds$", recursive = TRUE)
    )),
  samples_bt = readRDS(samplefile),
  param_names = readLines(file_in("param_names.txt")),
  params_matrix = BayesianTools::getSample(samples_bt) %>%
    `colnames<-`(param_names) %>%
    .[seq(floor(NROW(.) / 2), NROW(.)), ],
  ensemble_trait_list = preprocess_samples(
    samples_bt,
    param_names,
    other_posteriors,
    50,
    fix_allom2 = TRUE
  ),
  tidy_params = params_matrix %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sample = dplyr::row_number()) %>%
    tidyr::gather(parameter, value, -sample) %>%
    dplyr::mutate(
      pft = gsub("^.*\\.(.*)\\..*$", "\\1", parameter),
      pft = forcats::fct_inorder(pft),
      ipft = as.integer(pft),
      parameter = gsub("^.*\\..*\\.(.*)$", "\\1", parameter),
      parameter = forcats::fct_inorder(parameter)
    ),
  pft_posteriors = tidy_params %>%
    dplyr::filter(!grepl("soil|residual", parameter)) %>%
    ggplot() +
    aes(x = pft, y = value, color = pft, fill = pft) +
    geom_violin() +
    facet_wrap(vars(parameter), scales = "free_y") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()),
  site_files <- purrr::map_chr(
    selected_sites,
    ~head(list.files(file.path("sites", .x), "css$", full.names = TRUE), 1)
  ) %>%
    setNames(selected_sites),
  site_df = purrr::map_dfr(site_files, PEcAn.ED2::read_css, .id = "site") %>%
    tibble::as_tibble() %>%
    dplyr::select(site, year = time, dbh, ipft = pft, nplant = n) %>%
    dplyr::mutate(
      ipft = get_ipft(ipft),
      hite = dbh2h(dbh, ipft)
    ) %>%
    dplyr::group_by(site, year, ipft, hite) %>%
    dplyr::summarize(
      dbh = mean(dbh),
      nplant = sum(nplant)
    ) %>%
    dplyr::arrange(desc(hite)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(site, year) %>%
    dplyr::mutate(cohort = dplyr::row_number()) %>%
    dplyr::ungroup(),
  site_lai_samples = site_df %>%
    dplyr::left_join(
      tidy_params %>%
        dplyr::filter(parameter %in% c("b1Bl", "SLA", "clumping_factor")) %>%
        tidyr::spread(parameter, value) %>%
        dplyr::select(sample, ipft, SLA, b1Bl, clumping_factor)
    ) %>%
    dplyr::mutate(
      bleaf = size2bl(dbh, b1Bl, purrr::map_dbl(allom_mu, "b2Bl")[ipft]),
      lai = nplant * bleaf * SLA,
      elai = lai * clumping_factor
    ),
  site_lai_summary = site_lai_samples %>%
    dplyr::group_by(site, year, ipft, hite, dbh, nplant, cohort) %>%
    dplyr::summarize_at(
      c("lai", "elai"),
      dplyr::funs(
        mean = mean(.),
        sd = sd(.),
        q025 = quantile(., 0.025),
        q975 = quantile(., 0.975)
      )
    )
)

spec_validation_template <- drake_plan(
  sitespec_predicted = predict_site_spectra(
    params_matrix,
    "site__",
    nsamp = 1000,
    dedup = TRUE,
    progress = FALSE
  ),
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
      geom_ribbon(aes(ymin = pmax(albedo_r_mean - albedo_r_sd, 0), ymax = albedo_r_mean + albedo_r_sd,
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
  ),
  sitespec_valid_plot_png = ggsave("figures/validplot_site__.png", sitespec_valid_plot_site__)
)

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
make(current_plan, jobs = 8)

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
make(plan, jobs = parallel::detectCores())
