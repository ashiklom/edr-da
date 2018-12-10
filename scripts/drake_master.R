#!/usr/bin/env Rscript
library(drake)
library(redr)
library(magrittr)
library(ggplot2)
library(PEcAnRTM)
pkgconfig::set_config("drake::strings_in_dots" = "literals")

expose_imports(redr)
.datatable.aware <- TRUE

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
  site_files = selected_sites %>%
    purrr::map_chr(
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
    ## dplyr::group_by(site, year, ipft, hite) %>%
    ## dplyr::summarize(
    ##   dbh = mean(dbh),
    ##   nplant = mean(nplant)
    ## ) %>%
    ## dplyr::ungroup() %>%
    dplyr::group_by(site, year) %>%
    dplyr::arrange(desc(hite)) %>%
    dplyr::mutate(cohort = dplyr::row_number()) %>%
    dplyr::ungroup(),
  site_dt = data.table::as.data.table(site_df),
  tidy_params_dt = data.table::as.data.table(tidy_params),
  b2Bl = purrr::map_dbl(allom_mu, "b2Bl"),
  params_structure = data.table::dcast(
    tidy_params_dt[parameter %in% c("b1Bl", "SLA", "clumping_factor")],
    sample + pft + ipft ~ parameter,
    value.var = "value"
  )[
  , lapply(.SD, unname)
  ][
  , b2Bl := b2Bl[ipft]
  ][
    sample %in% sample.int(max(sample), 10000)
  ],
  site_lai_samples = site_dt[
    params_structure,
    on = "ipft",
    allow.cartesian = TRUE
  ][
  , bleaf := size2bl(dbh, b1Bl, b2Bl)
  ][
  , lai := nplant * bleaf * SLA
  ][
  , elai := lai * clumping_factor
  ],
  site_lai_summary_dt = site_lai_samples[
   ,
     list(lai_mean = mean(lai),
          lai_sd = sd(lai),
          lai_lo = quantile(lai, 0.025),
          lai_hi = quantile(lai, 0.975),
          elai_mean = mean(elai),
          elai_sd = sd(elai),
          elai_lo = quantile(elai, 0.025),
          elai_hi = quantile(elai, 0.975)),
     c("site", "year", "pft", "ipft", "hite", "dbh", "nplant", "cohort")
  ],
  site_lai_summary = tibble::as_tibble(site_lai_summary_dt) %>%
    purrr::modify(unname) %>%
    dplyr::group_by(site, year) %>%
    dplyr::arrange(hite) %>%
    dplyr::mutate(
      cum_lai = cumsum(lai_mean),
      cum_elai = cumsum(elai_mean)
    ) %>%
    dplyr::ungroup(),
  lai_observed = readr::read_csv(file_in("other_site_data/NASA_FFT_LAI_FPAR_Data.csv")) %>%
    dplyr::group_by(site = Site_Plot) %>%
    dplyr::summarize(
      obs_LAI = mean(LAI_CLX_10_60_m2_m2, na.rm = TRUE),
      obs_LAI_SD = sd(LAI_CLX_10_60_m2_m2, na.rm = TRUE),
      obs_LAI_lo = obs_LAI - obs_LAI_SD,
      obs_LAI_hi = obs_LAI + obs_LAI_SD
    ) %>%
    dplyr::filter_if(is.numeric, dplyr::all_vars(. > 0)),
  prior = create_prior(nsite = 1, heteroskedastic = FALSE, limits = TRUE),
  prior_samples = prior$sampler(1000),
  tidy_prior = tibble::as_tibble(prior_samples) %>%
    tidyr::gather(pft_param, value) %>%
    dplyr::mutate(
      pft = gsub("^.*\\.(.*)\\..*$", "\\1", pft_param),
      pft = forcats::fct_inorder(pft),
      ipft = as.integer(pft),
      parameter = gsub("^.*\\..*\\.(.*)$", "\\1", pft_param),
      parameter = forcats::fct_inorder(parameter),
      type = "prior"
    ),
  params_prior_posterior = tidy_prior %>%
    dplyr::mutate(type = "prior") %>%
    dplyr::bind_rows(tidy_params %>% dplyr::mutate(type = "posterior")) %>%
    dplyr::mutate(type = forcats::fct_inorder(type)),
  params_prior_posterior_plot = params_prior_posterior %>%
    dplyr::filter(!grepl("soil|residual", parameter)) %>%
    ggplot() +
    aes(x = interaction(type, pft), y = value, fill = pft, color = pft) +
    geom_violin() +
    facet_wrap(vars(parameter), scales = "free") +
    labs(x = "Prior vs. posterior", y = "Value") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
  params_prior_posterior_png = ggsave(
    file_out("figures/params_prior_posterior.png"),
    params_prior_posterior_plot,
    width = 8, height = 6, units = "in"
  )
)

site_lai_plots_template <- drake_plan(
  siteplot = site_lai_summary %>%
    dplyr::filter(site == "site__") %>%
    ggplot() +
    aes(y = hite, color = pft) +
    geom_segment(aes(x = 0, xend = xvar__, yend = hite)) +
    geom_point(aes(x = xvar__, y = hite)) +
    geom_vline(aes(xintercept = obs_LAI), linetype = "solid",
               data = lai_observed %>% dplyr::filter(site == "site__")) +
    geom_vline(aes(xintercept = obs_LAI_lo), linetype = "dashed",
               data = lai_observed %>% dplyr::filter(site == "site__")) +
    geom_vline(aes(xintercept = obs_LAI_hi), linetype = "dashed",
               data = lai_observed %>% dplyr::filter(site == "site__")) +
    labs(x = "Leaf area index", y = "Height (m)", color = "PFT") +
    ggtitle("site__"),
  siteplot_png = ggsave(
    "figures/siteplot_site___xvar__.png",
    siteplot_site___xvar__,
    width = 7, height = 7, units = "in"
  )
)

site_lai_plots_plan <- evaluate_plan(
  site_lai_plots_template,
  rules = list(
    site__ = selected_sites,
    xvar__ = c("cum_lai", "cum_elai", "lai_mean", "elai_mean")
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

current_plan <- bind_plans(
  pre_plan,
  spec_validation_plan,
  spec_validation_plot_plan,
  site_lai_plots_plan
)
current_config <- drake_config(current_plan)
make(current_plan, jobs = parallel::detectCores())

## make(current_plan, c("site_dt", "tidy_params_dt"))

## ed_plan_template <- drake_plan(
##   edresult = run_ed_site_ens(
##     site = "site__",
##     trait_values = ensemble_trait_list[[ens__]],
##     outdir_prefix = file.path("testsamples", ens__),
##     end_date = "2011-12-31"
##   )
## )

## ed_plan <- evaluate_plan(
##   ed_plan_template,
##   rules = list(
##     site__ = selected_sites,
##     ens__ = seq_len(nens)
##   )
## )

## plan <- bind_plans(pre_plan, ed_plan)
## dconf <- drake_config(plan)
## make(plan, jobs = parallel::detectCores())
