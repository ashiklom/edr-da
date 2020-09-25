#' Run ED for a particular site-ensemble combination
#'
#' @param site (Character) Site at which to run ED2
#' @param trait_values (List) ED parameter values
#' @param outdir_prefix (Character) Name of output directory. Usually
#'   associated with an individual ensemble member (e.g. `ens_01`)
#' @param outdir_base (Character) Path to base output directory, which
#'   will hold all ensemble outputs. Default = `"ensemble_outputs"`.
#' @param sitedir_base (Character) Path to directory containing all site information.
#' @param start_date (Character, formatted as date) Start date of
#'   simulation. Defaults to June 6 in the year of the site input file.
#' @param end_date (Character, formatted as date) End date of
#'   simulation. Default is 2017-12-31.
#' @param soil_info_file (Character) Path to file containing soil
#'   texture information. Default is `"other_site_data/soil_texture.csv"`.
#' @param hetero (Logical) Whether or not the parameters use a
#'   heteroskedastic variance model. (Default = `FALSE`)
#' @param fix_allom2 (Logical) Whether or not to fix the second
#'   allometry coefficient. (Default = `TRUE`)
#' @param is_initial (Logical) Whether the run is a bare ground run
#'   (`TRUE`) or initial condition run (`FALSE`, default).
#' @param ed2in_template_file (Character) Path to ED2IN template.
#'   Default = `"inst/ED2IN"`.
#' @param ed2in_custom (List) Custom ED2IN settings for ED run
#'   (default = `list()`).
#' @param ed2_exe (Character) Path to ED2 executable. Passed directly
#'   to [base::system2()], so if it is in `$PATH`, this can just be
#'   the name of the executable. Default = `"ed2"`.
#' @importFrom magrittr %>%
#' @return List summarizing ED2 execution status.
#' @author Alexey Shiklomanov
#' @export
run_ed_site_ens <- function(site, trait_values, outdir_prefix,
                            outdir_base = "ensemble_outputs",
                            sitedir_base = "sites",
                            start_date = NULL,
                            end_date = "2017-12-31",
                            soil_info_file = file.path("other_site_data", "soil_texture.csv"),
                            hetero = FALSE,
                            fix_allom2 = TRUE,
                            is_initial = FALSE,
                            ed2in_template_file = system.file("ED2IN", package = "redr"),
                            ed2in_custom = list(),
                            ed2_exe = "ed2") {
  stopifnot(
    requireNamespace("PEcAn.ED2"),
    is.null(start_date) || is.character(start_date),
    is.character(end_date),
    is.character(ed2_exe),
    length(ed2_exe) == 1,
    file.exists(sitedir_base),
    file.exists(soil_info_file),
    file.exists(ed2in_template_file),
    is.list(ed2in_custom)
  )

  # Directory where outputs are saved
  outdir <- file.path(outdir_base, outdir_prefix, site)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Directory where site inputs (css, pss, site) files are stored
  site_input_dir <- file.path(sitedir_base, site)
  stopifnot(file.exists(site_input_dir))

  # Build site information from file names
  # NOTE: If multiple files present, use oldest (first) one
  site_info <- tibble::tibble(site_files = head(list.files(site_input_dir, "\\.css$"), 1)) %>%
    dplyr::mutate(
      prefix = stringr::str_match(site_files, "(^.*)lat.*")[, 2],
      full_prefix = file.path(site_input_dir, prefix),
      year = stringr::str_extract(prefix, "[[:digit:]]{4}") %>% as.numeric(),
      year = dplyr::if_else(is_initial, 1983, year),
      veg = purrr::map(full_prefix, PEcAn.ED2::read_ed_veg),
      latitude = purrr::map_dbl(veg, "latitude"),
      longitude = purrr::map_dbl(veg, "longitude"),
      start_date = if (!is.null(start_date)) start_date else paste0(year, "-06-02"),
      end_date = end_date
    )

  # Directory where ED meteorology (ED_MET_DRIVER_HEADER) is stored
  ed_met_dir <- file.path(site_input_dir, "ED_NARR")
  stopifnot(file.exists(file.path(ed_met_dir, "ED_MET_DRIVER_HEADER")))

  # Build site meteorology info table
  met_info <- tibble::tibble(
    file = file.path(ed_met_dir, "ED_MET_DRIVER_HEADER"),
    host = "",
    mimetype = "text/plain",
    formatname = "ed.met_driver_header files format",
    startdate = as.POSIXct(site_info$start_date, tz = "UTC"),
    enddate = as.POSIXct(site_info$end_date, tz = "UTC"),
    dbfile.name = "ED_MET_DRIVER_HEADER"
  )

  # Add soil information to input vegetation file
  soil_all <- readr::read_csv(soil_info_file)
  soil_site <- soil_all %>% dplyr::filter(sites == site)

  split_soil <- function(soil_df) {
    toprow <- botrow <- soil_df
    toprow$depth_bottom_cm <- toprow$depth_bottom_cm / 2
    botrow$depth_top_cm <- toprow$depth_bottom_cm
    dplyr::bind_rows(toprow, botrow)
  }

  nsoil <- nrow(soil_site)
  # HACK: ED complains if there's only one soil layer.
  # So, if there's only one, split it into two identical layers.
  if (nsoil == 1) {
    soil_site <- split_soil(soil_site)
    nsoil <- nrow(soil_site)
  }

  add_site_soil <- function(veg, soil_df) {
    site <- veg$site
    nsite <- attr(site, "nsite")
    site$soil <- NULL
    texture_vec <- soil_df %>%
      dplyr::left_join(texture_codes, by = "texture_code") %>%
      dplyr::select(ed_texture_id) %>%
      purrr::transpose() %>%
      purrr::flatten() %>%
      rev() %>%
      setNames(paste0("soil", seq(length(.))))
    site <- tibble::add_column(site, !!!texture_vec)
    attr(site, "file_format") <- 3
    attr(site, "nsite") <- nsite
    veg$site <- site
    veg
  }

  # Write input vegetation file
  veg_mod <- add_site_soil(site_info$veg[[1]], soil_site)
  veg_prefix <- file.path(outdir, "FFT.ensemble")
  PEcAn.ED2::write_ed_veg(veg_mod, veg_prefix)

  # Configure directory for writing ED inputs and outputs 
  ed_out_dir <- file.path(outdir, "out")
  dir.create(ed_out_dir, showWarnings = FALSE, recursive = TRUE)

  # Prepare ED2IN file
  ed2in_default <- list(
    IED_INIT_MODE = if (is_initial) 0 else 3,
    NZG = nsoil,  # number of soil layers
    SLZ = rev(-soil_site$depth_bottom_cm / 100),  # soil layer depths
    SLMSTR = rep(0.65, nsoil),
    STGOFF = rep(0, nsoil),
    UNITSTATE = 1, FRQSTATE = 1, # Write daily history files
    DTLSM = 300,  # Lower these values to reduce the likelihood of integration failures
    RADFRQ = 300,
    RK4_TOLERANCE = 1e-4,
    INTEGRATION_SCHEME = 1,
    H2O_PLANT_LIM = 0,    # Turn of plant hydraulic limitation
    IOOUTPUT = 0,
    PLANT_HYDRO_SCHEME = 0,
    ISTOMATA_SCHEME = 0,
    ISTRUCT_GROWTH_SCHEME = 0,
    TRAIT_PLASTICITY_SCHEME = 0,
    INCLUDE_THESE_PFT = pft_lookup$pft_num  # Limit to PDA PFTs
  )
  ed2in_custom <- modifyList(ed2in_default, ed2in_custom)
  ed2in <- PEcAn.ED2::read_ed2in(ed2in_template_file) %>%
    PEcAn.ED2::modify_ed2in(
      veg_prefix = veg_prefix,
      latitude = site_info$latitude,
      longitude = site_info$longitude,
      met_driver = met_info$file,
      start_date = site_info$start_date,
      end_date = site_info$end_date,
      add_if_missing = TRUE,
      output_dir = ed_out_dir,
      run_dir = outdir,
      EDI_path = normalizePath(file.path("ed-inputs", "EDI")),
      output_types = c("instant", "monthly", "restart"),
      runtype = "INITIAL",
      .dots = ed2in_custom
    )
  ed2in_path <- file.path(outdir, "ED2IN")
  PEcAn.ED2::write_ed2in(ed2in, ed2in_path, barebones = TRUE)

  # Prepare ED config.xml (parameter file)
  saveRDS(trait_values, file.path(outdir, "trait_values.rds"))
  xml_path <- file.path(outdir, "config.xml")
  PEcAn.logger::logger.setLevel("INFO")
  xml <- PEcAn.ED2::write.config.xml.ED2(
    defaults = list(),
    settings = list(model = list(revision = "git"), config.header = NULL),
    trait.values = trait_values
  )
  PREFIX_XML <- "<?xml version=\"1.0\"?>\n<!DOCTYPE config SYSTEM \"ed.dtd\">\n"
  XML::saveXML(xml, file = xml_path, indent = TRUE, prefix = PREFIX_XML)

  # Run ED
  exec_start <- Sys.time()
  ed_raw_output <- system2(
    ed2_exe,
    c("-s", "-f", ed2in_path),
    stdout = TRUE,
    stderr = TRUE
  )
  exec_stop <- Sys.time()

  # Return result object
  list(
    outdir = outdir,
    exec_start = exec_start,
    exec_stop = exec_stop,
    ed_output = ed_raw_output
  )
}

#' Convert raw `BayesianTools` samples object to an ED-compatible
#' parameter list
#'
#' @param samples_bt `BayesianTools` samples object, as returned by
#'   [BayesianTools::runMCMC()].
#' @param param_names (Character) Vector of parameter names
#' @param other_posteriors Additional posterior samples to use
#' @param nens (Integer) Size of the resulting ensemble
#' @param fix_allom2
#' @param last_n (Integer) How many samples to grab from `samples_bt`.
#' @return List length `nens`, each element a named list of PFT
#'   parameters.
#' @inheritParams run_ed_site_ens
#' @importFrom magrittr %>%
#' @author Alexey Shiklomanov
#' @export
preprocess_samples <- function(samples_bt, param_names, other_posteriors, nens,
                               fix_allom2 = TRUE,
                               last_n = 5000) {
  samples_raw <- BayesianTools::getSample(samples_bt, numSamples = last_n)
  stopifnot(NCOL(samples_raw) == length(param_names))
  colnames(samples_raw) <- param_names

  samples <- filter_samples(samples_raw)
  edr_sub <- samples[sample.int(nrow(samples), nens), ]

  # Combine with other posteriors
  other_sub <- other_posteriors[sample.int(nrow(other_posteriors), nens), ]
  combined <- cbind(edr_sub, other_sub)

  ensemble_trait_list <- apply(combined, 1, convert_ed_params, fix_allom2 = fix_allom2)
  ensemble_trait_list
}

#' Remove illegal (e.g. less than zero) values, and optionally remove
#' duplicate rows
#'
#' @param samples_raw (Numeric matrix) Matrix of sample values (MCMC
#'   iteration rows x parameter columns) 
#' @param dedup (Logical) If `TRUE`, remove duplicate rows
#' @return Numeric matrix of processed samples
#' @author Alexey Shiklomanov
#' @export
filter_samples <- function(samples_raw, dedup = TRUE) {
  # Remove duplicate rows
  if (dedup) {
    sdup <- duplicated(samples_raw)
    samples <- samples_raw[!sdup, , drop = FALSE]
  } else {
    samples <- samples_raw
  }

  # Some parameters must be >0
  # filter them here
  noneg_cols <- grep("prospect", colnames(samples))
  edr_nneg <- samples[, noneg_cols, drop = FALSE] < 0
  isneg <- rowSums(edr_nneg) > 0
  samples[!isneg, ]
}

#' Convert raw parameter values to ED2 input format
#'
#' First, convert PROSPECT parameter values to spectra and aggregate
#' up to visible (VIS) and near-infrared (NIR) range. Second, if
#' `fix_allom2`, retrieve the allometry mean and add it to the sample list.
#'
#' @param params (Numeric) Long parameter vector, as sampled by
#'   `BayesianTools`.
#' @param vis (Numeric) Visible wavelength range, in nm (default = `400:700`)
#' @param nir (Numeric) Near-infrared wavelength range, in nm (default
#'   = `701:1300`)
#' @inheritParams run_ed_site_ens
#' @author Alexey Shiklomanov
convert_ed_params <- function(params, vis = 400:700, nir = 701:1300, fix_allom2 = TRUE) {
  params <- params[!grepl("sitesoil", names(params))]
  pedr <- PEcAnRTM::params2edr(params)
  pspec <- pedr$spectra_list
  traits <- pedr$trait.values
  if (fix_allom2) {
    # Set allometry exponent to posterior mean
    traits <- purrr::map2(
      traits,
      purrr::map(allom_mu, "b2Bl")[names(traits)],
      ~`[<-`(.x, "b2Bl", .y)
    )
  }
  pvis <- purrr::map(pspec, ~colMeans(.[[vis]])) %>%
    purrr::map(setNames, c("leaf_reflect_vis", "leaf_trans_vis")) %>%
    .[names(traits)]
  pnir <- purrr::map(pspec, ~colMeans(.[[nir]])) %>%
    purrr::map(setNames, c("leaf_reflect_nir", "leaf_trans_nir")) %>%
    .[names(traits)]
  stopifnot(
    all(names(traits) == names(pvis)),
    all(names(traits) == names(pnir))
  )
  traits %>%
    purrr::map2(pvis, c) %>%
    purrr::map2(pnir, c)
}

#' Generate EDR spectra (mean, SD, and confidence envelope)
#' predictions for a single site
#'
#' @param params_matrix
#' @param site
#' @param nsamp
#' @inheritParams filter_samples
#' @return
#' @author Alexey Shiklomanov
#' @export
predict_site_spectra <- function(params_matrix, site,
                                 nsamp = 500,
                                 dedup = FALSE,
                                 progress = FALSE,
                                 run_config = "homo-pooled") {
  site_list <- readLines("other_site_data/site_list")
  if (!site %in% site_list) {
    stop("Site ", site, " is not in site_list")
  }
  site_css_file <- list.files(
    file.path("sites", site),
    "\\.css$",
    full.names = TRUE
  ) %>% head(1)
  stopifnot(length(site_css_file) > 0, file.exists(site_css_file))
  site_data <- read.table(site_css_file, header = TRUE)
  isite <- which(site == site_list)

  params_filtered_all <- filter_samples(params_matrix, dedup = dedup)
  isamp <- sample.int(min(nsamp, NROW(params_filtered_all)))
  params_filtered <- params_filtered_all[isamp, ]
  waves <- seq(400, 1300, 5)
  pb <- if (progress) progress::progress_bar$new(total = NROW(params_filtered)) else NULL
  result <- apply(params_filtered, 1, run_edr_sample, isite = isite, site_data = site_data,
                  pb = pb, wavelengths = waves)
  nulls <- vapply(result, is.null, logical(1))
  result <- result[!nulls]
  if (grepl("homo-pooled", run_config)) {
    rint <- params_filtered[!nulls, "residual"]
    rslope <- 0
  } else if (grepl("hetero-pooled", run_config)) {
    rint <- params_filtered[!nulls, "residual_intercept"]
    rslope <- params_filtered[!nulls, "residual_slope"]
  } else if (grepl("homo-sitespecific", run_config)) {
    rint <- params_filtered[!nulls, paste0("residual", isite)]
    rslope <- 0
  } else if (grepl("hetero-sitespecific", run_config)) {
    rint <- params_filtered[!nulls, paste0("residual_intercept", isite)]
    rslope <- params_filtered[!nulls, paste0("residual_slope", isite)]
  } else {
    stop("Invalid run config: ", run_config)
  }
  alb_list <- purrr::map(result, "albedo")
  albedos <- purrr::pmap_dfr(
    list(alb_list, rslope, rint),
    ~tibble::tibble(
      wavelength = waves,
      albedo = .x,
      albedo_r = rnorm(length(waves), ..1, ..1 * ..2 + ..3)
    )
  )
  output <- albedos %>%
    dplyr::group_by(wavelength) %>%
    dplyr::summarize_if(
      is.numeric,
      list(
        mean = ~mean(.),
        sd = ~sd(.),
        q025 = ~quantile(., 0.025),
        q975 = ~quantile(., 0.975)
      )
    )
  output
}

#' Run EDR at a particular site, with a particular set of parameters
#'
#' @param params (Numeric) Long parameter vector, such as that used by
#'   `BayesianTools`
#' @param isite (Integer) Site index
#' @param site_data (data.frame) Site cohort data, as `data.frame`.
#'   See [PEcAn.ED2::read_css()].
#' @param npft (Integer) Number of PFTs. Defaults to number of rows in
#'   [pft_lookup] that are not "Southern Pine"
#' @param npft_param (Integer) Number of parameters per PFT. Default = 10.
#' @param direct_sky_frac (Numeric) Direct sky fraction. See
#'   [edr_r()]. Default = 0.9.
#' @param czen (Numeric) Cosine of solar zenith angle. See [edr_r()].
#'   Default = 1.
#' @param pb (Progress bar) Optional progress bar (default = `NULL`)
#' @inheritParams edr_r
#' @return Output of EDR model. See [edr_r()].
#' @author Alexey Shiklomanov
#' @export
run_edr_sample <- function(params, isite, site_data,
                           fix_allom2 = TRUE,
                           npft = NROW(
                             dplyr::filter(pft_lookup,
                                           !grepl("Southern_Pine", pft_name))
                           ),
                           npft_param = 10,
                           nsite = NULL,
                           direct_sky_frac = 0.9,
                           czen = 1,
                           wavelengths = 400:2500,
                           pb = NULL) {
  if (!fix_allom2) stop("Only fixed allometry currently supported.")
  if (!is.null(pb)) pb$tick()
  stopifnot(isite %% 1 == 0)
  has_names <- !is.null(names(params))
  stopifnot(has_names)
  pft <- get_ipft(site_data[["pft"]])
  dbh <- site_data[["dbh"]]
  nplant <- site_data[["n"]]
  ncohort <- length(dbh)

  # Calculate heights and height order (shortest first)
  hite <- dbh2h(dbh, pft)
  ihite <- order(hite, decreasing = FALSE)

  # Order cohorts by decreasing height (tallest first)
  dbh <- dbh[ihite]
  pft <- pft[ihite]
  nplant <- nplant[ihite]
  hite <- hite[ihite]

  # Extract site-specific soil moisture
  soil_moisture <-  params[grepl(paste0("sitesoil_", isite, "$"), names(params))]

  # Extract residuals
  rslope <- params[paste0("residual_intercept", isite)]
  rint <- params[paste0("residual_intercept", isite)]

  # Remaining parameters are PFT-specific
  pft_params_v <- params[!grepl("residual|sitesoil", names(params))]

  # Create a matrix nparam (rows) x npft (cols)
  pft_params <- matrix(pft_params_v, ncol = npft)

  # Extract parameters
  N <- pft_params[1, ]
  Cab <- pft_params[2, ]
  Car <- pft_params[3, ]
  Cw <- pft_params[4, ]
  Cm <- pft_params[5, ]
  SLA <- pft_params[6, ]
  b1Bl <- pft_params[7, ]
  b1Bw <- pft_params[8, ]
  clumping_factor <- pft_params[9, ]
  orient_factor <- pft_params[10, ]

  b2Bl <- purrr::map_dbl(allom_mu, "b2Bl")
  b2Bw <- purrr::map_dbl(wallom_mu, "b2Bw")

  # Calculate allometries
  bleaf <- size2bl(dbh, b1Bl[pft], b2Bl[pft])
  lai <- nplant * bleaf * SLA[pft]
  wai <- wai_allometry(dbh, nplant, b1Bw[pft], b2Bw[pft])

  # Cohort area index is constant (no crown radius model)
  cai <- rep(1, ncohort)

  # Run EDR
  tryCatch(
    error = function(e) NULL,
    edr_r(pft, lai, wai, cai,
          N, Cab, Car, Cw, Cm,
          orient_factor, clumping_factor,
          soil_moisture,
          direct_sky_frac,
          czen,
          wavelengths = wavelengths)
  )
}

pft_factor <- function(pft) {
  lvl <- c("Early_Hardwood", "North_Mid_Hardwood", "Late_Hardwood",
           "Northern_Pine", "Late_Conifer")
  lbl <- c("EH", "MH", "LH", "NP", "LC")
  factor(pft, lvl, lbl)
}
