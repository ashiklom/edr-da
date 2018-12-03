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
    requireNamespace("redr"),
    requireNamespace("PEcAn.ED2"),
    requireNamespace("PEcAnRTM"),
    requireNamespace("magrittr"),
    requireNamespace("dplyr"),
    requireNamespace("readr"),
    requireNamespace("tibble"),
    requireNamespace("stringr"),
    requireNamespace("purrr"),
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
      dplyr::left_join(redr::texture_codes, by = "texture_code") %>%
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
    INCLUDE_THESE_PFT = redr::pft_lookup$pft_num  # Limit to PDA PFTs
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


preprocess_samples <- function(samplefile, param_names, other_posteriors, nens,
                               fix_allom2 = TRUE,
                               last_n = 5000) {
  stopifnot(file.exists(samplefile))
  samples_bt <- readRDS(samplefile)
  samples_raw <- BayesianTools::getSample(samples_bt, numSamples = last_n)
  stopifnot(NCOL(samples_raw) == length(param_names))
  colnames(samples_raw) <- param_names

  # Remove duplicate rows
  sdup <- duplicated(samples_raw)
  samples <- samples_raw[!sdup, ]

  # Some parameters must be >0
  # filter them here
  noneg_cols <- grep("prospect", colnames(samples))
  edr_nneg <- samples[, noneg_cols] < 0
  isneg <- rowSums(edr_nneg) > 0
  edr_samples <- samples[!isneg, ]
  edr_sub <- samples[sample.int(nrow(samples), nens), ]

  # Combine with other posteriors
  other_sub <- other_posteriors[sample.int(nrow(other_posteriors), nens), ]
  combined <- cbind(edr_sub, other_sub)

  sum4ed <- function(params, vis = 400:700, nir = 701:1300) {
    params <- params[!grepl("sitesoil", names(params))]
    pedr <- PEcAnRTM::params2edr(params)
    pspec <- pedr$spectra_list
    traits <- pedr$trait.values
    if (fix_allom2) {
      # Set allometry exponent to posterior mean
      traits <- purrr::map2(
        traits,
        purrr::map(redr::allom_mu, "b2Bl")[names(traits)],
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

  ensemble_trait_list <- apply(combined, 1, sum4ed)
  ensemble_trait_list
}
