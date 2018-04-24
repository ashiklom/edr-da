#' Read and organize ED history file
#'
#' Reads an ED history file into a list containing `data.frames` with similar 
#' themes (e.g. PFT, cohort, soil).
#'
#' - `datetime` -- Date and time of history file
#'
#' - `scalar`: Two columns -- `variable` and `value`
#'
#' - `cohort`: Each row is one cohort. List columns are as follows:
#'  - `CB_*` -- Carbon budget. Maximum values over last 12 months + current 
#'  month / partial month integration
#'  - `RAD_PROFILE_CO` -- Radiation profile, 10 values, organized as follows: 
#'  (PAR, NIR, TIR) x (beam, diff) x (down, up). E.g. PAR beam down, PAR beam 
#'  up, PAR diff down, PAR diff up, NIR beam down, ...
#'  - `MORT_RATE_CO` -- Rate of mortality. Each value is a different mortality 
#'  type: (1) Aging, (2) Negative Carbon, (3) Treefall, (4) Cold, (5) Disturbance.
#'
#' - `pft`: Each row is one PFT. List columns are as follows:
#'  - Length 11 -- By DBH class, in 10 cm increments. See `dbh_class`.
#'  - Length 8 -- Height classes. See `hgt_class`
#'
#' - `soil`: Each row is a soil depth, sorted from deepest to shallowest
#'
#' - `other`: Additional data that doesn't fit into other categories.
#'
#' - `dbh_class` -- Vector of DBH classes
#'
#' - `hgt_class` -- Vector of height classes
#'
#' - `disturbance_types` -- Vector of disturbance types
#'
#' @param historyfile ED `history-S` file
#' @return List containing scalar, cohort, PFT, and soil data, and other metadata. See Details.
#' @export
read_ed_history <- function(historyfile) {
  stopifnot(file.exists(historyfile))
  datetime_str <- stringr::str_extract(
    basename(historyfile),
    "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}-[[:digit:]]{6}"
  )
  datetime <- lubridate::ymd_hms(datetime_str)
  # First group of values are hard-coded in ED source
  # DBH classes
  dbh_class_min <- seq(0, 100, by = 10)
  dbh_class_max <- c(dbh_class_min[-1], Inf)
  dbh_class <- paste(dbh_class_min, dbh_class_max, sep = "-")
  ndbh <- length(dbh_class)
  stopifnot(ndbh == 11)

  # Disturbance types
  disturbance_types <- c("Clear cut", "Forest plantation", "Tree fall", "Fire",
                         "Forest regrowth", "Logged forest")

  # Load the history file
  hf <- hdf5r::H5File$new(historyfile, "r")
  on.exit(hf$close_all())

  # Extract "metadata" objects
  ncohort <- hf[["NCOHORTS_GLOBAL"]][]
  nsoil <- hf[["NZG"]][]

  # Height classes
  hgt_class_min <- hf[["HGT_CLASS"]][]
  hgt_class_max <- c(hgt_class_min[-1], Inf)
  hgt_class <- paste(hgt_class_min, hgt_class_max, sep = "-")

  # Extract dimensions
  dims <- hf$ls()$dataset.dims %>%
    strsplit(" x ") %>%
    purrr::map(as.numeric) %>%
    setNames(names(hf))

  message("Reading cohort data")
  # Simple vector data
  cohort_vec <- dims %>%
    purrr::keep(~length(.) == 1 & .[1] == ncohort) %>%
    purrr::imap_dfc(~read_hf_var(hf, .y))
  # Matrix data, coerced to list columns
  cohort_mat <- dims %>%
    purrr::keep(~length(.) == 2) %>%
    purrr::keep(~.[2] == ncohort) %>%
    purrr::imap(~read_hf_var(hf, .y)) %>%
    purrr::map(~split(., col(.))) %>%
    purrr::transpose() %>%
    purrr::modify_depth(2, list) %>%
    dplyr::bind_rows()
  cohort_all <- dplyr::bind_cols(cohort_vec, cohort_mat) %>%
    dplyr::mutate(cohort_id = row_number()) %>%
    dplyr::select(cohort_id, dplyr::everything())

  message("Reading scalar data")
  scalars <- dims %>%
    purrr::keep(~all(. == 1)) %>%
    purrr::imap_dbl(~read_hf_var(hf, .y)) %>%
    tibble::enframe("variable", "value")

  message("Reading PFT data")
  pft_vec <- dims %>%
    purrr::keep(~.[1] == npft) %>%
    purrr::keep(~prod(.) == npft) %>%
    purrr::imap_dfc(~read_hf_var(hf, .y))
  pft_mat <- dims %>%
    purrr::keep(~.[1] == npft) %>%
    purrr::keep(~prod(.) > npft) %>%
    purrr::imap(~read_hf_var(hf, .y)) %>%
    purrr::map(t) %>%
    purrr::map(~split(., col(.))) %>%
    purrr::transpose() %>%
    purrr::modify_depth(2, list) %>%
    dplyr::bind_rows()
  pft_all <- dplyr::bind_cols(pft_vec, pft_mat) %>%
    dplyr::mutate(pft_id = row_number()) %>%
    dplyr::select(pft_id, dplyr::everything())

  message("Reading soil data")
  soil_all <- dims %>%
    purrr::keep(~.[1] == nsoil) %>%
    purrr::imap_dfc(~read_hf_var(hf, .y)) %>%
    dplyr::mutate(soil_id = row_number()) %>%
    dplyr::select(soil_id, dplyr::everything())

  left <- dims %>%
    .[!names(.) %in% colnames(cohort_all)] %>%
    .[!names(.) %in% scalars$variable] %>%
    .[!names(.) %in% colnames(pft_all)] %>%
    .[!names(.) %in% colnames(soil_all)]

  other <- purrr::imap(left, ~read_hf_var(hf, .y))

  list(
    datetime = datetime,
    scalar = scalars,
    cohort = cohort_all,
    pft = pft_all,
    soil = soil_all,
    other = other,
    dbh_class = dbh_class,
    hgt_class = hgt_class,
    disturbance_types = disturbance_types
  )
}

#' Map-friendly way to read a variable from an HDF5 file object
#'
#' @param hf `H5File` object
#' @param variable Character string of variable to read
read_hf_var <- function(hf, variable) hf[[variable]]$read()
