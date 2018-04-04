#' Mapping to ED PFTs
#'
#' @export
ed_pft_mapping <- tibble::tribble(
  ~pft_name, ~ed_pft_number,
  "temperate.Early_Hardwood", 9,
  "temperate.North_Mid_Hardwood", 10,
  "temperate.Late_Hardwood", 11,
  "temperate.Northern_Pine", 6,
  "temperate.Late_Conifer", 8
) %>%
  dplyr::mutate(pft_name = factor(pft_name, pft_name))
