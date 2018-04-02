#' Extract PFTs from css file
#'
#' @param css_file_path Path to CSS file
#' @param delim File delimiter. Default = " "
#' @param reference PFT reference table (default = pft_lookup)
#' @export
get_pfts <- function(css_file_path, delim = " ", reference = pft_lookup) {
  css_data <- read.delim(css_file_path, header = T, sep = delim)
  css_pfts <- as.vector(unique(css_data["pft"]))
  dplyr::filter(reference, pft_num %in% css_pfts$pft)
}

#' ED PFT lookup table
#'
#' @export
pft_lookup <- tibble::tribble(
  ~pft_num, ~pft_name,
  6, "temperate.Northern_Pine",
  7, "temperate.Southern_Pine",
  8, "temperate.Late_Conifer",
  9, "temperate.Early_Hardwood",
  10, "temperate.North_Mid_Hardwood",
  11, "temperate.Late_Hardwood"
)
