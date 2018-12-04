#' Get PFT index (starting from 1) of a PFT vector
#'
#' Thin wrapper around [base::match] with some additional error checking.
#' 
#' @param pft Integer vector of original (ED) PFT numbers
#' @param pft_dict Integer vector describing the order of PFTs (first
#'   PFT is index 1, second PFT is index 2, etc.)
#' @return Vector of PFTs
#' @author Alexey Shiklomanov
#' @export
get_ipft <- function(pft, pft_dict = c(9, 10, 11, 6, 8)) {
  ipft <- match(pft, pft_dict)
  if (any(is.na(ipft))) {
    stop(
      "Problem with PFTs.\n",
      "Original PFT vector: ", paste(pft, collapse = ", "), "\n",
      "PFT dictionary: ", paste(pft_dict, collapse = ", "), "\n",
      "Resulting PFTs: ", paste(ipft, collapse = ", ")
    )
  }
  ipft
}
