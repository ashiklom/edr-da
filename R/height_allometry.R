#' Given DBH and PFT identity, calculate cohort height
#'
#' @param dbh (Numeric) Vector of cohort DBH values (must be same
#'   length as `pft`)
#' @param pft (Integer) Vector of cohort PFT identities (must be same
#'   length as `dbh`)
#' @param b1Ht (Numeric) PFT-specific height allometry intercept
#'   (default defined in function)
#' @param b2Ht (Numeric) PFT-specific
#'   height allometry slope (default defined in function)
#' @return Numeric vector of cohort heights
#' @author Alexey Shiklomanov
#' @export
dbh2h <- function(dbh, pft, b1Ht = default_b1Ht, b2Ht = default_b2Ht) {
  # Order is: Early, Mid, Late Hardwood, North pine, Late conifer
  # Based on ED/src/init/ed_params.f90:2576
  default_b1Ht <- c(22.6799, 25.18, 23.3874, 27.14, 22.79)
  default_b2Ht <- c(-0.06534, -0.04964, -0.05404, -0.03884, -0.04445)
  stopifnot(
    all(pft <= 5), all(pft > 0),
    length(dbh) == length(pft),
    length(b1Ht) == 5, length(b2Ht) == 5
  )
  hgt_ref <- 1.3                        # Reference height
  hgt_ref + b1Ht[pft] * (1 - exp(b2Ht[pft] * dbh))
}
