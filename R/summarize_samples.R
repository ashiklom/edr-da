#' Summarize samples as a tidy data frame
#'
#' @param samples MCMC samples object
#' @export
summarize_samples <- function(samples) {
  samples_summary_raw <- summary(samples)
  samples_summary_raw %>%
    .[c("quantiles", "statistics")] %>%
    purrr::map(~tibble::as_tibble(., rownames = "parameter")) %>%
    purrr::reduce(dplyr::left_join) %>%
    split_params("parameter")
}

#' Split parameters column into categories
#'
#' @param .data Data frame containing parameter colum
#' @param param_col Name of parameter column
#' @param ... Additional arguments to [tidyr::separate]
#' @export
split_params <- function(.data, param_col, ...) {
  tidyr::separate(.data, !!param_col, c("biome", "pft", "variable"),
                  sep = "\\.", fill = "left", ...)
}

#' Prior and posterior density plots
#'
#' @param post Posterior density object
#' @param pri Prior density object
#' @param param Parameter name
#' @export
prior_posterior_density <- function(post, pri, param, ...) {
  xrange <- range(c(post$x, pri$x))
  yrange <- range(c(post$y, pri$y))
  plot(0, 0, type = "n", xlim = xrange, ylim = yrange,
       main = param, xlab = "value", ylab = "density")
  lines(post$x, post$y, col = "black")
  lines(pri$x, pri$y, col = "red")
  legend("topright", legend = c("posterior", "prior"),
         lty = "solid", col = c("black", "red"))
}

#' Compute column kernel densities
#'
#' @param mat Matrix on whose columns to compute densities
#' @param alpha Confidence level for trimming. Default = 0.05
#' @export
col_densities <- function(mat, alpha = 0.05, ...) {
  mat_list <- split(mat, col(mat)) %>% setNames(colnames(mat))
  los <- purrr::map(mat_list, quantile, alpha / 2)
  his <- purrr::map(mat_list, quantile, 1 - (alpha / 2))
  mat_list2 <- purrr::pmap(
    list(mat_list, los, his),
    ~..1[..1 > ..2 & ..1 < ..3]
  )
  purrr::map(mat_list2, density, ...)
}
