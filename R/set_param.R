#' @export
set_param <- function(params_df, pft, param, value) {
    params_df[[param]][params_df[['pft']] == pft] <- value
    return(params_df)
}
