#' Extend a one-row data.frame with more elements
#' 
#' @param df data.frame to extend. Must have column names and only have one row.
#' @param ... Values with which to extend `df`. Argument names are column names in `df`.
#' @export
extend_df <- function(df, ...) {
    stopifnot(nrow(df) == 1, !is.null(colnames(df)))
    collist <- list(...)
    df_list <- as.list(df)
    for (i in seq_along(collist)) {
        df_list[[names(collist)[i]]] <- collist[[i]]
    }
    return(data.frame(df_list))
}
