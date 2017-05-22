#' @export
params_df2list <- function(params_df, prospect_version = 4, pftcol = 'pft', datetime = as.POSIXlt('2004-07-01 12:00:00')) {
    stopifnot(pftcol %in% colnames(params_df))
    prospect_params_rxp <- 'N|Cab|Cbrown|Canth|Car|Cw|Cm'
    prospect_cols <- grep(sprintf('^(%s)', prospect_params_rxp), colnames(params_df), value = TRUE)
    stopifnot(length(prospect_cols) >= 4)
    edr_cols_rxp <- sprintf('^(%1$s|prospect_params|%2$s)', pftcol, prospect_params_rxp)
    edr_cols <- grep(edr_cols_rxp, colnames(params_df), value = TRUE, invert = TRUE)
    edr_df <- params_df %>% 
        mutate(prospect_params = pmap(.[prospect_cols], list),
               spectra_list = map(prospect_params, PEcAnRTM::prospect, version = prospect_version, include.wl = TRUE),
               trait.values = pmap(.[edr_cols], list))
    spectra_list <- edr_df[['spectra_list']]
    trait.values <- edr_df[['trait.values']]
    names(spectra_list) <- names(trait.values) <- edr_df[[pftcol]]
    out_list <- list(spectra_list = spectra_list, trait.values = trait.values, datetime = datetime)
    return(out_list)
}

