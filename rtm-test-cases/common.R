library(redr)

inputs_path <- normalizePath('../ed-inputs')
common_inputs_dir <- file.path(inputs_path, 'EDI')
site_met_dir <- file.path(inputs_path, 'met3', 'US-WCr')
ed_exe_path <- '/home/ashiklom/Projects/ED2/ED/build/ed_2.1'
edr_exe_path <- '/home/ashiklom/Projects/ED2/EDR/build/ed_2.1'

site_lat <- 45.5
site_lon <- -90.5
datetime <- as.POSIXlt('2004-07-01 12:00:00')

data(css_ex1)
data(pss_ex1)
data(site_ex1)

params_df2list <- function(params_df, pftcol = 'pft', spectra_col = 'spectra_list', trait_col = 'trait.values') {
    edr_df <- params_df %>% 
        mutate(params = pmap(list(N, Cab, Cw, Cm), c),
               spectra_list = map(params, PEcAnRTM::prospect, version = 4, include.wl = TRUE),
               trait.values = pmap(list(orient_factor = orient_factor, 
                                        clumping_factor = clumping_factor), list))
    spectra_list <- edr_df[[spectra_col]]
    trait.values <- edr_df[[trait_col]]
    names(spectra_list) <- names(trait.values) <- edr_df[[pftcol]]
    list(spectra_list = spectra_list, trait.values = trait.values, datetime = datetime)
}

set_param <- function(params_df, pft, param, value) {
    params_df[[param]][params_df[['pft']] == pft] <- value
    return(params_df)
}
