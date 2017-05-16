source('common.R')

library(tidyverse)

prefix <- 'bl_same_dbh'

css_bl_same_dbh <- extend_df(css_df, cohort = 1:3, dbh = 30, pft = c(9, 10, 11))
genrun <- generate_run(prefix = prefix,
                       site_lat = site_lat,
                       site_lon = site_lon,
                       site_df = site_df,
                       pss_df = pss_df,
                       css_df = css_bl_same_dbh,
                       common_inputs_dir = common_inputs_dir,
                       site_met_dir = site_met_dir,
                       ed_exe_path = ed_exe_path,
                       RMDIR = TRUE)

run_ed(prefix)

edr_setup <- setup_edr(prefix, edr_exe_path = edr_exe_path)

params <- tribble(
    ~pft, ~N, ~Cab, ~Cw, ~Cm, ~orient_factor, ~clumping_factor,
    'temperate.Early_Hardwood', 1.4, 40, 0.01, 0.01, 0.5, 0.5, 
    'temperate.North_Mid_Hardwood', 1.4, 30, 0.01, 0.01, 0.5, 0.5,
    'temperate.Late_Hardwood', 1.4, 20, 0.01, 0.01, 0.5, 0.5
                                   )

# TODO: Sample from the prior distribution

sensitivity <- function(base_params, pft, parameter, ...) {
    param_seq <- seq(...)
    nparams <- length(param_seq)
    albedo_mat <- matrix(0, nrow = 2101, ncol = nparams)
    for (i in seq_len(nparams)) {
        temp_params <- set_param(params, pft, parameter, param_seq[i])
        albedo_mat[,i] <- run_edr(prefix, edr_args = params_df2list(temp_params))
    }
    return(albedo_mat)
}

# Some examples
early_N <- sensitivity(params, 'temperate.Early_Hardwood', 'N', from = 1.1, to = 2.0, length.out = 10)
mid_N <- sensitivity(params, 'temperate.North_Mid_Hardwood', 'N', from = 1.1, to = 2.0, length.out = 10)

