#--------------------------------------------------------------------------------------------------#
source('config.R')

library(tidyverse)
library(redr)

data(css_ex1)
data(pss_ex1)
data(site_ex1)

save_plot <- TRUE #TRUE/FALSE
hidden <- FALSE #TRUE/FALSE
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (hidden) {
  prefix <- '.edr_sensitivity'
} else {
  prefix <- 'edr_sensitivity'
}

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

message('Running ED...')
runed <- run_ed(prefix)
tail(runed)
message('Done!')

edr_setup <- setup_edr(prefix, edr_exe_path = edr_exe_path)

params <- tribble(
    ~pft, ~N, ~Cab, ~Cw, ~Cm, ~orient_factor, ~clumping_factor,
    'temperate.Early_Hardwood', 1.4, 40, 0.01, 0.01, 0.5, 0.8, 
    'temperate.North_Mid_Hardwood', 1.4, 30, 0.01, 0.01, 0.5, 0.8,
    'temperate.Late_Hardwood', 1.4, 20, 0.01, 0.01, 0.5, 0.8
    )

sensitivity <- function(base_params, pft, parameter, param_seq) {
    nparams <- length(param_seq)
    albedo_mat <- matrix(0, nrow = 2101, ncol = nparams)
    for (i in seq_len(nparams)) {
        # This changes the value of one parameter for one PFT from a properly 
        # formatted parameter data.frame (like the `params` data.frame above)
        temp_params <- set_param(base_params, pft, parameter, param_seq[i])
        arg_list <- params_df2list(temp_params)
        #arg_list$par.wl <- 400:2499
        #arg_list$nir.wl <- 2500
        albedo_mat[,i] <- run_edr(prefix, edr_args = arg_list)
    }
    return(albedo_mat)
}

# Run sensitivity analysis for PROSPECT N parameter
N_seq <- seq(from = 1.1, to = 3, length.out = 20)
early_N <- sensitivity(params, 'temperate.Early_Hardwood', 'N', N_seq)
mid_N <- sensitivity(params, 'temperate.North_Mid_Hardwood', 'N', N_seq)
late_N <- sensitivity(params, 'temperate.North_Mid_Hardwood', 'N', N_seq)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (save_plot) { 
  png(file.path(prefix,'Albedo_sensitivity_analysis_PROSPECT_N_param.png'),width=3500, height=1300, res=200)
  par(mfrow=c(1,3), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
  matplot(early_N, type='l', main = 'Early Hardwood')
  box(lwd=2.2)
  matplot(mid_N, type='l', main = 'Mid Hardwood')
  box(lwd=2.2)
  matplot(late_N, type='l', main = 'Late Hardwood')
  box(lwd=2.2)
  dev.off()
} else {
  par(mfrow = c(1,3))
  matplot(early_N, type='l', main = 'Early Hardwood')
  matplot(mid_N, type='l', main = 'Mid Hardwood')
  matplot(late_N, type='l', main = 'Late Hardwood')
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF
