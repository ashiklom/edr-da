#--------------------------------------------------------------------------------------------------#
source('config.R')

library(tidyverse)
library(redr)

data(css_ex1)
data(pss_ex1)
data(site_ex1)

save_plot <- FALSE #TRUE/FALSE
hidden <- TRUE #TRUE/FALSE
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (hidden) {
  prefix <- '.edr_prior_sensitivity'
} else {
  prefix <- 'edr_prior_sensitivity'
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

############################################################
# Load priors
load('priors/sunlit_meanjp.RData')

clumping_factor_prior <- function(n = 1) runif(n, 0, 0.5)
orient_factor_prior <- function(n = 1) runif(n, -0.5, 0.5)

sample_params_pft <- function(pfts) {
    npft <- length(pfts)
    nmeans <- ncol(means)
    drawmeans <- matrix(numeric(), npft, nmeans)
    dimnames(drawmeans) <- list(pfts, colnames(means))
    for (i in seq_len(npft)) {
        drawmeans[i,] <- mvtnorm::rmvnorm(1, means[pfts[i],], covars[pfts[i],,])[1,]
    }
    draw_df <- as_data_frame(drawmeans) %>% 
        mutate(pft = rownames(drawmeans),
               clumping_factor = clumping_factor_prior(npft),
               orient_factor = orient_factor_prior(npft),
               sla = 1 / leaf_mass_per_area) %>% 
        select(pft, everything()) %>% 
        select(-leaf_mass_per_area)
    return(draw_df)
}

nsamp <- 100
pfts <- c('temperate.Early_Hardwood', 'temperate.North_Mid_Hardwood', 'temperate.Late_Hardwood')
albedo <- matrix(numeric(), 2101, nsamp)
for (i in seq_len(nsamp)) {
    print(i)
    params <- sample_params_pft(pfts)
    arg_list <- params_df2list(params, prospect_version = 5)
    #arg_list$par.wl <- 400:2499
    #arg_list$nir.wl <- 2500
    albedo[,i] <- run_edr(prefix, edr_args = arg_list)
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
wl <- 400:2500
mu <- rowMeans(albedo)
lo <- apply(albedo, 1, quantile, 0.025)
hi <- apply(albedo, 1, quantile, 0.975)

if (save_plot) {
  png(file.path(prefix,'prior_sensitivity_test_albedo.png'),width=4900, height =2700,res=400)
  par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
  matplot(wl, albedo, type='l', lty = 'dashed', col = 'grey')
  lines(wl, hi, col = 'red', lwd = 2)
  lines(wl, mu, col = 'black', lwd = 2)
  lines(wl, lo, col = 'red', lwd = 2)
  box(lwd=2.2)
  dev.off()
} else {
  matplot(wl, albedo, type='l', lty = 'dashed', col = 'grey')
  lines(wl, hi, col = 'red', lwd = 2)
  lines(wl, mu, col = 'black', lwd = 2)
  lines(wl, lo, col = 'red', lwd = 2)
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF