#--------------------------------------------------------------------------------------------------#
source('config.R')

library(magrittr)
library(redr)

data(css_ex1)
data(pss_ex1)
data(site_ex1)

plot_albedo <- TRUE #TRUE/FALSE
generate_summary_figs <- TRUE #TRUE/FALSE
hidden <- FALSE  #TRUE/FALSE
PEcAn.logger::logger.setLevel("INFO")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (hidden) {
  prefix <- '.edr_inversion'
} else {
  prefix <- paste("edr_inversion", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
  #prefix <- 'edr_inversion'
}
PEcAn.logger::logger.info(paste0("Running inversion in dir: ",prefix))

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

# Params vector
# Temperate.Early_Hardwood: 1:8
#   1:5 -- PROSPECT params
#   6 -- sla
#   7 -- clumping_factor
#   8 -- orient_factor
# Temperate.North_Mid_Hardwood: 9:16
# Temperate.Late_Hardwood: 17:24
pft_end <- c(temperate.Early_Hardwood = 8,
             temperate.North_Mid_Hardwood = 16,
             temperate.Late_Hardwood = 24)

param_sub <- function(i, params) {
        param_seq <- (pft_end[i] - 7):pft_end[i]
        param_sub <- params[param_seq]
        return(param_sub)
}

prior_function <- function(params) {
    prior <- numeric(1)
    for (i in seq_along(pft_end)) {
        pft <- names(pft_end)[i]
        param_sub <- param_sub(i, params)
        # PROSPECT prior
        prior <- prior + mvtnorm::dmvnorm(c(param_sub[1:5], (1/param_sub[6])*1000), means[pft,], covars[pft,,],
                                          log = TRUE)
        # ED priors
        #prior <- prior + dunif(param_sub[7], 0, 1, TRUE) + dunif(param_sub[8], -1, 1, TRUE)
        prior <- prior + dunif(param_sub[7], 0, 1, TRUE) + dunif(param_sub[8], -0.5, 0.5, TRUE)
    }
    return(prior)
}

# Static initial conditions
# Can be replaced with a function that draws these values from distributions
#inits_function <- function() {
    # N, Cab, Car, Cw, Cm, SLA, clumping, orient
#    c(1, 35, 5, 0.006, 0.005, 15, 0.5, 0,     # Early
#      1, 35, 5, 0.006, 0.005, 15, 0.5, 0,     # Mid
#      1, 35, 5, 0.006, 0.005, 15, 0.5, 0)    # Late
#}

# inits_function <- function() {
#                   ## N, Cab, Car, Cw, Cm, SLA, clumping, orient
# vals <- rnorm(24, c(1.4, 35, 5, 0.006, 0.005, 15, 0.5, 0,     # Early
#                     1.4, 35, 5, 0.006, 0.005, 15, 0.5, 0,     # Mid
#                     1.4, 35, 5, 0.006, 0.005, 15, 0.5, 0), 0.001) # Late
# names(vals) <- rep(c('N', 'Cab', 'Car', 'Cw', 'Cm', 'SLA', 'clumping_factor', 'orient_factor'),3)
# return(vals)
# }

inits_function <- function() {
  samples <- numeric()
  for (i in seq_along(pft_end)) {
    pft <- names(pft_end)[i]
    # PROSPECT prior
    samples <- c(samples, mvtnorm::rmvnorm(1, means[pft,], covars[pft,,]))
    # ED priors
    samples <- c(samples, runif(1, 0, 1), runif(1, -0.5, 0.5)) # clumping and orient factor
  }
  name_vec <- rep(c("N","Cab","Car","Cw","Cm","leaf_mass_per_area","clumping_factor","orient_factor"),length(pft_end))
  samples[which(name_vec=="leaf_mass_per_area")] <- 1/samples[which(name_vec=="leaf_mass_per_area")]*1000 # convert to SLA
  names(samples) <- rep(c('N', 'Cab', 'Car', 'Cw', 'Cm', 'SLA', 'clumping_factor', 'orient_factor'),3)
  return(samples)
}

prior_bt <- BayesianTools::createPrior(density = prior_function, sampler = inits_function)

# Test observation param values
obs_params <- function() {
         #N, Cab, Car, Cw, Cm, SLA, clumping, orient
vals <- c(1.8, 47, 8.7, 0.009, 0.007, (1/66.3)*1000, 0.8, 0.12,     # Early
    1.4, 47, 8.8, 0.01, 0.009, (1/128.3)*1000, 0.82, 0.12,     # Mid
    1.9, 45, 8.5, 0.007, 0.008, (1/65.35)*1000, 0.86, 0.12)    # Late
names(vals) <- rep(c('N', 'Cab', 'Car', 'Cw', 'Cm', 'SLA', 'clumping_factor', 'orient_factor'),3)
return(vals)
}

vec2list <- function(params, ...) {
    spectra_list <- list()
    ed_list <- list()
    for (i in seq_along(pft_end)) {
        pft <- names(pft_end)[i]
        param_sub <- param_sub(i, params)
        spectra_list[[pft]] <- PEcAnRTM::prospect(param_sub[1:5], 5, TRUE)
        ed_list[[pft]] <- list(SLA = param_sub[6],
                               clumping_factor = param_sub[7],
                               orient_factor = param_sub[8])
    }
    outlist <- list(spectra_list = spectra_list, trait.values = ed_list, ...)
    return(outlist)
}

# Quick test to make sure everything works
testrun <- inits_function() %>%
    vec2list(datetime = datetime, par.wl = 400:2499, nir.wl = 2500) %>%
    run_edr(prefix, edr_args = .)

head(testrun)

model <- function(params) {
    edr_dir <- 'edr'
    args_list <- vec2list(params,
                          paths = list(ed2in = NA, history = file.path(prefix, 'outputs')),
                          par.wl = 400:2499,
                          nir.wl = 2500,
                          datetime = datetime,
                          change.history.time = FALSE)
    albedo <- run_edr(prefix, args_list, edr_dir)

    if (plot_albedo) {
      # Create quick figure
      waves <- seq(400,2500,1)
      png(file.path(prefix,edr_dir,'simulated_albedo.png'),width=4900, height =2700,res=400)
      par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
      plot(waves,unlist(albedo)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",
           ylab="Reflectance (%)",
      cex.axis=1.5, cex.lab=1.7,col="black")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey80")
      lines(waves,unlist(albedo)*100,lwd=3, col="black")
      dev.off()
    }

    return(albedo)
}

# Test that the model works
test_model <- model(inits_function())
head(test_model)

# Simulate observations
# Add small amount of noise so first fit isn't perfect (this breaks neff calculation)
observed <- model(obs_params()) + PEcAnRTM::generate.noise()
#--------------------------------------------------------------------------------------------------#

invert_options <- list(
  init = list(iterations = 200),
  loop = list(iterations = 100),
  other = list(max_iter = 1e6,
               save_progress = file.path(prefix, "inversion_samples_inprogress.rds")))

#--------------------------------------------------------------------------------------------------#
samples <- PEcAnRTM::invert_bt(observed = observed, model = model, prior = prior_bt,
                               custom_settings = invert_options)
save(samples, file = file.path(prefix,'inversion_samples_finished.RData'))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (generate_summary_figs) {
  
  main_out <- prefix
  samples_mcmc <- BayesianTools::getSample(samples, coda = TRUE)
  
  coda::niter(samples_mcmc)
  coda::nvar(samples_mcmc)
  coda::nchain(samples_mcmc)
  
  
  par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  png(file.path(main_out,"final_trace_plot.png"), width = 1500, height = 1600, res=150)
  plot(samples_mcmc)
  dev.off()
  
  rawsamps <- do.call(rbind, samples_mcmc)
  par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  png(file.path(main_out,"final_pairs_plot.png"), width = 1500, height = 1600, res=150)
  pairs(rawsamps)
  dev.off()
  
  # par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  # png(file.path(main_out,"final_deviance_plot.png"), width = 1500, height = 1600, res=150)
  # plot(PEcAn.assim.batch::makeMCMCList(input.pda.data$deviance))
  # dev.off()
  # 
  # par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  # png(file.path(main_out,"finale_neff_plot.png"), width = 1500, height = 1600, res=150)
  # plot(PEcAn.assim.batch::makeMCMCList(input.pda.data$n_eff_list))
  # dev.off()
}
#--------------------------------------------------------------------------------------------------#
### EOF