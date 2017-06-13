#--------------------------------------------------------------------------------------------------#
source('config.R')

library(magrittr)
library(redr)

data(css_ex1)
data(pss_ex1)
data(site_ex1)

plot_albedo <- FALSE #TRUE/FALSE
generate_summary_figs <- TRUE #TRUE/FALSE
hidden <- TRUE  #TRUE/FALSE

nchains <- 3 #3
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (hidden) {
  prefix <- '.edr_inversion'
} else {
  prefix <- paste("edr_inversion", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
  #prefix <- 'edr_inversion'
}
PEcAn.utils::logger.info(paste0("Running inversion in dir: ",prefix))

css_bl_same_dbh <- extend_df(css_df, cohort = 1, dbh = 30, pft = c(11))
PEcAn.utils::logger.info(css_bl_same_dbh)
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
# temperate.Early_Hardwood: 1:8
#   1:5 -- PROSPECT params
#   6 -- sla
#   7 -- clumping_factor
#   8 -- orient_factor
# temperate.North_Mid_Hardwood: 9:16
# temperate.Late_Hardwood: 17:24
pft_end <- c(temperate.Late_Hardwood = 8)
PEcAn.utils::logger.info(" pft_end ")
PEcAn.utils::logger.info(pft_end)
PEcAn.utils::logger.info(names(pft_end))

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

        #names(prior) <- c('N', 'Cab', 'Car', 'Cw', 'Cm', 'SLA', 'clumping_factor', 'orient_factor')
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
inits_function <- function() {
                  # N, Cab, Car, Cw, Cm, SLA, clumping, orient
vals <- rnorm(8, c(2, 35, 5, 0.006, 0.005, 15, 0.5, 0), c(0.2, 2, 1, 0.001, 0.001, 3, 0.1, 0.1)) # Late
names(vals) <- rep(c('N', 'Cab', 'Car', 'Cw', 'Cm', 'SLA', 'clumping_factor', 'orient_factor'),1)
return(vals)
}

# Test observation param values
obs_params <- function() {
         #N, Cab, Car, Cw, Cm, SLA, clumping, orient
vals <- c(1.9, 45, 8.5, 0.007, 0.008, (1/65.35)*1000, 0.86, 0.12)    # Late
names(vals) <- rep(c('N', 'Cab', 'Car', 'Cw', 'Cm', 'SLA', 'clumping_factor', 'orient_factor'),1)
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

run_first <- function(inputs) {
    edr_dir <- paste('edr', inputs$runID, sep = '.')

    # Setup chain-specific directory for EDR
    setup_edr(prefix, edr_exe_path, edr_dir)

    # Create dummy params from initial conditions
    args_list <- vec2list(inits_function(), datetime = datetime)

    # Run EDR
    albedo <- run_edr(prefix, args_list, edr_dir)
    return(albedo)
}

model <- function(params, runID = 'test') {
    edr_dir <- paste('edr', runID, sep = '.')
    args_list <- vec2list(params, 
                          paths = list(ed2in = NA, history = file.path(prefix, 'outputs')),
                          par.wl = 400:2499,
                          nir.wl = 2500,
                          datetime = datetime, 
                          change.history.time = FALSE)
    albedo <- run_edr(prefix, args_list, edr_dir)
    print(head(albedo))
    
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
test_runfirst <- run_first(list(runID = 'test'))
test_model <- model(inits_function(), runID = 'test')
head(test_model)

# Other inversion parameters
#param_mins <- rep(c(1, rep(0, 7)), 3)
param_mins <- rep(c(N = 1, Cab = 1, Car = 0, Cw = 0.0001, Cm = 0.0001, SLA = 1, clumping_factor = 0.001, 
		orient_factor = -0.5),1)

invert_options <- list(model = model,
                       run_first = run_first,
                       nchains = nchains,
                       inits.function = inits_function,
                       prior.function = prior_function,
                       ngibbs.max = 1e6,
                       ngibbs.min = 500,
                       ngibbs.step = 1000,
                       param.mins = param_mins)

# Simulate observations
# Add small amount of noise so first fit isn't perfect (this breaks neff calculation)
observation <- run_first(list(runID = 'observation'))
observed <- model(obs_params(), runID = 'observation') + PEcAnRTM::generate.noise()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
logfile <- "invert.auto_log.txt"
samples <- PEcAnRTM::invert.auto(observed = observed,
                                 invert.options = invert_options,
                                 parallel = TRUE,
                                 parallel.output = file.path(prefix,logfile),
                                 save.samples = file.path(prefix,'inversion_samples_inprogress.rds'))
save(samples, file = file.path(prefix,'inversion_samples_finished.RData'))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (generate_summary_figs) {

  main_out <- prefix
  samples.bt <- PEcAn.assim.batch::autoburnin(samples$samples)
  samples.bt <- PEcAn.assim.batch::makeMCMCList(samples.bt)
  
  par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  png(file.path(main_out,"final_trace_plot.png"), width = 1500, height = 1600, res=150)
  plot(samples.bt)
  dev.off()
  
  rawsamps <- do.call(rbind, samples.bt)
  par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  png(file.path(main_out,"final_pairs_plot.png"), width = 1500, height = 1600, res=150)
  pairs(rawsamps)
  dev.off()
  
  par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  png(file.path(main_out,"final_deviance_plot.png"), width = 1500, height = 1600, res=150)
  plot(PEcAn.assim.batch::makeMCMCList(input.pda.data$deviance))
  dev.off()
  
  par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  png(file.path(main_out,"finale_neff_plot.png"), width = 1500, height = 1600, res=150)
  plot(PEcAn.assim.batch::makeMCMCList(input.pda.data$n_eff_list))
  dev.off()
}
#--------------------------------------------------------------------------------------------------#
### EOF
