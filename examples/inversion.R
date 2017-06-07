source('config.R')

library(magrittr)
library(redr)

data(css_ex1)
data(pss_ex1)
data(site_ex1)

prefix <- '.edr_inversion'

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
        prior <- prior + mvtnorm::dmvnorm(c(param_sub[1:5], 1/param_sub[6]), means[pft,], covars[pft,,], 
                                          log = TRUE)
        # ED priors
        prior <- prior + dunif(param_sub[7], 0, 1, TRUE) + dunif(param_sub[8], -1, 1, TRUE)
    }
    return(prior)
}

# Static initial conditions
# Can be replaced with a function that draws these values from distributions
inits_function <- function() {
    # N, Cab, Car, Cw, Cm, SLA, clumping, orient
    c(1.8, 47, 8.7, 0.009, 0.007, 1/66.3, 0.5, 0,     # Early
      1.4, 47, 8.8, 0.01, 0.009, 1/128.3, 0.5, 0,     # Mid
      1.9, 45, 8.5, 0.007, 0.008, 1/65.35, 0.5, 0)    # Late
}

vec2list <- function(params, ...) {
    spectra_list <- list()
    ed_list <- list()
    for (i in seq_along(pft_end)) {
        pft <- names(pft_end)[i]
        param_sub <- param_sub(i, params)
        spectra_list[[pft]] <- PEcAnRTM::prospect(param_sub[1:5], 5, TRUE)
        ed_list[[pft]] <- list(sla = param_sub[6],
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
    return(albedo)
}

# Test that the model works
test_runfirst <- run_first(list(runID = 'test'))
test_model <- model(inits_function(), runID = 'test')
head(test_model)

# Other inversion parameters
param_mins <- rep(c(1, rep(0, 7)), 3)

invert_options <- list(model = model,
                       run_first = run_first,
                       nchains = 3,
                       inits.function = inits_function,
                       prior.function = prior_function,
                       ngibbs.max = 1e6,
                       ngibbs.min = 500,
                       ngibbs.step = 1000,
                       param.mins = param_mins)

# Simulate observations
# Add small amount of noise so first fit isn't perfect (this breaks neff calculation)
observed <- model(inits_function(), runID = 'test') + PEcAnRTM::generate.noise()

samples <- PEcAnRTM::invert.auto(observed = observed,
                                 invert.options = invert_options,
                                 parallel = TRUE,
                                 save.samples = 'inversion_samples_inprogress.rds')
save(samples, file = 'inversion_samples_finished.RData')
