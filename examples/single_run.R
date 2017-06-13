#--------------------------------------------------------------------------------------------------#
source('config.R')

library(redr)

save_plot <- FALSE #TRUE/FALSE
hidden <- TRUE #TRUE/FALSE
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# This is the directory for storing all required elements for a single run.
# It provides a shortcut for executing related instances of ED2 and EDR.
# NOTE: Any directories starting with `.edr` will be automatically gitignored
if (hidden) {
  prefix <- '.edr-single-run'
} else {
  prefix <- 'edr-single-run'
}

# Load default css, pss, and site files
data(css_ex1)   # css_df
data(pss_ex1)   # pss_df
data(site_ex1)  # site_df

# Generate ED2 run configuration
genrun <- generate_run(prefix = prefix,
                       # Site coordinates
                       site_lat = site_lat, 
                       site_lon = site_lon,
                       # CSS, PSS, and SITE files
                       # NOTE: These can be valid data.frames OR paths to files
                       site_df = site_df, 
                       pss_df = pss_df, 
                       css_df = css_df,
                       # Other required paths (set in config.R)
                       common_inputs_dir = common_inputs_dir, 
                       site_met_dir = site_met_dir,
                       ed_exe_path = ed_exe_path,
                       # Purge prefix directory if already exists (useful for debugging)
                       RMDIR = TRUE)

# Run ED2 using configuration 
# This will generate the output files necessary for EDR
# NOTE: ED output is stored in return object as a character
message('Running ED2...')
run_log <- run_ed(prefix)
print(tail(run_log))
message('Done!')

# Setup EDR in a special subdirectory of the ED2 directory used above
edr_setup <- setup_edr(prefix, edr_exe_path = edr_exe_path)

# Place EDR parameters into a data.frame, and use a custom function to convert 
# to list to be passed into EDR.
# The `tibble::tribble` function is a convenient way to quickly create simple data frames.
edr_params_df <- tibble::tribble(
    ~pft, ~N, ~Cab, ~Cw, ~Cm, ~orient_factor, ~clumping_factor,
    'temperate.Early_Hardwood', 1.4, 40, 0.01, 0.01, 0.5, 0.5, 
    'temperate.North_Mid_Hardwood', 1.4, 30, 0.01, 0.01, 0.5, 0.5,
    'temperate.Late_Hardwood', 1.4, 20, 0.01, 0.01, 0.5, 0.5
    )
edr_params_list <- params_df2list(edr_params_df, prospect_version = 4, pftcol = "pft", datetime = datetime)
edr_params_list$par.wl <- 400:2499
edr_params_list$nir.wl <- 2500
str(edr_params_list)


# Alternatively, can just pass a list in directly.
# The list needs to contain, at a minimum:
#   - `$spectra_list` -- Leaf spectra for each PFT, as a named list, with names 
#   corresponding to PFTs 
#   - `$trait.values` -- EDR trait values, as a named list, like above. 
#       - NOTE: Each PFT in spectra_list also has to exist in `trait.values`, 
#       even if no parameters are passed (just use a blank list for that PFT)
#   - `$datetime` -- POSIXlt object containing the date and time of execution

# Run EDR
albedo <- run_edr(prefix, edr_params_list)
head(albedo)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (save_plot) { 
  waves <- seq(400,2500,1)
  png(file.path(prefix,'simulated_albedo.png'),width=4900, height =2700,res=400)
  par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
  plot(waves,unlist(albedo)*100,type="l",lwd=3,ylim=c(0,65),xlab="Wavelength (nm)",ylab="Reflectance (%)",
       cex.axis=1.5, cex.lab=1.7,col="black")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
         "grey80")
  lines(waves,unlist(albedo)*100,lwd=5, col="black")
  dev.off()
}

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF