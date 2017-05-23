############################################################
# EDR data assimilation
# Configuration script
#
# Copy this script to 'config.R' and `source` it at the 
# beginning of any EDR script
############################################################

############################################################
# Set paths
############################################################
# NOTE: `normalizePath` is called to generate a complete 
# path, which should make this script working 
# directory-agnostic.

# This directory should contain the following:
#   - Degree-day inputs -- ed_inputs/{chd,dgd}
#   - OGE data -- oge2OLD/{*.h5, OGE2_HEADER}
common_inputs_dir <- file.path('ed-inputs', 'EDI')
common_inputs_dir <- normalizePath(common_inputs_dir)

# This directory should contain COMPLETE site-specific met data:
#   - Met driver header (complete, NOT a temp file) -- ED_MET_DRIVER_HEADER
#   - Met data -- <site-tag>_2004{JUN,JUL,AUG}.h5
site_met_dir <- file.path('ed-inputs', 'met3', 'US-WCr')
common_inputs_dir <- normalizePath(site_met_dir)

# These are the paths to the executables for the full ED2 
# model and EDR, respectively
ed_exe_path <- '/home/ashiklom/Projects/ED2/ED/build/ed_2.1'
edr_exe_path <- '/home/ashiklom/Projects/ED2/EDR/build/ed_2.1'

############################################################
# Additional configuration
############################################################

# Site coordinates
site_lat <- 45.5
site_lon <- -90.5

# Date-time for EDR runs
datetime <- as.POSIXlt('2004-07-01 12:00:00')
