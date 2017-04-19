library(redr)
edr_da_dir <- '/home/ashiklom/Projects/nasa-rtm/edr-da'
ed_exe_path <- '/home/ashiklom/Projects/ED2/ED/build/ed_2.1'
#devtools::load_all(file.path(edr_da_dir, 'redr'))

test_outdir <- 'edr-testthat-outdir'

css_df <- data.frame(year = 2000, patch = 1, cohort = 1, dbh = 20, ht = 0, pft = 9,
                     den = 0.05, bdead = 0, balive = 0, lai = -999)

data(pss_ex1)
data(site_ex1)

ed2in_changes <- list(IMONTHA = 06, IDATEA = 29, IYEARA = 2004,
                      IMONTHZ = 07, IDATEZ = 04, IYEARZ = 2004)

try_genrun <- generate_run(prefix = test_outdir,
                           site_lat = 45.5,
                           site_lon = -90.5,
                           css_df = css_df, 
                           pss_df = pss_df, 
                           site_df = site_df,
                           output_dir = test_outdir,
                           common_inputs_dir = file.path(edr_da_dir, 'ed-inputs/EDI'),
                           site_met_dir = file.path(edr_da_dir, 'ed-inputs/met3/US-WCr/'),
                           ed_exe_path = ed_exe_path,
                           ed2in_template = file.path(edr_da_dir, 'run-ed/template/ED2IN'),
                           ed2in_changes = ed2in_changes,
                           RMDIR = TRUE)

test_run <- run_ed(test_outdir)

# Test EDR
dir <- 'edr-testthat-outdir'
edr_exe_path <- '~/Projects/ED2/EDR/build/ed_2.1'
edr_args <- list(spectra_list = list(temperate.Early_Hardwood = 
                                     PEcAnRTM::prospect(c(1.4, 30, 0.01, 0.01), 4, TRUE)),
                 trait.values = list(temperate.Early_Hardwood = list()),
                 datetime = as.POSIXlt('2004-06-29 12:00:00'))

setup_edr(dir, edr_exe_path)
albedo <- run_edr(dir, edr_args)
