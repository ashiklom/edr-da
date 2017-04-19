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

css_case1 <- extend_df(css_df, cohort = 1:3, dbh = c(20, 30, 40), pft = c(9, 10, 11))
genrun <- generate_run(prefix = 'case1',
                       site_lat = site_lat,
                       site_lon = site_lon,
                       site_df = site_df,
                       pss_df = pss_df,
                       css_df = css_case1,
                       common_inputs_dir = common_inputs_dir,
                       site_met_dir = site_met_dir,
                       ed_exe_path = ed_exe_path,
                       RMDIR = TRUE)

case1_run <- run_ed('case1')

case1_setup <- setup_edr('case1', edr_exe_path = edr_exe_path)
case1_edr_args <- list(spectra_list = 
                       list('temperate.Early_Hardwood' = 
                            PEcAnRTM::prospect(c(1.4, 20, 0.01, 0.01), 4, TRUE)),
                       trait.values = list('temperate.Early_Hardwood' = list()),
                       datetime = datetime)
case1_albedo <- run_edr('case1', edr_args = case1_edr_args)

