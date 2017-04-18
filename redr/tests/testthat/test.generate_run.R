library(testthat)
#library(redr)
edr_da_dir <- '/home/ashiklom/Projects/nasa-rtm/edr-da'
ed_exe_path <- '/home/ashiklom/Projects/ED2/ED/build/ed_2.1'
devtools::load_all(file.path(edr_da_dir, 'redr'))

test_outdir <- 'edr-testthat-outdir'

context('Check that a working test run is created')

css_df <- data.frame(year = 2000, patch = 1, cohort = 1, dbh = 20, ht = 0, pft = 9,
                     den = 0.05, bdead = 0, balive = 0, lai = -999)

pss_df <- data.frame(site = 1, year = 2000, patch = 1, dst = 1, age = 70, area = 1,
                     water = 0.1, fsc = 5, stsc = 2.15, stsl = 2.15, ssc = 0.03,
                     psc = 0, msn = 0.16, fsn = 1.14)

site_df <- data.frame(sitenum = 1, area = 1.0, TCI = -6.3, elev = 540.0, slope = 0,
                      aspect = 0, soil1 = 3)

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


#ed2_text <- readLines(ed2in_template)
#tag <- 'VEG_DATABASE'
#regex <- sprintf("(^[[:blank:]]*NL%%%s)[[:blank:]]+=.*", tag)
#grep(regex, ed2_text, value = TRUE)
