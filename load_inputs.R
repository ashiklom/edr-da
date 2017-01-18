library(PEcAn.DB)
library(dplyr)
library(dbhelpers)

db <- src_postgres(dbname = 'bety', user = 'bety', password = 'bety')

# Create inputs
site_id <- 676          # Willow Creek (US-WCr)
format_id <- 1073       # ED site/css/pss file
machine_id <- 77        # My laptop machine ID
file_path <- normalizePath('ed-inputs/sites/US-WCr/rtm/1cohort/dbh20/early_hardwood')
met_path <- normalizePath('ed-inputs/met3/US-WCr')

dbGetQuery(db$con,
'SELECT id from dbfiles limit 1
'
)

new_inputs_vec <- c('EDR_test.1cohort.cohort',
                    'EDR_test.1cohort.patch',
                    'EDR_test.1cohort.site',
                    'EDR_test.met_driver')

file_prefix <- 'early_hardwood.lat45.5lon-90.5'

new_inputs <- tibble(name = new_inputs_vec) %>%
  mutate(site_id = site_id,
         format_id = case_when(grepl('\\.site$', .$name) ~ 10,
                               grepl('\\.patch$', .$name) ~ 15,
                               grepl('\\.cohort$', .$name) ~ 11,
                               grepl('\\.met_driver$', .$name) ~ 12),
         container_type = 'Input',
         machine_id = 77) %>%
  mutate(file_path = recode(format_id, `12` = met_path, .default = file_path),
         file_name = recode(format_id, 
                            `10` = paste0(file_prefix, '.site'),
                            `15` = paste0(file_prefix, '.pss'),
                            `11` = paste0(file_prefix, '.css'),
                            `12` = 'ED_MET_DRIVER_HEADER'))

inputs <- db_merge_into(db = db,
                        table = 'inputs',
                        values = new_inputs,
                        by = 'name',
                        id_colname = 'id') 

dbfiles <- db_merge_into(db = db,
                         table = 'dbfiles',
                         values = new_inputs,
                         by = c('file_name', 'file_path'),
                         id_colname = 'id')

#insert <- insert_dbfile(file_path = file1_path,
                        #file_prefix = file1_prefix,
                        #siteid = siteid,
                        #startdate = startdate,
                        #enddate = enddate,
                        #machineid = machineid,
                        #db = db,
                        #hostname = hostname)
