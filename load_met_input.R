library(dbhelpers)
library(dplyr)

source('common.R')

met_driver <- tibble(format_id = 12,   # ED_MET_DRIVER_HEADER
                     file_path = normalizePath('ed-inputs/met3/US-WCr'),
                     file_name = 'ED_MET_DRIVER_HEADER',
                     notes = 'EDR_testing.ED_MET_DRIVER_HEADER',
                     site_id = site_id,
                     machine_id = machine_id,
                     container_type = 'Input')

met_driver$name <- met_driver$notes

met_input <- db_merge_into(db = db,
                           table = 'inputs',
                           values = met_driver,
                           by = 'name',
                           id_colname = 'id')

met_driver <- left_join(met_driver, select_(met_input, container_id = 'id', 'name'))

met_dbfiles <- db_merge_into(db = db,
                             table = 'dbfiles',
                             values = met_driver,
                             by = c('file_name', 'file_path'),
                             id_colname = 'id')

