library(dbhelpers)
library(dplyr)
library(tibble)

source('common.R')

land_data <- tribble(
  ~name,                    ~format_id,   ~file_path,                                 ~file_name,
  'EDR_testing.glu',        41,           normalizePath('ed-inputs'),                 'glu',
  'EDR_testing.chd_dgd',    42,           normalizePath('ed-inputs/EDI/ed_inputs'),   '',
  'EDR_testing.oge2OLD',    43,           normalizePath('ed-inputs/EDI/oge2OLD'),     'OGE2_',
  'EDR_testing.fao',        44,           normalizePath('ed-inputs'),                 'fao') %>%
  mutate(notes = 'EDR_testing',
         site_id = site_id,
         machine_id = machine_id,
         container_type = 'Input')

land_input <- db_merge_into(db = db,
                           table = 'inputs',
                           values = land_data,
                           by = 'name',
                           id_colname = 'id')

land_data <- left_join(land_data, select_(land_input, container_id = 'id', 'name'))

land_dbfiles <- db_merge_into(db = db,
                             table = 'dbfiles',
                             values = land_data,
                             by = c('file_name', 'file_path'),
                             id_colname = 'id')

