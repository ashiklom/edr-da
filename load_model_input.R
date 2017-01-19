library(dbhelpers)
library(dplyr)

source('common.R')

tbl(db, 'models') %>% glimpse()
tbl(db, 'dbfiles') %>% glimpse()

model_id <- tbl(db, 'models') %>%
    filter(model_name %like% '%modularized%') %>%
    select(id) %>%
    collect() %>% .[['id']]

model <- tibble(file_name = 'ed-2.1',
                file_path = normalizePath('~/Projects/ED2/ED/build'),
                machine_id = machine_id,
                container_type = 'Model',
                container_id = model_id)

model_dbfiles <- db_merge_into(db = db,
                               table = 'dbfiles',
                               values = model,
                               by = c('file_name', 'file_path'),
                               id_colname = 'id')
