# Figure dir abbreviations:
# h - Heteroskedastic variance (variance ~ a * reflectance + b)
# ss - Site-specific variance (each site has a separate variance calculation)

figdir <- dir_create(here("text", "figures", "h-ss"))
site_list_file <- here("other_site_data", "site_list")
site_structure_file <- here("other_site_data", "site_structure.csv")
fft_lai_file <- here("other_site_data", "NASA_FFT_LAI_FPAR_Data.csv")

pda_result_file <- last_result_file("multi_site_pda_results-exp")
pda_start <- 6000
param_names_file <- path(path_dir(pda_result_file), "param_names.txt")
