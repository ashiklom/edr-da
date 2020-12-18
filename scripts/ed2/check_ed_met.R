library(hdf5r)
library(lubridate)

ed_metfiles_inorder <- function(ed_met_dir) {
  file_list <- list.files(ed_met_dir, "\\.h5$")
  frxp <- ".*([[:digit:]]{4})([[:upper:]]{3})\\.h5"
  years <- as.numeric(gsub(frxp, "\\1", file_list))
  month_levels <- toupper(substr(month.name, 0, 3))
  month_names <- factor(gsub(frxp, "\\2", file_list), month_levels)
  dates <- paste(years, as.integer(month_names), "01", sep = "-") %>%
    lubridate::as_date()
  tibble(
    file_base = file_list,
    file_full = file.path(ed_met_dir, file_list),
    dates = dates
  ) %>%
    arrange(dates)
}

h5tibble <- function(hfile, pb = NULL) {
  stopifnot(file.exists(hfile))
  on.exit(if (!is.null(pb)) pb$tick())
  h <- H5File$new(hfile, "r")
  on.exit(h$close_all(), add = TRUE)
  variables <- names(h)
  map(variables, ~h[[.]][,,]) %>%
    setNames(variables) %>%
    dplyr::bind_cols()
}

read_metfiles <- function(ed_met_dir) {
  files_df <- ed_metfiles_inorder(ed_met_dir)
  pb <- progress::progress_bar$new(total = nrow(files_df))
  files_df2 <- files_df %>%
    mutate(data = map(file_full, h5tibble, pb = pb))
  files_df3 <- files_df2 %>%
    mutate(
      data = map(data, ~mutate(., i = seq_len(nrow(.)))),
      data = map2(data, dates, ~mutate(.x, datetime = .y + (i - 1) * lubridate::dhours(3)))
    )
  tidyr::unnest(files_df3, data)
}

ed_met_dir <- "sites/SF03_site_1-25721/ED_NARR"
full_data <- read_metfiles(ed_met_dir)

sub_data <- full_data %>%
  filter(datetime >= as.Date("2009-07-01"), datetime < as.Date("2009-08-01")) %>%
  mutate(d = mday(datetime), h = hour(datetime))

xint <- as.POSIXct(c("2009-07-18", "2009-07-19"))
i <- imguR::imgur("pdf")
p <- ggplot(sub_data) +
  aes(x = h, y = dlwrf, group = d) +
  geom_point()
  #geom_vline(xintercept = xint, color = "red", linetype = "dashed")
p
p %+% aes(y = nbdsf)
p %+% aes(y = nddsf)
p %+% aes(y = prate)
p %+% aes(y = pres)
p %+% aes(y = sh)
p %+% aes(y = tmp)
imguR::imgur_off(i)

long_data <- sub_data %>%
  select(datetime, dlwrf:vgrd) %>%
  gather(variable, value, -datetime)

ggplot(long_data) +
  aes(x = datetime, y = value) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y")
