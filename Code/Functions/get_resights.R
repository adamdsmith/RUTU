get_resights <- function(rutu_dep = rutu_deploy) {
  fn <- list.files("Data", "Smith_all_RUTU", full = TRUE)
  bbl <- lapply(fn, function(i) {
    tmp <- readr::read_csv(i) %>%
      janitor::clean_names() %>%
      select(band_num, lat = e_lat_decimal_degrees,
             lon = e_lon_decimal_degrees,
             rdate = encounter_date) %>%
      mutate(rdate = as_date(mdy(rdate, tz = my_tz)),
             source = "BBL")
  })
  bborg <- readxl::read_excel("Data/RUTU encounters bandedbirds SC April 2020.xlsx") %>%
    select(band_num = MetalID, lat = Resighting.Latitude,
           lon = Resighting.Longitude, rdate = ResightDate)  %>%
    mutate(band_num = as.numeric(gsub("-", "", band_num)),
           rdate = as_date(rdate, tz = my_tz),
           source = "bandedbird.org")
  resight <- bind_rows(bbl, bborg) %>%
    mutate(lat = round(lat, 3),
           lon = round(lon, 3)) %>%
    distinct() %>%
    inner_join(rutu_dep, by = "band_num") %>%
    mutate(dtime = time_length(rdate - bdate, "year")) %>%
    filter(dtime < 0.5)
  return(resight)
}