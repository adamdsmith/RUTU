retrieve_path_wind_components <- function(speed_dat, plevel = 1000,
                                          tailwind = TRUE) {
  if (!requireNamespace("RNCEP", quietly = TRUE)) install.packages("RNCEP", quiet = TRUE)
  if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate", quiet = TRUE)
  lat_rng <- c(floor(min(c(speed_dat$start_recvDeployLat, speed_dat$end_recvDeployLat))),
                  ceiling(max(c(speed_dat$start_recvDeployLat, speed_dat$end_recvDeployLat))))
  lon_rng <- c(floor(min(c(speed_dat$start_recvDeployLon, speed_dat$end_recvDeployLon))),
                  ceiling(max(c(speed_dat$start_recvDeployLon, speed_dat$end_recvDeployLon))))
  years_rng <- range(lubridate::year(speed_dat$start_time))
  months_rng <- range(lubridate::month(speed_dat$start_time))
  used_m_y <- tibble(year = lubridate::year(speed_dat$start_time),
                     month = lubridate::month(speed_dat$start_time)) %>%
    distinct()
  message("Downloading relevant NCEP Reanalysis II grids...")
  uwind <- RNCEP::NCEP.gather("uwnd", level = plevel,  months.minmax = months_rng,
                              years.minmax = years_rng, lat.southnort = lat_rng,
                              lon.westeast = lon_rng, reanalysis2 = TRUE, 
                              status.bar = FALSE, return.units = FALSE)
  vwind <- RNCEP::NCEP.gather("vwnd", level = plevel,  months.minmax = months_rng,
                              years.minmax = years_rng, lat.southnort = lat_rng,
                              lon.westeast = lon_rng, reanalysis2 = TRUE, 
                              status.bar = FALSE, return.units = FALSE)
  wind_df <- RNCEP::NCEP.array2df(wx.data = list(uwind, vwind),
                                  var.names = c("uwind", "vwind")) %>%
    mutate(datetime = lubridate::ymd_h(datetime),
           longitude = longitude - 360, # Won't work if crosses PM or IDL
           year = lubridate::year(datetime),
           month = lubridate::month(datetime)) %>%
    inner_join(used_m_y, by = c("year", "month"))

  grid_xy <- select(wind_df, longitude, latitude) %>% distinct()
  grid_z <- select(wind_df, datetime) %>% distinct()

  message("Appending u and v wind components (m/s)...")
  speed_dat_grid <- pbapply::pblapply(seq_len(nrow(speed_dat)), function(i) {
    i_xy <- as.matrix(speed_dat[i, c("start_recvDeployLon", "start_recvDeployLat")])
    i_z <- pull(speed_dat[i, ], start_time)
    grid_dists <- geosphere::distGeo(i_xy, as.matrix(grid_xy))
    grid_tdiffs <- as.numeric(abs(difftime(i_z, grid_z$datetime, units = "hours")))
    out <- bind_cols(grid_xy[which.min(grid_dists), ],
                     grid_z[which.min(grid_tdiffs), , drop = FALSE])
    out <- inner_join(wind_df, out, by = c("datetime", "latitude", "longitude"))
    stopifnot(nrow(out) == 1)
    bind_cols(speed_dat[i, ], out[, c("uwind", "vwind")])
  })
  speed_dat <- bind_rows(speed_dat_grid)
  
  # Add tailwind, if requested
  if (tailwind) {
    message("Appending tailwind estimates...")
    speed_sub <- speed_dat %>%
      select(motusTagID, start_recvDeployID, end_recvDeployID, brng_init, uwind, vwind) %>%
      rowwise() %>%
      mutate(calcs = purrr::map(uwind, RNCEP::NCEP.Tailwind, vwind, brng_init)) %>%
      tidyr::unnest(calcs) %>%
      select(motusTagID:end_recvDeployID, tailwind)
    if (!tailwind) speed_sub <- select(speed_sub, -tailwind)
    speed_dat <- left_join(speed_dat, speed_sub, 
                           by = c("motusTagID", "start_recvDeployID", "end_recvDeployID"))
  }
  return(speed_dat)
}
