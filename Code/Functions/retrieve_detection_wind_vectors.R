retrieve_detection_wind_vectors <- function(det_dat, plevel = 1000) {
  if (!requireNamespace("RNCEP", quietly = TRUE)) install.packages("RNCEP", quiet = TRUE)
  if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate", quiet = TRUE)
  lat_rng <- c(floor(min(det_dat$recvDeployLat)), ceiling(max(det_dat$recvDeployLat)))
  lon_rng <- c(floor(min(det_dat$recvDeployLon)), ceiling(max(det_dat$recvDeployLon)))
  years_rng <- range(lubridate::year(det_dat$first_det))
  months_rng <- range(lubridate::month(det_dat$first_det))
  used_m_y <- tibble(year = lubridate::year(det_dat$first_det),
                     month = lubridate::month(det_dat$first_det)) %>%
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
  speed_dat_grid <- pbapply::pblapply(seq_len(nrow(det_dat)), function(i) {
    i_xy <- as.matrix(det_dat[i, c("recvDeployLon", "recvDeployLat")])
    i_z <- pull(det_dat[i, ], first_det)
    grid_dists <- geosphere::distGeo(i_xy, as.matrix(grid_xy))
    grid_tdiffs <- as.numeric(abs(difftime(i_z, grid_z$datetime, units = "hours")))
    out <- bind_cols(grid_xy[which.min(grid_dists), ],
                     grid_z[which.min(grid_tdiffs), , drop = FALSE])
    out <- inner_join(wind_df, out, by = c("datetime", "latitude", "longitude"))
    stopifnot(nrow(out) == 1)
    bind_cols(det_dat[i, ], out[, c("uwind", "vwind")])
  })
  det_dat <- bind_rows(speed_dat_grid)

  return(det_dat)
}
