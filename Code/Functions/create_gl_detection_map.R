create_gl_detection_map <- function(gl_det_data, gl_sf, legend_pos = c(0.708, 0.9), 
                                    xlim = c(-92.5, -73.5), ylim = c(40, 50.5),
                                    det_color = "#1f78b4", stop_color = "#8adf90",
                                    stop_buffer = 15, map_crs = NULL) {
  
  if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc", quiet = TRUE)
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) install.packages("rnaturalearth", quiet = TRUE)
  if (!requireNamespace("ggspatial", quietly = TRUE)) install.packages("ggspatial", quiet = TRUE)
  
  dets <- select(gl_det_data, motusTagID, recvDeployID, recvSiteName, lat = recvDeployLat, lon = recvDeployLon, dt = first_det)
  
  gl_stops_solo <- filter(gl_det_data, stop_solo) %>%
    st_as_sf(coords = c("recvDeployLon", "recvDeployLat"), crs = 4326)
  gl_stops_adj <- filter(gl_det_data, stop_adj) %>%
    st_as_sf(coords = c("recvDeployLon", "recvDeployLat"), crs = 4326)
  gl_all_stops <- rbind(gl_stops_solo, gl_stops_adj)
  gl_all_stops <- buffer_ll_pts(gl_all_stops, stop_buffer)
  n_stop_recvs <- unlist(strsplit(paste(c(gl_all_stops$recvDeployID, na.omit(gl_all_stops$adj_recv_grp)), collapse = ","), ",")) %>%
    n_distinct()
  message("Stopovers of ", n_distinct(gl_all_stops$motusTagID), " individuals will be dispayed ",
          "representing ", n_stop_recvs, " receivers.")

  station_n <- gl_det_data %>%
    group_by(recvSiteName) %>%
    summarize(n_indivs = n_distinct(motusTagID), .groups = "drop")
  
  gl_station_locs <- select(gl_det_data, recvSiteName, lat = recvDeployLat, lon = recvDeployLon) %>%
    distinct() %>%
    group_by(recvSiteName) %>%
    # Keep only one location in case of slight variation in position among deployments
    slice(1) %>%
    ungroup()
  
  station_n <- left_join(station_n, gl_station_locs, by = "recvSiteName") %>%
    mutate(pct_indivs = round(n_indivs / n_distinct(dets$motusTagID) * 100, 1))
  
  # Gather active receivers in network
  recvs <- get_motus_recvs(recv_path = "Data/receiver-deployments.csv", tz = attr(gl_det_data$first_det, "tzone")) %>%
    filter(tsStart < max(gl_det_data$last_det),
           (is.na(tsEnd) | tsEnd > min(gl_det_data$tagDeployStart)))
  
  # Gather geographic data
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  na <- rnaturalearth::ne_states(country = c("united states of america", "canada", "mexico"), returnclass = "sf")
  if (!is.null(map_crs)) {
    na <- st_transform(na, map_crs)
    # Convert limits to new crs
    lims <- data.frame(lon = xlim, lat = ylim) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_transform(map_crs) %>%
      st_coordinates() %>%
      as.data.frame()
    xlim <- lims$X
    ylim <- lims$Y
  }
  
  recvs_sf <- st_as_sf(recvs, coords = c("longitude", "latitude"), crs = 4326) %>%
    st_intersection(st_transform(gl_sf, crs = 4326))
  station_sf <- st_as_sf(station_n, coords = c("lon", "lat"), crs = 4326)
  gl_sf <- summarize(gl_sf)
  species <- unique(gl_det_data$spp)
  
  p <- ggplot(data = na) +
    geom_sf(fill = "gray95", color = "black") +
    geom_sf(data = world, fill = NA, color = "black", lwd = 1) +
    geom_sf(data = gl_sf , fill = NA, color = "grey50", lwd = 1.5) +
    geom_sf(data = gl_all_stops, shape = 21, color = NA, fill = stop_color) +
    geom_sf(data = recvs_sf, shape = 21, color = "grey5", 
            fill = "white", stroke = 0.25, size = 2) +
    
    geom_sf(data = station_sf, aes(size = pct_indivs), shape = 21, stroke = 0.5, color = "white", fill = det_color) +
    coord_sf(xlim = xlim, ylim = ylim) +
    # scale_size(range = c(2, 10)) +
    scale_size_binned(n.breaks = 5, range = c(2, 6)) +
    labs(x = NULL, y = NULL, size = paste0("% ", if (length(species) == 1) species else "individuals",
                                           "\ndetected (n = ", n_distinct(dets$motusTagID), ")")) + 
    ggspatial::annotation_scale(location = "bl", width_hint = 0.35, text_cex = 0.9) +
    # annotation_north_arrow(location = "bl", which_north = "true", 
    #                        pad_x = unit(0.1, "in"), pad_y = unit(0.5, "in"),
    #                        style = north_arrow_fancy_orienteering) +
    theme_bw(base_size = 14) + 
    theme(legend.background = element_rect(color = "black"),
          legend.direction = "horizontal",
          legend.position = legend_pos)
  return(p)
}
