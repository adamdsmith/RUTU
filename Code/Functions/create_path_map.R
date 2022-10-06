create_path_map <- function(det_data, xlim = c(-100, -60), ylim = c(25, 60), 
                            legend_pos = c(0.28, 0.085), backdrop_sf = NULL,
                            backdrop_label_sf = NULL,
                            path_size = 2, path_color = "#1f78b4",
                            scale_loc = "bl",
                            map_crs = NULL) {
  
  if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc", quiet = TRUE)
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) install.packages("rnaturalearth", quiet = TRUE)
  if (!requireNamespace("ggspatial", quietly = TRUE)) install.packages("ggspatial", quiet = TRUE)
  
  deps <- select(det_data, motusTagID, lat = tagDepLat, lon = tagDepLon, dt = tagDeployStart) %>%
    distinct() %>%
    mutate(recvDeployID = NA_integer_, recvSiteName = NA_character_)
  dets <- select(det_data, motusTagID, recvDeployID, recvSiteName, lat = recvDeployLat, lon = recvDeployLon, dt = first_det)
  tag_loc <- tibble(lon = -80.89, lat = 32.0665) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  tag_labs <- tibble(lon = -79.8, lat = 31) %>%
    mutate(label = "Capture location") %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  tag_lines <- st_sfc(mapply(function(a,b) {st_cast(st_union(a,b),"LINESTRING")}, 
                             st_geometry(tag_labs), st_geometry(tag_loc), SIMPLIFY=FALSE)) %>%
    st_as_sf(crs = 4326) 
  
  # Restrict departures to those with detections
  deps <- filter(deps, motusTagID %in% unique(dets$motusTagID))
  
  path_dat <- bind_rows(dets, deps) %>%
    group_by(motusTagID) %>%
    arrange(motusTagID, dt) %>%
    ungroup()
  
  station_n <- path_dat %>%
    filter(!is.na(recvDeployID)) %>%
    group_by(recvSiteName) %>%
    summarize(n_indivs = n_distinct(motusTagID), .groups = "drop")
  
  station_locs <- select(path_dat, recvSiteName, lat, lon) %>%
    distinct() %>%
    group_by(recvSiteName) %>%
    # Keep only one location in case of slight variation in position among deployments
    slice(1) %>%
    ungroup()
  
  station_n <- left_join(station_n, station_locs, by = "recvSiteName") %>%
    mutate(pct_indivs = round(n_indivs / n_distinct(path_dat$motusTagID) * 100, 1))
  
  # Gather active receivers in network
  recvs <- get_motus_recvs(recv_path = "Data/receiver-deployments.csv", tz = attr(det_data$first_det, "tzone")) %>%
    filter(tsStart < max(det_data$last_det),
           (is.na(tsEnd) | tsEnd > min(det_data$tagDeployStart)))
  
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
  
  recvs_sf <- st_as_sf(recvs, coords = c("longitude", "latitude"), crs = 4326)
  station_sf <- st_as_sf(station_n, coords = c("lon", "lat"), crs = 4326)
  det_paths_sf <- st_as_sf(path_dat, coords = c("lon", "lat"), crs = 4326) %>%
    group_by(motusTagID) %>%
    summarize(do_union=FALSE, .groups = "keep") %>%
    st_cast("LINESTRING") %>%
    ungroup()

  species <- unique(det_data$spp)
  
  p <- ggplot(data = na) +
    geom_sf(fill = "gray95", color = "black") +
    geom_sf(data = world, fill = NA, color = "black", lwd = 1)
  
  if (!is.null(backdrop_sf)) {
    p <- p +
      geom_sf(data = backdrop_sf, fill = "grey50", alpha = 0.5, color = NA)
  }
  
  p <- p +
    geom_sf(data = recvs_sf, pch = 21, color = "grey5", 
               fill = "white", stroke = 0.5, size = 1.5) +
    geom_sf(data = det_paths_sf, color = path_color, alpha = 0.5, size = path_size) +
    geom_sf(data = station_sf, aes(size = pct_indivs), shape = 21, stroke = 0.5, color = "white", fill = path_color)+
    geom_sf(data = tag_lines, lwd = 0.5) +
    geom_sf_label(data = tag_labs, aes(label = label), hjust = 0, size = 3) 
  
  if (!is.null(backdrop_label_sf)) {
    p <- p +
      geom_sf_label(data = backdrop_label_sf, aes(label = label))
  }

  p <- p +   
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_size_binned(n.breaks = 5, range = c(1, 5)) +
    labs(x = NULL, y = NULL, size = paste0("% ", if (length(species) == 1) species else "individuals",
                                           "\ndetected (n = ", n_distinct(path_dat$motusTagID), ")")) + 
    ggspatial::annotation_scale(location = scale_loc, width_hint = 0.55, text_cex = 0.9) +
    theme_bw(base_size = 14) + 
    theme(legend.background = element_rect(color = "black"),
          legend.direction = "horizontal",
          legend.position = legend_pos)
  return(p)
}
