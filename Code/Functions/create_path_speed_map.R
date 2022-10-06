create_path_speed_map <- function(speed_dat, xlim = c(-90, -70), ylim = c(32.5, 52.5), 
                            legend_pos = c(0.22, 0.1), path_size = 2, scale_loc = "br",
                            map_crs = NULL, base_size = 14) {
  
  if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc", quiet = TRUE)
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) install.packages("rnaturalearth", quiet = TRUE)
  if (!requireNamespace("ggspatial", quietly = TRUE)) install.packages("ggspatial", quiet = TRUE)
  
  speed_paths <- pbapply::pblapply(seq(nrow(speed_dat)), function(i) {
    tmp_spd <- speed_dat[i, ]
    id <- pull(tmp_spd, motusTagID)
    start <- pull(tmp_spd, start_recvDeployID)
    end <- pull(tmp_spd, end_recvDeployID)
    pt_df <- filter(rutu_det, motusTagID == id) %>%
      arrange(first_det)
    start_row <- which(pt_df$recvDeployID == start)
    end_row <- which(pt_df$recvDeployID == end)
    pt_df <- slice(pt_df, start_row:end_row)
    i_path <- pt_df %>%
      st_as_sf(coords = c("recvDeployLon", "recvDeployLat"), crs = 4326) %>%
      group_by(motusTagID) %>%
      summarize(do_union=FALSE, .groups = "drop") %>%
      st_cast("LINESTRING") %>% ungroup() %>%
      left_join(tmp_spd, by = "motusTagID")
  })
  speed_paths <- do.call("rbind", speed_paths)

  # Gather active receivers in network
  recvs <- get_motus_recvs(recv_path = "Data/receiver-deployments.csv", tz = attr(speed_dat$start_time, "tzone")) %>%
    filter(tsStart < max(speed_dat$end_time),
           (is.na(tsEnd) | tsEnd > min(speed_dat$start_time)))
  
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
  det_paths_sf <- speed_paths

  p <- ggplot(data = na) +
    geom_sf(fill = "gray95", color = "black") +
    geom_sf(data = world, fill = NA, color = "black", lwd = 1) +
    geom_sf(data = recvs_sf, pch = 21, color = "grey5", 
               fill = "white", stroke = 0.5, size = 1.5) +
    geom_sf(data = det_paths_sf, aes(color = groundspeed_net), alpha = 0.75, size = path_size) +   
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_color_gradient("Net ground\nspeed (m/s)", low = "#c7e9b4", high = "#081d58") +
    labs(x = NULL, y = NULL) + 
    ggspatial::annotation_scale(location = scale_loc, width_hint = 0.4, text_cex = 0.9) +
    theme_bw(base_size = base_size) + 
    theme(legend.background = element_rect(color = "black"),
          legend.direction = "horizontal",
          legend.position = legend_pos) +
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0))
  return(p)
}
