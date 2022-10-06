create_detection_maps <- function(det_data, outfile = "detection_maps.pdf",
                                  lat_rng = c(30, 55), lon_rng = c(-90, -60),
                                  add_stops = FALSE, stopover_dat = NULL, 
                                  verbose = TRUE) {
  old_options <- options(dplyr.summarise.inform = FALSE)
  on.exit(options(old_options))
  
  tmp_map_dat <- det_data %>%
    select(mfgID, motusTagID, label, dep_lat = tagDepLat, dep_lon = tagDepLon, tagDeployStart,
           date, recvDeployID:recvSiteName, antBin, age) %>%
    distinct()
  tmp_det_dat <- det_data %>%
    mutate(site_label = paste0(recvSiteName, "\n(", recv, ")")) %>%
    select(site_label, mfgID, motusTagID, ts, sig, runID, rlClass)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  na <- ne_states(country = c("united states of america", "canada"), returnclass = "sf")
  timez <- st_read("~/FWS_Projects/BLRA_SALS_SESP/SC_GA_SALS/SALS_analysis/Resources/geodata/ne_10m_time_zones.shp", quiet = TRUE) %>%
    st_transform(st_crs(world))
  
  ## Retrieve receivers in network
  recvs <- get_motus_recvs(tz = my_tz) %>%
    select(tsStart:isMobile)
  
  pdf(outfile, width = 10, height = 8.5)
  tags <- get_motusTagIDs(det_data)
  for (i in tags) {
    if (verbose) message("Processing Motus tag ", i)
    tmp_map <- filter(tmp_map_dat, motusTagID == i)
    tmp_det <- filter(tmp_det_dat, motusTagID == i) %>%
      mutate(site_sorted = forcats::fct_reorder(site_label, ts))
    tmp_site_lab <- tmp_map %>%
      group_by(recvDeployLat, recvDeployLon, recvSiteName) %>%
      summarize(label = paste0(max(recvSiteName), ":\n", dates_by_month(date))) %>%
      select(lat = recvDeployLat, lon = recvDeployLon, label) %>%
      mutate(type = "detection")
    tmp_tag_loc <- tmp_map %>% select(dep_lon, dep_lat, dep_date = tagDeployStart) %>% 
      distinct() %>%
      mutate(label = paste("Deployed:", format(as.Date(dep_date, tz = my_tz), format = "%d %b %Y"))) %>%
      select(lat = dep_lat, lon = dep_lon, label) %>%
      mutate(type = "deploy")
    all_labels <- bind_rows(tmp_site_lab, tmp_tag_loc) %>% 
      ungroup() %>%
      mutate(type = factor(type, levels = c("deploy", "detection"))) %>%
      arrange(lat)
    
    # Apparently active and functioning receivers from deployment to last detection
    tmp_recvs <- recvs %>%
      filter(tsStart < max(tmp_det$ts),
             (is.na(tsEnd) | tsEnd > min(tmp_map$tagDeployStart)))
    
    p <- ggplot(data = world) +
      geom_sf(fill = "black", color = NA) +
      geom_sf(data = na, fill = NA, color = "grey25") + 
      geom_sf(data = timez, fill = NA, color = "white", lty = "dashed") +
      geom_point(data = tmp_recvs, aes(longitude, latitude), pch = 21, color = "grey25", 
                 fill = "grey70", stroke = 0.25, size = 1.5) +
      geom_point(data = tmp_tag_loc, aes(lon, lat), 
                 shape = 21, fill = "#984ea3", color = "white", stroke = 0.5, size = 3)
    
    if (add_stops) {
      if (is.null(stopover_dat)) {
        add_stops <- FALSE
        warning("No stopover data provided. Skipping addition of stopovers to map...")
      } else {
        i_stop_solo <- filter(stopover_dat, motusTagID == i, stop_solo)
        i_stop_adj <- filter(stopover_dat, motusTagID == i, stop_adj)
        if (nrow(i_stop_adj) > 0) {
          p <- p +
            geom_point(data = i_stop_adj, aes(recvDeployLon, recvDeployLat), 
                       shape = 16, color = "#fdc086", alpha = 0.75, size = 7)
        }
        if (nrow(i_stop_solo) > 0) {
          p <- p +
            geom_point(data = i_stop_solo, aes(recvDeployLon, recvDeployLat), 
                       shape = 16, color = "#7fc97f", alpha = 0.75, size = 7)
        }
      }
    }

    p <- p + 
      geom_point(data = tmp_map, aes(x = recvDeployLon, y = recvDeployLat, shape = antBin), 
                 fill = "#ff7f00", color = "white", stroke = 0.5, size = 3) +
      scale_shape_manual(values=c(24, 22, 25, 23, 21), drop = FALSE)
    
    n_labels <- nrow(all_labels)
    if (n_labels < 15) {
      p <- p + 
        geom_label_repel(data = all_labels, aes(lon, lat, label = label, fill = type), 
                         size = 2.5, 
                         direction = "y", hjust = 0, 
                         nudge_x = pmin(max(pull(filter(all_labels, lat >= 25), lon)) + 6 - all_labels$lon, 12),
                         segment.color = "white", segment.size = 0.5, force = 5)
    } else {
      # Dividing labels to left and right side, with more southern on right
      mid_row <- ceiling(n_labels / 2)
      first <- filter(all_labels, row_number() <= mid_row)
      last <- filter(all_labels, row_number() > mid_row)
      
      p <- p +
        geom_label_repel(data = first,
                         aes(lon, lat, label = label, fill = type), 
                         size = 2.5, ylim = c(NA, 50),
                         direction = "y", hjust = 0, 
                         nudge_x = pmin(max(pull(filter(first, lat >= 25), lon)) + 6 - first$lon, 12),
                         segment.color = "white", segment.size = 0.5, force = 5) +
        geom_label_repel(data = last,
                         aes(lon, lat, label = label, fill = type), 
                         size = 2.5, ylim = c(35, NA),
                         direction = "y", hjust = 1, 
                         nudge_x = pmin(median(pull(filter(last, lat >= 25), lon)) - 6 - last$lon, -3),
                         segment.color = "white", segment.size = 0.5, force = 5)
      
    }
    p <- p +
      scale_fill_manual(values = c("#984ea3","#ff7f00"), 
                        drop = FALSE, guide = "none") +
      coord_sf(xlim = lon_rng + c(-2, 2),
               ylim = lat_rng + c(-2, 2),
               expand = FALSE) +
      ggtitle(gsub("NA", "Unk", unique(paste(tmp_map$label, tmp_map$age, sep = ", ")))) + 
      guides(shape = "none") +
      labs(x = NULL, y = NULL) + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "grey50", color = NA))
    print(p)
  }
  dev.off()
  system(paste('open', shQuote(outfile)))
}