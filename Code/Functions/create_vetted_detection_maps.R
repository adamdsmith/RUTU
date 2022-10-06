# Detection history maps to help with vetting
rutu_map_dat <- rutu %>%
  select(mfgID, motusTagID, label, marker, age, dep_lat = tagDepLat, dep_lon = tagDepLon, tagDeployStart,
         date, recvDeployID:recvSiteName, antBin) %>%
  distinct()

add_weight <- TRUE
if (add_weight) {
  stopifnot(exists("rutu_deploy"))
  rutu_map_dat <- left_join(rutu_map_dat, select(rutu_deploy, motusTagID, wt_g), by = "motusTagID")
}

rutu_det_dat <- rutu %>%
  mutate(site_label = paste0(recvSiteName, "\n(", recv, ")")) %>%
  select(site_label, mfgID, motusTagID, ts, sig, runID, rlClass)
world <- ne_countries(scale = "medium", returnclass = "sf")
na <- ne_states(country = c("united states of america", "canada"), returnclass = "sf")
timez <- st_read("~/FWS_Projects/GIS/ne_10m_time_zones.shp", quiet = TRUE) %>%
  st_transform(st_crs(world))
tags <- sort(unique(rutu$motusTagID))

# Get band/flag resightings
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
  inner_join(rutu_deploy, by = "band_num") %>%
  mutate(dtime = time_length(rdate - bdate, "year")) %>%
  filter(dtime < 0.5)

## Retrieve receivers in network
recvs <- readr::read_csv("https://motus.org/data/downloads/api-proxy/receivers/deployments?fmt=csv") %>%
  mutate(tsStart = as_datetime(tsStart, tz = "UTC"),
         tsEnd = as_datetime(tsEnd, tz = "UTC")) %>%
  filter(!isMobile,
         !is.na(latitude),
         !is.na(longitude))
attr(recvs$tsStart, "tzone") <- attr(recvs$tsEnd, "tzone") <- my_tz

pdf(paste0("Output/rutu_detection_maps_vetted_", Sys.Date(), ".pdf"), width = 10, height = 8.5)
for (i in tags) {
  tmp_map <- filter(rutu_map_dat, motusTagID == i)
  tmp_det <- filter(rutu_det_dat, motusTagID == i) %>%
    mutate(site_sorted = forcats::fct_reorder(site_label, ts))
  tmp_res <- filter(resight, motusTagID == i)
  has_res <- nrow(tmp_res) > 0
  if (has_res) {
    tmp_res_lab <- tmp_res %>%
      group_by(lat, lon) %>%
      summarize(label = paste("Resight:", dates_by_month(rdate))) %>%
      mutate(type = "resight")
  }
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
  all_labels <- bind_rows(tmp_site_lab, tmp_tag_loc) %>% ungroup()
  if (has_res) all_labels <- bind_rows(all_labels, tmp_res_lab)
  all_labels <- mutate(all_labels, 
                       type = factor(type, levels = c("deploy", "detection", "resight"))) %>%
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
               shape = 21, fill = "#984ea3", color = "white", stroke = 0.5, size = 3) + 
    geom_point(data = tmp_map, aes(x = recvDeployLon, y = recvDeployLat, shape = antBin), 
               fill = "#ff7f00", color = "white", stroke = 0.5, size = 3) +
    scale_shape_manual(values=c(24, 22, 25, 23, 21), drop = FALSE)
  
  if (has_res) {
    p <- p +
      geom_point(data = tmp_res, aes(x = lon, y = lat), pch = 21, 
                 fill = "#4daf4a", color = "white", stroke = 0.5, size = 3)
  }
  
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
    scale_fill_manual(values = c("#984ea3","#ff7f00","#4daf4a"), 
                      drop = FALSE, guide = "none") +
    coord_sf(xlim = c(-105, -50) + c(-1, 1),
             ylim = c(27, 58) + c(-2, 2),
             expand = FALSE) +
    ggtitle(unique(paste(tmp_map$label, tmp_map$age, if(add_weight) paste("wt:", tmp_map$wt_g, "g"), tmp_map$marker, sep = ", "))) + 
    guides(shape = "none") +
    labs(x = NULL, y = NULL) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey50", color = NA))
  print(p)
}
dev.off()
system(paste('open', shQuote(paste0("Output/rutu_detection_maps_vetted_", Sys.Date(), ".pdf"))))
