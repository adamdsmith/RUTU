identify_stopovers <- function(det_data, adj_dist = 30, solo_h = 8, adj_h = 10) {
  message("Identifying stopovers...")
  tags <- get_motusTagIDs(det_data)
  # Calculate distance matrix (km) between all receiver deployments in detection data set
  used_recvs <- get_uniq_recvs(det_data)
  recv_dist_matrix <- used_recvs %>% 
    sf::st_as_sf(coords = c("recvDeployLon", "recvDeployLat"), crs = 4326) %>%
    sf::st_distance() %>% units::drop_units()
  recv_dist_matrix <- recv_dist_matrix / 1000
  rownames(recv_dist_matrix) <- colnames(recv_dist_matrix) <- used_recvs$recvDeployID
  stops <- vector(mode = "list", length = length(tags))
  names(stops) <- tags
  current_stop_id <- 0
  pb = txtProgressBar(min = 0, max = length(stops), initial = 0, style = 3, width = 50) 
  for (t in tags) {
    # message(t)
    tmp_home <- filter(det_data, motusTagID == t, home)
    tmp_away <- filter(det_data, motusTagID == t, !home)
    away_recvs <- as.character(unique(tmp_away$recvDeployID))
    
    if (n_distinct(away_recvs) > 0) {
      tmp_recv_dist_matrix <- recv_dist_matrix[away_recvs, away_recvs, drop = FALSE]
      away_recvs <- as.integer(away_recvs)
      # Cycle through receivers
      r_tmp <- lapply(away_recvs, function (r) {
        # message(r)
        # Is there a stopover based on individual receiver?
        solo_det <- filter(tmp_away, recvDeployID == r)
        if (nrow(tmp_home) == 0)
          return_home <- FALSE
        else
          return_home <- max(tmp_home$last_det) > solo_det$last_det
        stop_solo <- solo_det$det_wind_h >= solo_h
        if (stop_solo) {
          # Only consider completely uninterupted detections 
          other_det <- filter(tmp_away, 
                              recvDeployID != r,
                              (between(first_det, solo_det$first_det, solo_det$last_det) |
                                 between(last_det, solo_det$first_det, solo_det$last_det)))
          if (nrow(other_det) > 0) stop_solo <- FALSE
        }
        solo_out <- select(solo_det, motusTagID, recvDeployID) %>%
          mutate(stop_solo = stop_solo,
                 stop_solo_h = ifelse(stop_solo, solo_det$det_wind_h, NA_real_))
        
        # Any detections at adjacent receivers, as defined by `adj_dist`?
        other_away <- as.character(setdiff(away_recvs, r))
        any_other_away <- n_distinct(other_away) > 0
        if (any_other_away) {
          r_recv_dist_matrix <- tmp_recv_dist_matrix[as.character(r), other_away]
          if (length(r_recv_dist_matrix) == 1) names(r_recv_dist_matrix) <- other_away
          is_adj <- r_recv_dist_matrix <= adj_dist
          any_adj <- any(is_adj)
        } else any_adj <- FALSE

        # Is there a stopover shared amongst any adjacent receivers?
        if (any_adj && !stop_solo) {
          adj_recvs <- as.integer(names(r_recv_dist_matrix)[is_adj])
          adj_det <- filter(tmp_away, recvDeployID %in% c(r, adj_recvs))
          adj_out <- adj_det %>%
            group_by(motusTagID) %>%
            summarize(adj_recv_grp = paste(sort(recvDeployID), collapse = ","),
                      first_stop_det = min(first_det),
                      last_stop_det = max(last_det),
                      adj_det_wind_h = as.numeric(difftime(last_stop_det, first_stop_det, units = "hours")),
                      .groups = "drop") %>%
            mutate(stop_adj = adj_det_wind_h >= adj_h,
                   stop_adj_h = ifelse(stop_adj, adj_det_wind_h, NA_real_))
        } else {
          adj_out <- tibble(motusTagID = t,
                            adj_recv_grp = NA_character_,
                            first_stop_det = if (stop_solo) solo_det$first_det else NA_Date_,
                            last_stop_det = if (stop_solo) solo_det$last_det else NA_Date_,
                            adj_det_wind_h = NA_real_,
                            stop_adj = FALSE,
                            stop_adj_h = NA_real_)
        }
        r_det <- left_join(solo_out, adj_out, by = "motusTagID") %>%
          mutate(home_after = return_home)
        r_det
      })
      r_tmp <- bind_rows(r_tmp) %>%
        mutate(is_stop = stop_solo | stop_adj)
      any_stops <- any(r_tmp$is_stop)
      if (any_stops) {
        stop_tz <- attr(r_tmp$first_stop_det, "tz")
        stop_ids <- filter(r_tmp, is_stop) %>%
          arrange(first_stop_det, last_stop_det) %>%
          mutate(stop_int = lubridate::interval(first_stop_det, last_stop_det, stop_tz),
                 stop_id = id_intervals(stop_int, current_stop_id)) %>%
          select(recvDeployID, stop_id)
        current_stop_id <- max(stop_ids$stop_id, na.rm = TRUE)
        r_tmp <- left_join(r_tmp, stop_ids, by = "recvDeployID")
      } else r_tmp <- mutate(r_tmp, stop_id = NA_integer_)
    } else {
      # Detections only at "home" stations, so it does not contribute here
      r_tmp <- NULL
    }
    stops[[t]] <- r_tmp
    setTxtProgressBar(pb, which(tags == t))
  }
  close(pb)
  stops <- bind_rows(stops)
  out <- left_join(det_data, stops, by = c("motusTagID", "recvDeployID")) %>%
    # Convert missing values for "home" receivers to FALSE
    mutate(stop_solo = ifelse(home, FALSE, stop_solo),
           stop_adj = ifelse(home, FALSE, stop_adj),
           home_after = ifelse(home, FALSE, home_after),
           is_stop = ifelse(home, FALSE, is_stop))
  return(out)
}