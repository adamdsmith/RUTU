calculate_ground_speed <- function(det_data, min_stn_sep = 150, max_time_sep = 12) {
  
  if (!requireNamespace("trajr", quietly = TRUE)) install.packages("trajr", quiet = TRUE)
  
  message("Estimating ground speed of flights between detections...")
  tags <- get_motusTagIDs(det_data)
  
  # Calculate distance matrix (km) between all receiver deployments in detection data set
  used_recvs <- get_uniq_recvs(det_data)
  recv_dist_matrix <- used_recvs %>% 
    sf::st_as_sf(coords = c("recvDeployLon", "recvDeployLat"), crs = 4326) %>%
    sf::st_distance() %>% units::drop_units()
  rownames(recv_dist_matrix) <- colnames(recv_dist_matrix) <- used_recvs$recvDeployID
  
  ground_speeds <- pbapply::pblapply(tags, function(t) {
    # message(t)
    tmp <- filter(det_data, motusTagID == t) %>% arrange(first_det)
    t_recvs <- tmp %>% pull(recvDeployID) %>% unique() %>% as.character()
    n_recvs <- length(t_recvs)
    
    if (n_recvs < 2) {
      r_tmp <- NULL
    } else {
      tmp_recv_dist_matrix <- recv_dist_matrix[t_recvs, t_recvs]

      # Cycle through receivers
      r_tmp <- lapply(head(t_recvs, -1), function(r) {
        # message(r)
        # Only need receivers encountered after the current one
        r_col <- which(colnames(tmp_recv_dist_matrix) == r)
        r_recv_dist_matrix <- tmp_recv_dist_matrix[r_col, (r_col + 1):n_recvs, drop = FALSE]
        away_recvs <- colnames(r_recv_dist_matrix)

        # Last detection at beginning receiver
        r_det <- filter(tmp, recvDeployID == as.integer(r)) %>%
          select(motusTagID:tagDeployID, recvSiteName:recvDeployLon, det_time = last_det)
        # First detection at arriving receiver
        a_det <- filter(tmp, recvDeployID %in% as.integer(away_recvs)) %>%
          select(motusTagID:tagDeployID, recvSiteName:recvDeployLon, det_time = first_det) %>%
          # Exclude detections beyond acceptable time separation
          mutate(keep = difftime(det_time, r_det$det_time, units = "hours") <= max_time_sep) %>%
          filter(keep) %>% select(-keep)
        all_det <- bind_rows(r_det, a_det)
        over_dist <- tmp_recv_dist_matrix[as.character(head(all_det$recvDeployID, 1)),
                                          as.character(tail(all_det$recvDeployID, 1))] >= (min_stn_sep * 1000)
        if (nrow(all_det) < 2 | !over_dist) {
          r_tmp <- NULL
        } else {
          
          # Create trajectory
          traj_crs <- sprintf("+proj=aeqd +lat_0=%1.0f +lon_0=%1.0f", 
                              mean(all_det$recvDeployLat), mean(all_det$recvDeployLon))
          # Calculate trip metrics
          r_traj <- all_det %>%
            st_as_sf(coords = c("recvDeployLon", "recvDeployLat"), crs = 4326) %>%
            st_transform(traj_crs) %>% st_coordinates() %>% as.data.frame() %>%
            mutate(det_time = as.numeric(all_det$det_time)) %>%
            trajr::TrajFromCoords(timeCol = 3)
          traj_net_disp <- trajr::TrajDistance(r_traj)
          traj_len <- trajr::TrajLength(r_traj)
          traj_dur_s <- trajr::TrajDuration(r_traj)
          # traj_gspd_net <- traj_net_disp / (traj_dur_h * 60 * 60)
          # traj_gspd <- traj_len / (traj_dur_h * 60 * 60)
          traj_straight <- trajr::TrajStraightness(r_traj)
          
          # Build output
          r_det <- select(r_det, motusTagID:tagDeployID, recvSiteName:recvDeployLon, start_time = det_time) %>%
            rename_with(~ paste0("start_", .x), starts_with("recv"))
          traj_meta <- a_det %>% slice(nrow(.)) %>%
            select(motusTagID, recvSiteName:recvDeployLon, end_time = det_time) %>%
            rename_with(~ paste0("end_", .x), starts_with("recv")) %>%
            mutate(traj_len = traj_len,
                   traj_net_disp = traj_net_disp,
                   traj_dur_h = traj_dur_s / 60 / 60,
                   traj_straightness = traj_straight,
                   groundspeed_net = traj_net_disp / traj_dur_s,
                   groundspeed_traj = traj_len / traj_dur_s,
                   traj = list(r_traj))
          r_det <- left_join(r_det, traj_meta, by = "motusTagID") %>%
            mutate(brng_init = (geosphere::bearing(c(start_recvDeployLon, start_recvDeployLat),
                                                   c(end_recvDeployLon, end_recvDeployLat)) + 360) %% 360)
        }
      })
      if (all(sapply(r_tmp, is.null))) 
        r_tmp <- NULL
      else
        r_tmp <- bind_rows(r_tmp)
    }
    r_tmp
  })
  if (all(sapply(ground_speeds, is.null))) {
    message("No detections met your criteria to calculate ground speeds.")
    ground_speeds <- NULL
  } else ground_speeds <- bind_rows(ground_speeds)
  return(ground_speeds)
}
