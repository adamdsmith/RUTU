find_stations_near_deploy <- function(det_data, recvs = NULL, n_stns = NULL, coords = c("tagDepLon", "tagDepLat"),
                                 max_dist = 10, time_window = c(0, 0), tz = "America/New_York") {
  tags <- get_motusTagIDs(det_data)
  if (is.null(recvs)) recvs <- get_motus_recvs(tz = tz)
  if (is.null(n_stns)) n_stns <- nrow(recvs)
  deps <- get_deployment_info(det_data) %>% arrange(tagDeployID)
  near <- find_deploy_recvs(deps, recvs, max_dist, time_window = time_window)
  out <- left_join(deps, near, by = "tagDeployID") %>%
    left_join(recvs, by = "recvDeployID") %>%
    group_by(tagDeployID) %>%
    arrange(tagDeployID, distkm) %>%
    slice(seq(n_stns)) %>%
    ungroup()
  return(out)
}
