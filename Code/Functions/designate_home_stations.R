designate_home_stations <- function(det_data, nearby_stns, max_dist = 10, min_stns = 1, 
                                    enforce_season = FALSE, season_marg = 5,
                                    input = c("detection", "nearby")) {
  has_season <- "season" %in% names(det_data)
  if (enforce_season) {
    if (has_season) has_season <- all(det_data$season %in% c("spring", "fall"))
    if (!has_season) stop("Expecting a variable `season` in `det_data` indicating ",
                          "season of tag deployment: 'spring' or 'fall'")
    dep_season <- select(det_data, tagDeployID, season) %>%
      distinct()
    nearby_stns <- left_join(nearby_stns, dep_season, by = "tagDeployID")
    
    season_marg <- season_marg / 111 # convert kilometers into approx degrees of latitude
  }
  
  home_stns <- nearby_stns %>%
    arrange(tagDeployID, distkm) %>%
    group_by(tagDeployID)
  
  if (enforce_season)
    home_stns <- home_stns %>%
    mutate(diff = latitude - tagDepLat,
           home = FALSE,
           home = ifelse(season == "fall", latitude - tagDepLat <= season_marg & distkm <= max_dist, home),
           home = ifelse(season == "spring", latitude - tagDepLat >= -1 * season_marg & distkm <= max_dist, home),
           home = ifelse(row_number() <= min_stns, TRUE, home))
  else
    home_stns <- home_stns %>%
    mutate(diff = latitude - tagDepLat,
           home = distkm <= max_dist,
           home = ifelse(row_number() <= min_stns, TRUE, home))
  
  home_stns <- home_stns %>%
    ungroup() %>%
    filter(home) %>%
    select(tagDeployID, recvDeployID, home)
  
  input <- match.arg(input)
  join_dat <- if (identical(input, "detection")) det_data else nearby_stns
  det_data <- left_join(join_dat, home_stns, by = c("tagDeployID", "recvDeployID")) %>%
    mutate(home = ifelse(is.na(home), FALSE, home))
  return(det_data)
}
