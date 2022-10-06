source("Code/00_Process_RUTU_Motus_data.R")
pacman::p_load(ggmap)
# Google API key required for creating nearby station maps
register_google(Sys.getenv("goog_key"))

source("Code/Functions/find_stations_near_deploy.R")
source("Code/Functions/create_nearby_station_maps.R")
source("Code/Functions/designate_home_stations.R")
source("Code/Functions/identify_stopovers.R")
source('Code/Functions/retrieve_detection_wind_vectors.R')
source("Code/Functions/calculate_ground_speed.R")
source("Code/Functions/retrieve_path_wind_components.R")
source("Code/Functions/create_detection_maps.R")

# Summarize detections
if (!file.exists("Data/Derived/rutu_detections.rds")) {
  # Purely academic, there were no nearby "home" stations for northbound movement
  nearby_stns <- find_stations_near_deploy(rutu, max_dist = 20, time_window = c(0, 60))
  
  rutu <- rutu %>%
    designate_home_stations(nearby_stns) %>%
    group_by(tagDeployID) %>%
    mutate(only_home = all(home) & n_distinct(recvSiteName) == 1,
           only_away = all(!home)) %>%
    ungroup()

  rutu_det <- rutu %>%
    group_by(motusTagID, spp, age, 
             tagDeployID, tagDeployStart, tagDepLat, tagDepLon, 
             recvSiteName, recvDeployID, recvDeployLat, recvDeployLon, home, only_away) %>%
    summarize(first_det = min(ts),
              last_det = max(ts), .groups = "drop") %>%
    arrange(motusTagID, first_det) %>%
    mutate(det_wind_h = as.numeric(difftime(last_det, first_det, units = "hours"))) %>%
    identify_stopovers(solo_h = 4, adj_h = 6) %>%
    retrieve_detection_wind_vectors()
  saveRDS(rutu_det, file = "Data/Derived/rutu_detections.rds")
} else rutu_det <- readRDS("Data/Derived/rutu_detections.rds")

# Update detection maps to include stopovers
has_stop <- filter(rutu_det, stop_solo | stop_adj) %>% pull(motusTagID) %>% unique()
if (!file.exists("Output/rutu_detection_maps_w_stopovers.pdf"))
  create_detection_maps(det_data = filter(rutu, motusTagID %in% has_stop),
                        outfile = "Output/rutu_detection_maps_w_stopovers.pdf",
                        lat_rng = c(27, 58), lon_rng = c(-105, -50),
                        add_stops = TRUE, stopover_dat = rutu_det)

# Estimate orthodromic flight (ground) speed where appropriate
if (!file.exists("Data/Derived/rutu_speed_calculations.rds")) {
  rutu_speeds <- calculate_ground_speed(rutu_det, min_stn_sep = 150, max_time_sep = 18) %>%
    retrieve_path_wind_components()
  saveRDS(rutu_speeds, file = "Data/Derived/rutu_speed_calculations.rds")
} 
if (!file.exists("Data/Derived/rutu_speed_checks.csv")) {
  rutu_speeds %>% 
    select(-traj) %>%
    mutate(keep = TRUE) %>%
    readr::write_csv("Data/Derived/rutu_speed_checks.csv")
}

# Prepare csv for manual assignment of departure windows
rutu_dep <- rutu_det %>%
  group_by(motusTagID, spp, age, tagDeployID, tagDeployStart, tagDepLat, tagDepLon) %>%
  summarize(min_depart = "YYYY-MM-DD",
            max_depart = "YYYY-MM-DD",
            departure_notes = NA_character_)
if (!file.exists("Data/Derived/rutu_departures.csv")) {
  readr::write_csv(rutu_dep, "Data/Derived/rutu_departures.csv")
}