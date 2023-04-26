source("Code/Functions/create_path_map.R")
source("Code/Functions/utils.R")
source("Code/Functions/plot_gl_durations.R")
source("Code/Functions/plot_capture_weights.R")
source("Code/Functions/plot_first_gl_histogram.R")
source("Code/Functions/plot_gl_wind_cor.R")
source("Code/Functions/create_gl_detection_map.R")

rutu_det <- readRDS("Data/Derived/rutu_detections.rds") %>%
  mutate(doy = yday(first_det)) %>%
  filter(doy <= 200)

detected_tags <- unique(rutu_det$motusTagID)
p <- plot_capture_weights(rutu_deploy, detected_tags, has_stop)
ggsave("Output/rutu_capture_weights.png", dpi = 600, height = 5, width = 5)

map_crs = "+proj=lcc +lat_1=30 +lat_2=55 +lat_0=42.5 +lon_0=-85 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs"
gl <- st_read("Resources/geodata/greatlakes_subbasins.shp", quiet = TRUE) %>%
  summarize()
db <- st_read("Resources/geodata/delaware-bay_HUC02040204.shp", quiet = TRUE) %>%
  st_buffer(30000) %>% summarize() %>% st_transform(st_crs(gl))
backdrop <- rbind(gl, db)
backdrop_labs <- tibble(lat = c(38, 49.5, 58.5, 53.75), 
                        lon = c(-71.25, -91.1, -85, -77), 
                        label = c("Delaware Bay", "Great Lakes\nBasin", "Hudson Bay", "James Bay")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
p <- create_path_map(rutu_det, map_crs = map_crs, backdrop_sf = backdrop, backdrop_label_sf = backdrop_labs)
ggsave("Output/Fig1_rutu_path_map.png", dpi = 600, height = 8.6, width = 6.8)
ggsave("Output/Fig1_rutu_path_map.pdf", height = 8.6, width = 6.8)

rutu_det_sf <- mutate(rutu_det,
                      lat = recvDeployLat,
                      lon = recvDeployLon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

rutu_det_gl <- st_intersection(rutu_det_sf, st_transform(gl, crs = 4326)) %>%
  as.data.frame()

# Histogram of minimum passage duration through Great Lakes Basin
gl_dur_p <- plot_gl_durations(rutu_det_gl)

# Plot first detection in GL Basin
rutu_gl_summary <- rutu_det_gl %>%
  group_by(motusTagID, tagDeployStart) %>%
  summarize(first_detection = min(doy),
            .groups = "drop")
gl_sum_p <- plot_first_gl_histogram(rutu_gl_summary)

all <- gl_dur_p | gl_sum_p
all + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 22))
ggsave("Output/Fig5_rutu_gl_passage.png", dpi = 600, height = 5, width = 10)

# Elapsed time between last GL detection and first James/Hudson Bay detection
last_gl <- rutu_det_gl %>%
  group_by(motusTagID) %>% 
  arrange(desc(last_det)) %>%
  slice(1)
first_breed <- rutu_det %>%
  filter(recvDeployLat > 50) %>%
  group_by(motusTagID) %>%
  arrange(first_det) %>%
  slice(1) %>%
  filter(motusTagID %in% unique(last_gl$motusTagID))
rutu_gl_breed <- filter(last_gl, motusTagID %in% unique(first_breed$motusTagID)) %>%
  bind_rows(first_breed) %>%
  group_by(motusTagID) %>%
  arrange(motusTagID, first_det) %>%
  summarize(tdiff = difftime(max(first_det), min(last_det), units = "hours"))
median(rutu_gl_breed$tdiff)
range(rutu_gl_breed$tdiff)

# Correlation between longitude of first detection in GL Basin and E-W wind component
p <- plot_gl_wind_cor(rutu_det_gl)
ggsave("Output/Fig4_detection_location_uwind_correlation.png", dpi = 600, height = 5, width = 5)

p <- create_gl_detection_map(rutu_det_gl, gl_sf = gl, map_crs = map_crs, stop_buffer = 40)
ggsave("Output/Fig3_rutu_gl_detection_map.png", dpi = 600, height = 5.3, width = 6.8)
ggsave("Output/Fig3_rutu_gl_detection_map.pdf", height = 5.3, width = 6.8)

