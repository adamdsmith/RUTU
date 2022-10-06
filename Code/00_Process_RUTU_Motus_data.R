if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", quiet = TRUE)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", quiet = TRUE)
if (!requireNamespace("nrsmisc", quietly = TRUE)) remotes::install_github("adamdsmith/nrsmisc")
# Need development version of ggrepel
if (compareVersion("0.9.0", as.character(packageVersion("ggrepel"))) > 0)
  remotes::install_github("slowkow/ggrepel")
# If you get an error during this installation, restart R and run 
# `Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")` prior to the install attempt

pacman::p_load(dplyr, purrr, lubridate, motus, nrsmisc, ggplot2, ggforce, 
               gganimate, ggrepel, sf, rnaturalearth, patchwork)
Sys.setenv(TZ = "UTC")
my_tz <- "America/New_York"

source("Code/Functions/compile_detections.R")  
source("Code/Functions/get_resights.R")
source('Code/Functions/create_detection_profiles.R')
source("Code/Functions/utils.R")

# SE rutu projects
projects <- c(4, 140)
update_detections = FALSE
if (update_detections) {
  rutu <- compile_detections(projects, dir = "Data/Motus_raw",
                             minRunLen = 3, speciesID = 4630, 
                             tagProjID = as.list(projects))
  rutu4 <- filter(rutu, runLen >= 4)
  rutu3 <- filter(rutu, runLen == 3, longRun == 3)
  rutu <- bind_rows(rutu4, rutu3)
  saveRDS(rutu, file = "Data/Derived/rutu.rds")
} else rutu <- readRDS("Data/Derived/rutu.rds")

# Retrieve tag deployment metadata
rutu_deploy <- list.files("Data", "tag-deploy", full.names = TRUE) %>%
  map_df(readr::read_csv) %>%
  filter(grepl("Arenaria", speciesName)) %>%
  mutate(markerType = gsub("green", "GRN", markerType),
         markerType = gsub(" flag", "", markerType),
         markerType = gsub("dark", "DK", markerType),
         markerType = gsub("light", "LT", markerType),
         marker = paste(markerType, markerNumber, sep = ": "),
         marker = ifelse(marker == "NA: NA", "", marker),
         age = ifelse(is.na(age), "", age),
         bandNumber = as.numeric(gsub("-", "", bandNumber)),
         bdate = as_date(
           ymd(paste(utcYearStart, utcMonthStart, utcDayStart, sep = "-"), tz = my_tz))) %>%
  # Drop 4 early May birds since this is about late May RUTU
  filter(yday(bdate) > 140) %>%
  select(motusTagID = tagID, band_num = bandNumber, marker, age, wt_g = weightGrams, bdate)

if (!file.exists("Data/Derived/rutu_preprocessed.rds")) {
  # Data inspection cutoff (restrict to data collected through Nov 2019)
  rutu <- filter(rutu, ts <= as_datetime("2019-12-01 00:00:00", tz = my_tz))
  # DROP ALL TORRANCE, ORCHARD HILL, KOFFLER, RUTHVEN, JOHNSON'S MILLS,
  # NORTHERN MONTEZUMA WMA, RUSHTON FARMS, and JAY DRASHER FOR NOW
  # WAY TOO MUCH NOISE OR ALIASING AT THOSE SITES
  rutu <- filter(rutu, !grepl("Torrance|Orchard|Koffl|Ruth|Johnson's Mills|Northern Montezuma|Drasher|Rushton", recvSiteName))
  rutu <- left_join(rutu, rutu_deploy, by = "motusTagID")
  saveRDS(rutu, "Data/Derived/rutu_preprocessed.rds")
} else rutu <- readRDS("Data/Derived/rutu_preprocessed.rds")

# Create inspection sheet for logging accept/reject status of marginal detections
if (!file.exists("Data/rutu_inspection_log.csv")) create_inspection_sheet(rutu, out_dir = "Data")

reproduce_pdfs <- FALSE
if (reproduce_pdfs) {
  
  # Get band/flag resightings
  resight <- get_resights()

  create_detection_profiles(rutu, out_pdf = paste0("Output/rutu_detection_validation_", Sys.Date(), ".pdf"))
  
  # Detection history maps to help with vetting
  rutu_map_dat <- rutu %>%
    select(mfgID, motusTagID, label, marker, age, dep_lat = tagDepLat, dep_lon = tagDepLon, 
           tagDeployStart, date, recvDeployID:recvSiteName, antBin) %>%
    distinct()
  rutu_det_dat <- rutu %>%
    mutate(site_label = paste0(recvSiteName, "\n(", recv, ")")) %>%
    select(site_label, mfgID, motusTagID, ts, sig, runID, rlClass)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  na <- ne_states(country = c("united states of america", "canada"), returnclass = "sf")
  timez <- st_read("~/FWS_Projects/GIS/ne_10m_time_zones.shp", quiet = TRUE) %>%
    st_transform(st_crs(world))
  tags <- sort(unique(rutu$motusTagID))
  
  ## Retrieve receivers in network
  recvs <- readr::read_csv("https://motus.org/data/downloads/api-proxy/receivers/deployments?fmt=csv") %>%
    mutate(tsStart = as_datetime(tsStart, tz = "UTC"),
           tsEnd = as_datetime(tsEnd, tz = "UTC")) %>%
    filter(!isMobile,
           !is.na(latitude),
           !is.na(longitude))
  attr(recvs$tsStart, "tzone") <- attr(recvs$tsEnd, "tzone") <- my_tz
  
  pdf(paste0("Output/rutu_detection_maps_", Sys.Date(), ".pdf"), width = 10, height = 8.5)
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
      ggtitle(unique(paste(tmp_map$label, tmp_map$age, tmp_map$marker, sep = ", "))) + 
      guides(shape = "none") +
      labs(x = NULL, y = NULL) + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "grey50", color = NA))
    print(p)
  }
  dev.off()
  system(paste('open', shQuote(paste0("Output/rutu_detection_maps_", Sys.Date(), ".pdf"))))
}

# Example vetting for RUTU 28461 at Point Mouillee 
# peruse(rutu, 28461, "mouill")

if (!file.exists("Data/Derived/rutu_vetted.rds")) {
  # Filter records based on inspection
  failed_insp <- readr::read_csv("Data/rutu_inspection_log.csv") %>%
    filter(!keep) %>%
    pull(runID) %>% unique()
  rutu <- filter(rutu, !(runID %in% failed_insp))
  
  # Several deployments need time fixes
  # Apparently the following deployments (all Lotek receivers) were downloaded in timezone GMT-4
  # rather than GMT as Motus expects/requires
  # Thus, all detections from these deployments will have 4 hours added to the time
  recv_adj_time <- c(4686L, 4687L, 4688L, 4689L, 4690L, 4691L, 
                     4692L, 4693L, 4694L, 4696L, 4697L, 5390L)
  rutu <- mutate(rutu,
                 ts = ifelse(recvDeployID %in% recv_adj_time, ts + 4 * 60 * 60, ts),
                 ts = as_datetime(ts, tz = my_tz))
  saveRDS(rutu, "Data/Derived/rutu_vetted.rds")
} else rutu <- readRDS("Data/Derived/rutu_vetted.rds")

# Create vetted detection maps
if (reproduce_pdfs) source("Code/Functions/create_vetted_detection_maps.R")
