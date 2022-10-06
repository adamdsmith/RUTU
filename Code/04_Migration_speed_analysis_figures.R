source("Code/Functions/create_path_speed_map.R")
source("Code/Functions/plot_gspd_twnd_cor.R")

rutu_speeds <- readRDS("Data/Derived/rutu_speed_calculations.rds")
rutu_speeds_redux <- readr::read_csv("Data/Derived/rutu_speed_checks.csv") %>%
  filter(keep) 

# Number of flights
nrow(rutu_speeds_redux)

# Number of individuals
n_distinct(rutu_speeds_redux$motusTagID)

# Distance summary
median(rutu_speeds_redux$traj_len) / 1000
range(rutu_speeds_redux$traj_len) / 1000
median(rutu_speeds_redux$traj_net_disp) / 1000
range(rutu_speeds_redux$traj_net_disp) / 1000

# Create map of paths highlighting speed
p <- create_path_speed_map(rutu_speeds_redux, map_crs = map_crs)
ggsave("Output/rutu_path_speed_map.png", dpi = 600, height = 8.6, width = 5.9)
ggsave("Output/rutu_path_speed_map.pdf", height = 8.6, width = 5.9)
p_sm <- create_path_speed_map(rutu_speeds_redux, map_crs = map_crs, base_size = 14,
                              legend_pos = c(0.24, 0.89), scale_loc = "bl")
# Create additional component figures and assemble
gspd_fig <- 
  ggplot(rutu_speeds_redux, aes(x = groundspeed_net)) +
  geom_histogram(breaks = seq(2, 34, by = 2), fill = "grey75", color = "black") +
  annotate("text", x = 2, y = 7, label = paste("n =", nrow(rutu_speeds_redux)),
           hjust = 0, vjust = 1, size = 6) +
  scale_y_continuous("# Ruddy Turnstone", breaks = seq(0, 10, 2)) +
  scale_x_continuous(NULL, breaks = seq(4, 32, 4)) +
  theme_bw(base_size = 16)

# Groundspeed and tailwind correlation, plus plot
twnd_fig <- plot_gpsd_twnd_cor(rutu_speeds_redux)

# Assemble
all <- (gspd_fig / twnd_fig) | p_sm
all + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 22))
ggsave("Output/rutu_speed_fig.png", dpi = 600, height = 7, width = 10)

# ggsave("Output/rutu_path_speed_map.pdf", height = 8.6, width = 5.9)

sum(rutu_speeds_redux$tailwind > 0) / nrow(rutu_speeds_redux)
