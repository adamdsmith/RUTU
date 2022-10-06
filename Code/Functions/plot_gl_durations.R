plot_gl_durations <- function(gl_det_data) {
 
  gl_summary <- gl_det_data %>%
    group_by(motusTagID) %>%
    summarize(gl_stay_d = as.numeric(difftime(max(last_det), min(first_det), units = "days")), .groups = "drop")

  gl_stays <- gl_summary$gl_stay_d
  bins <- seq(round_any(min(gl_stays), 1, floor), round_any(max(gl_stays), 1, ceiling), by = 1)
  
  species <- unique(gl_det_data$spp)
  
  p <- ggplot(gl_summary, aes(x = gl_stay_d)) +
    geom_histogram(breaks = bins, fill = "grey75", color = "black") +
    scale_x_continuous("Minimum passage duration through\nthe Great Lakes Basin (days)", breaks = bins, labels = bins) + 
    labs(y = paste("#", if (length(species) == 1) species else "individuals")) +
    theme_bw(base_size = 16) +
    theme(panel.grid.minor.x = element_blank())
  return(p)
}