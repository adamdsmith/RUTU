plot_first_gl_histogram <- function(det_summary_dat) {
  doy_rng <- c(round_any(min(det_summary_dat$first_detection), 2, floor),
               round_any(max(det_summary_dat$first_detection), 2, ceiling))
  doy_breaks <- seq(doy_rng[1], doy_rng[2], by = 3)
  doy_labels <- date_seq()[doy_breaks]
  p <- ggplot(det_summary_dat, aes(x = first_detection)) +
    geom_histogram(binwidth = 1, fill = "grey75", color = "black") +
    scale_x_continuous("Date of first detection in Great Lakes Basin", 
                       breaks = doy_breaks, labels = doy_labels, 
                       minor_breaks = seq(doy_rng[1], doy_rng[2], by = 1),
                       limits = doy_rng) +
    # Hack Y axis to make dotplot reflect number of individuals
    scale_y_continuous("# Ruddy Turnstone") +
    theme_bw(base_size = 16) +
    theme(legend.position = c(0.975, 0.975),
          legend.justification = c(1, 1),
          legend.background = element_rect(color = "black"))
  return(p)
}