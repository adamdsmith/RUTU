plot_gpsd_twnd_cor <- function(speed_dat, method = c("pearson", "spearman")) {
  method <- match.arg(method)

  # Trajectory groundspeed association with initial tailwind to destination receiver
  gspd_twnd_cor <- with(speed_dat, cor.test(groundspeed_net, tailwind, method = method))
  r_val <- round(gspd_twnd_cor$estimate, 2)
  p_val <- round(gspd_twnd_cor$p.value, 3)

  cor_label <- bquote(italic(r) == .(r_val))

  y_rng <- c(min(floor(speed_dat$tailwind)), max(ceiling(speed_dat$tailwind)))
  x_rng <- range(c(round_any(speed_dat$groundspeed_net, 2, floor),
                   round_any(speed_dat$groundspeed_net, 2, ceiling)))
  x_breaks <- seq(min(x_rng) + 2, max(x_rng) - 2, by = 4)
  
  p <-  
    ggplot(speed_dat, aes(x = groundspeed_net, y = tailwind)) +
    geom_point(shape = 21, fill = "grey50", color = "black", size = 5, alpha = 0.5) +
    scale_y_continuous("Initial tailwind (m/s)") +
    scale_x_continuous("Net ground speed (m/s)",
                       breaks = x_breaks,
                       limits = x_rng) +
    annotate(x = min(x_rng), y = max(y_rng), label = cor_label, 
             hjust = 0, vjust = 1, geom = "text", size = 6) +
    annotate(x = min(x_rng), y = max(y_rng) - diff(y_rng)/10, label = bquote(italic(p) == .(p_val)),
             hjust = 0, vjust = 1, geom = "text", size = 6) +
    theme_bw(base_size = 16)
  print(gspd_twnd_cor)
  return(p)
}
