plot_gl_wind_cor <- function(gl_det_data, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  gl_first_det <- gl_det_data %>%
    group_by(motusTagID) %>%
    arrange(first_det) %>%
    slice(1)
  
  # East-west wind association with longitude of first detection...
  wind_loc_cor <- with(gl_first_det, cor.test(uwind, recvDeployLon, method = method))
  r_val <- round(wind_loc_cor$estimate, 2)
  p_val <- round(wind_loc_cor$p.value, 3)
  cor_label <- bquote(italic(r) == .(r_val))
  
  # Make Y-axis symmetrical for labelling purposes
  y_rng <- max(ceiling(abs(gl_first_det$uwind))) * c(-1, 1)
  x_breaks <- seq(min(floor(gl_first_det$recvDeployLon)), max(ceiling(gl_first_det$recvDeployLon)), by = 2)
  make_w_long <- function(x) {
    if (all(x < 0)) x <- x * -1
    paste0(x, "\u00b0W")
   }
  p <- ggplot(gl_first_det, aes(recvDeployLon, uwind)) +
    geom_jitter(shape = 21, fill = "grey50", color = "black", size = 5, alpha = 0.5,
                width = diff(range(x_breaks)) / 100, height = diff(y_rng) / 100) +
    scale_y_continuous("East        Wind component (m/s)        West",
                       limits = y_rng) +
    scale_x_continuous("Longitude of initial detection", 
                       breaks = x_breaks,
                       limits = range(x_breaks),
                       labels = make_w_long(x_breaks)) +
    annotate(x = min(x_breaks), y = max(y_rng) - diff(y_rng)/25, label = cor_label, 
             hjust = 0, vjust = -0.25, geom = "text", size = 6) +
    annotate(x = min(x_breaks), y = max(y_rng) - diff(y_rng)/25, label = bquote(italic(p) == .(p_val)),
             hjust = 0, vjust = 1.25, geom = "text", size = 6) +
    annotate(x = min(x_breaks), y = max(y_rng) - diff(y_rng)/25, label = paste("n =", nrow(gl_first_det)), 
             hjust = 0, vjust = 2.75, geom = "text", size = 6) +
    theme_bw(base_size = 16)
  print(wind_loc_cor)
  return(p)
}
