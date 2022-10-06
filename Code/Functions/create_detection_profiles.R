create_detection_profiles <- function(det_data, 
                                      out_pdf = "detection_profiles.pdf") {
  if (!requireNamespace("forcats", quietly = TRUE)) install.packages("forcats", quiet = TRUE)
  tags <- get_motusTagIDs(det_data)
  view_dets <- select(det_data, motusTagID, label, lab_day_site, tagDeployStart, ts, sig, key, 
                      runLen, rlClass, antBin, recvSiteName, antBearing, inspect)
  tf_colors <- c("grey50", "#b2182b"); names(tf_colors) <- c(FALSE,TRUE)
  pdf(file = out_pdf, height = 10.5, width = 8, onefile = TRUE)
  for (i in tags) {
    tmp_view_dets <- filter(view_dets, motusTagID == i) %>%
      mutate(lab_day_site = forcats::fct_reorder(lab_day_site, ts))
    n_pages <- ceiling(n_distinct(tmp_view_dets$lab_day_site) / 8)
    for (j in seq_len(n_pages)) {
      p <- ggplot(tmp_view_dets, aes(x = ts, y = sig)) + 
        geom_point(aes(fill = rlClass, shape = antBin, color = inspect, stroke = inspect), size = 2.5) +
        scale_fill_manual("Run length", values = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), drop = FALSE) + 
        scale_color_manual(values = tf_colors) + 
        scale_shape_manual(values=c(24, 22, 25, 23, 21), drop = FALSE) +
        scale_discrete_manual(aesthetics = "stroke", values = c(`TRUE` = 1.5, `FALSE` = 0), guide = "none") +
        scale_x_datetime(NULL, date_labels = "%H:%M:%S") +
        guides(shape = "none", color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
        facet_wrap_paginate(~ lab_day_site, scales = "free", nrow = 4, ncol = 2, page = j) +
        theme_black() +  
        theme(legend.position = "top", legend.direction = "horizontal",
              axis.title = element_blank()) +
        labs(y = "Signal strength") +
        ggtitle(unique(tmp_view_dets$label), 
                subtitle = paste("Deployed:", format(unique(tmp_view_dets$tagDeployStart), 
                                                     format = "%d %b %Y")))
      print(p)
    }
  }
  dev.off()
  system(paste('open', shQuote(out_pdf)))
}