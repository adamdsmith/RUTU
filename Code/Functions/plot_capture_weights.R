plot_capture_weights <- function(deploy_data, detectedIDs = NULL, stopoverIDs = NULL, binwidth = 10) {
  wts <- deploy_data$wt_g
  bins <- seq(round_any(min(wts), binwidth, floor), round_any(max(wts), binwidth, ceiling), by = binwidth)
  breaks <- nrsmisc::every_nth(bins, 2, empty = FALSE, inverse = TRUE)
  
  # Base plot of all tag deployments
  p <- ggplot(deploy_data, aes(x = wt_g)) +
    geom_histogram(breaks = bins, fill = "grey75", color = "black")
  
  # Superimpose distribution of those detected during migration
  if (!is.null(detectedIDs)) {
    deploy_data_sub <- filter(deploy_data, motusTagID %in% detectedIDs)
    p <- p + geom_histogram(data = deploy_data_sub, fill = "grey25", breaks = bins,
                         color = "black")
  }

  # Superimpose distribution of those with identified stopovers during migration
  if (!is.null(stopoverIDs)) {
    deploy_data_sub <- filter(deploy_data, motusTagID %in% stopoverIDs)
    p <- p + geom_histogram(data = deploy_data_sub, fill = "white", breaks = bins,
                            color = "black")
  }

  p <- p +
    scale_x_continuous("Capture weight (g)", breaks = breaks) +
    labs(y = "# Ruddy Turnstones") +
    theme_bw(base_size = 16)
  return(p)
}