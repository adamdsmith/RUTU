create_nearby_station_maps <- function(nearby_stns, deps = NULL, outfile = NULL, recv_labels = TRUE, enforce_season = FALSE) {
  if (is.null(deps)) deps <- get_deployment_info(nearby_stns)
  uniq_dep_locs <- select(deps, tagDepLat, tagDepLon) %>% 
    distinct()
  uniq_dep_locs$dep_loc_id <- seq(nrow(uniq_dep_locs))
  deps <- left_join(deps, uniq_dep_locs, by = c("tagDepLat", "tagDepLon"))
  has_season <- "season" %in% names(deps)
  
  nearby_stns <- designate_home_stations(deps, nearby_stns, enforce_season = enforce_season, input = "nearby")
  
  out <- vector(mode = "list", length = nrow(deps))
  out_nms <- vector(mode = "integer", length = nrow(deps))
  i <- 0
  
  for (l in sort(unique(deps$dep_loc_id))) {
    tmp_deps <- filter(deps, dep_loc_id == l)
    tmp_dats <- filter(nearby_stns, tagDeployID %in% tmp_deps$tagDeployID)
    dep_loc <- c(lon = min(tmp_deps$tagDepLon), lat = min(tmp_deps$tagDepLat))
    max_dist <- max(tmp_dats$distkm)
    zm_adj <- as.integer(cut(max_dist, breaks = c(-1, 10, 20, 75, 165, 400, 1000)))
    zm <- 12 - zm_adj
    ggm <- ggmap::get_map(dep_loc, zoom = zm)
    tags <- sort(get_motusTagIDs(tmp_deps))
    
    for (t in tags) {
      i <- i + 1
      out_nms[i] <- t
      tmp_dat <- filter(tmp_dats, motusTagID == t)
      tmp_dep <- filter(tmp_deps, motusTagID == t)
      
      p <- ggmap(ggm)
      has_home <- "home" %in% names(tmp_dat)
      
      if (has_home) {
        tmp_dat$home <- factor(tmp_dat$home, levels = c(FALSE, TRUE))
        p <- p + geom_point(data = tmp_dat, aes(longitude, latitude, color = home), 
                            shape = 7, size = if (is.null(outfile)) 5 else 3) +
          scale_color_manual("Home station", values = c("#b2182b", "#2166ac"), drop=FALSE) +
          theme(legend.justification = c(0, 1), 
                legend.position = c(0.01, 0.99),
                legend.key = element_blank())
      } else {
        p <- p + geom_point(data = tmp_dat, aes(longitude, latitude), shape = 7, size = if (is.null(outfile)) 5 else 3)
      }
      
      if (recv_labels)
        p <- p + geom_text(data = tmp_dat, aes(longitude, latitude, label = siteName), size = if (is.null(outfile)) 4 else 2,
                           hjust = 0, vjust = 0.5, position = position_nudge(x = 0.0075 * zm_adj))
      
      p <- p +
        geom_point(data = tmp_dep, aes(tagDepLon, tagDepLat), shape = 23, fill = "deepskyblue1", size = if (is.null(outfile)) 2.5 else 1) +
        labs(x = NULL, y = NULL) +
        ggtitle(paste0("Motus tag: ", t, "; deployed ", format(tmp_dep$tagDeployStart, format = "%d %b %Y")))
      out[[i]] <- p
    }
  }
  names(out) <- out_nms
  out <- out[as.character(sort(out_nms))]
  if (!is.null(outfile)) {
    pdf(outfile, width = 5, height = 5.25)
    for (i in seq_along(out)) {
      print(out[[i]])
    }
    dev.off()
    system(paste('open', shQuote(outfile)))
  } else {
    for (i in seq_along(out)) {
      print(out[[i]])
    }
  }
}