compile_detections <- function(projRecv, update = TRUE, dir = getwd(), 
                               tz = "America/New_York",
                               speciesID = NULL, tagProjID = NULL,
                               minRunLen = 3, medFreqsd = 0.05) {
  
  # Set up data directories and speciesID and tagProjectID filters
  stopifnot(length(dir) == 1 | identical(length(projRecv), length(dir)))
  if (!is.list(dir)) {
    tmpdir <- vector("list", length(projRecv))
    tmpdir[] <- dir
    dir <- tmpdir
  }
  stopifnot(is.null(speciesID) | length(speciesID) == 1 | identical(length(projRecv), length(speciesID)))
  if (!is.list(speciesID)) {
    spID <- vector("list", length(projRecv))
    spID[] <- speciesID
  } else spID <- speciesID
  stopifnot(is.null(tagProjID) | length(tagProjID) == 1 | identical(length(projRecv), length(tagProjID)))
  if (!is.list(tagProjID)) {
    tpID <- vector("list", length(projRecv))
    tpID[] <- tagProjID
  } else tpID <- tagProjID
  
  app_fun <- if(length(projRecv) > 1) pbapply::pblapply else lapply
  out <- app_fun(projRecv, function(p) {
    i <- which(projRecv == p)
    data_path <- dir[[i]]
    isnew <- !file.exists(file.path(data_path, 
                                    if (is.numeric(p)) 
                                      paste0("project-", p, ".motus") 
                                    else 
                                      paste0(p, ".motus")))
    db <- tagme(p, new = isnew, update = update, forceMeta = TRUE, dir = data_path)
    tmp <- db %>%
      tbl("alltags")
    if (!is.null(speciesID)) tmp <- filter(tmp, speciesID %in% !!spID[[i]])
    if (!is.null(tagProjID)) tmp <- filter(tmp, tagProjID %in% !!tpID[[i]])
    tmp <- tmp %>%
      filter(ts >= tagDeployStart,
             runLen >= minRunLen,
             # Filter out detections with missing receiver metadata, for now...
             !is.na(recvDeployID)) %>%
      select(tagProjID, fullID, spp = speciesEN, markerNumber, mfgID, motusTagID:runLen, 
             hitID, runID, batchID, ts = tsCorrected, sig, sigsd, noise, freq, freqsd, slop, burstSlop, tagBI,
             tagDeployID:tagDepLon, recvDeployID:recvDeployLon, recv, recvSiteName, 
             recvProjID, deviceID, antType, antBearing) %>%
      collect()
    
    if (nrow(tmp) == 0) {
      message("No detections as defined.")
      tmp <- NULL
    } else {
      tmp <- tmp %>%
        group_by(runID) %>%
        mutate(key = paste(motusTagID, recv, port, min(ts), sep = "_")) %>%
        ungroup() %>%
        mutate(ts = as_datetime(ts, tz = "UTC"),
               tagDeployEnd = as_datetime(tagDeployEnd, tz = "UTC"),
               tagDeployStart = as_datetime(tagDeployStart, tz = "UTC"),
               date = as.Date(ts, tz = tz),
               label = paste0(spp, ": ", mfgID, " (", motusTagID, ")"),
               lab_day_site = paste0(recvSiteName, ": ", date)) %>%
        # Drop any duplicate detections (that may vary in hitID for some unknown reason)
        distinct(motusTagID, runID, port, ts, sig, .keep_all = TRUE)
      
      # Join noise/run activity 
      tmp_runs <- db %>% tbl("runs") %>% 
        mutate(hourBin = as.integer(tsBegin/3600)) %>%
        select(runID, hourBin, motusFilter) %>%
        collect()
      tmp <- left_join(tmp, tmp_runs, by = "runID")
      tmp_activity <- db %>% tbl("activity") %>% 
        select(batchID, motusDeviceID, ant, hourBin:run7plus) %>% 
        distinct() %>% collect()
      tmp <- left_join(tmp, tmp_activity, by = c("deviceID" = "motusDeviceID", "batchID", "hourBin", "port" = "ant"))
      attr(tmp$ts, "tzone") <- attr(tmp$tagDeployStart, "tzone") <- attr(tmp$tagDeployEnd, "tzone") <- tz
    }
    tmp
  })
  out <- bind_rows(out) %>%
    mutate(recv_type = gsub("-.*$", "", recv))
  
  # Set up some hit labeling options
  prop_breaks <- c(0, 5, 10, 20, 10000)
  new_breaks <- prop_breaks[prop_breaks > minRunLen]
  brk_labs <- c(paste0(c(minRunLen, head(new_breaks, -2) + 1), "-"), ">")
  brk_labs <- paste0(brk_labs, c(head(new_breaks, -1), new_breaks[length(new_breaks)-1]))
  
  # Calculate mediate freqsd of runs
  runfreqsd <- select(out, runID, freqsd) %>%
    group_by(runID) %>%
    summarize(runfreqsd = median(freqsd))
  
  # Spin up the output
  out <- out %>%
    group_by(runID) %>% 
    arrange(ts) %>% 
    mutate(posInRun = as.integer(round((as.numeric(ts) - as.numeric(min(ts))) / tagBI)) + 1L, 
           maxburstSlop = (4 + c(-3, diff(posInRun)) - 1) / 1000,
           longRun = max(rle(cumsum(c(1, diff(posInRun) - 1)))$lengths), 
           maxRunLength = max(posInRun)) %>%
    ungroup() %>%
    left_join(runfreqsd, by = "runID") %>%
    mutate(antBin = forcats::fct_explicit_na(cut(antBearing, 
                                                 breaks = c(0, 46, 136, 226, 316, 360), right = FALSE,
                                                 labels = c("N", "E", "S", "W", "N"))),
           rlClass = cut(runLen, breaks = c(0, new_breaks),
                         labels = brk_labs),
           inspectFlag = case_when(!is.na(ambigID)               ~ "Ambiguous tag",
                                   runfreqsd > medFreqsd         ~ paste("median freqsd of run >", medFreqsd, "kHz"),
                                   !is.na(pulseCount) & pulseCount > 12000 ~ "High hourly pulse count",
                                   runLen == 3                   ~ "run length of 3",
                                   motusFilter < 1               ~ "flagged by motusFilter",
                                   recv_type %in% c("Lotek", "CTT") & longRun < 2 ~ "Lotek/CTT: <2 consecutive bursts",
                                   recv_type == "SG" & longRun < 3 ~ "SG: <3 consecutive bursts",
                                   freqsd > 0.1 & longRun < 10   ~ "hit with high freqSD",
                                   abs(burstSlop) > maxburstSlop ~ "excessive burstSlop",
                                   TRUE                          ~ NA_character_),
           inspect = !is.na(inspectFlag))
  out
}

peruse <- function(data, mt, rcv_name = "", inspect_only = TRUE) {
	runs <- filter(data, 
	               motusTagID %in% mt, 
	               grepl(rcv_name, recvSiteName, ignore.case = TRUE))
	if (inspect_only) runs <- filter(runs, inspect)
	runs <- runs %>% 
	  pull(runID) %>% unique()
	out <- filter(data, runID %in% runs) %>%
	  mutate(percRun2 = round(run2/numRuns*100, 1)) %>%
	  select(tagProjID, tagDeployID, mfgID, motusTagID, ambigID, spp, lab_day_site, 
	         ts, antBin, runID, sig:burstSlop, maxburstSlop, motusFilter:numHits, 
	         percRun2, runLen, runfreqsd, posInRun:maxRunLength, inspect, inspectFlag) %>%
	  arrange(runID, ts) %>%
	  as.data.frame()
	out
}

create_inspection_sheet <- function(data, out_dir = NULL) {
  if (is.null(out_dir)) 
    out_dir <- getwd()
  else out_dir <- normalizePath(out_dir)
  inspect_tags <- filter(data, inspect) %>% pull(motusTagID) %>% unique()
  inspect_data <- peruse(data, inspect_tags) %>%
    filter(!is.na(inspectFlag)) %>%
    arrange(motusTagID, ts) %>%
    select(mfgID, motusTagID, tagDeployID, site_date = lab_day_site, runID, antBin, 
           pulseCount, numHits, runLen, longRun, inspectFlag) %>%
    distinct() %>%
    mutate(keep = as.logical(NA))
  out_fn <- file.path(out_dir, paste0(deparse(substitute(data)), "_inspection_log.csv"))
  overwrite <- ifelse(file.exists(out_fn), readline("Output file exists. Overwrite? y/n: "), "y")
  if (identical(overwrite, "y")) 
    readr::write_csv(inspect_data, na = "", path = out_fn) else message("Output not saved.")
  return(inspect_data)
}

tag_map <- function(tagDeployID) {
  url <- paste0("https://motus.org/data/track?tagDeploymentId=", tagDeployID)
  utils::browseURL(url)
} 
