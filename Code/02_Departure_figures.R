source("Code/01_Calculate_parameters_of_interest.R")

rutu_dep <- readr::read_csv("Data/Derived/rutu_departures.csv",
                            col_types = "icciTddDDc") %>%
  filter(!is.na(min_depart)) %>%
  mutate(tdiff = yday(max_depart) - yday(min_depart) + 1)

# ANy association between capture weight and how long it took to detect birds elsewhere
# in the Motus network?
rutu_dep <- left_join(rutu_dep, select(rutu_deploy, motusTagID, wt_g))
with(rutu_dep, plot(wt_g, tdiff))
with(rutu_dep, cor.test(wt_g, tdiff, method = "spearman"))

# Two approaches to SC departure date figure
# 1 - use only departures with a precisely known departure date
# 2 - use all departures, but construct density that reflects our uncertainty in
#     that each individual contributes equally, but dates of departure for those
#     individuals with a possible range of dates only contribute fractionally to
#     the constructed "density" (i.e., each possible date contributes only
#     [1 / # possible days of departure])
# In both scenarios, we summarize by week of departure since data are fairly sparse

days_of_year <- seq.Date(from = as.Date("2018-01-01"), to = as.Date("2018-12-31"), by = 1) %>%
  format(format = "%d %b")
doy_labels <- format(days_of_year, format = "%d %b")
wk_starts <- seq.Date(from = as.Date("2015-01-01"), by = 7, length.out = 52)
wk_breaks <- week(wk_starts)
wk_labels <- format(wk_starts, format = "%d %b")
label_days  <- function(x, nth = 1) every_nth(doy_labels[x], nth, inverse = TRUE)
label_weeks  <- function(x, nth = 1) every_nth(wk_labels[x], nth, inverse = TRUE)
make_breaks <- function(x, by, ...) seq(floor(x[1]), ceiling(x[2]), by = by)

# Range of all possible departures for plotting later
dep_range <- range(yday(c(pull(rutu_dep, min_depart), pull(rutu_dep, max_depart))))

# Approach 1
dep_days <- filter(rutu_dep, tdiff == 1) %>%
  mutate(dep_doy = yday(max_depart),
         dep_wk = week(max_depart))

# calculate median departure date for display
dep_days_sum <- dep_days %>%
  summarize(n = n(),
            med_dep_doy = median(dep_doy),
            med_dep_wk = med_dep_doy / 7 + 1,
            med_dep_date = days_of_year[round(med_dep_doy)])

dep_fig_precise <- 
  ggplot(dep_days, aes(dep_doy)) +
  geom_vline(xintercept = dep_days_sum$med_dep_doy, 
               color = "gray50", lwd = 4) + 
  # geom_text(data = dep_all_sum, aes(x = wt_med_dep_doy + .5, y = max(dep_all_wt$dep_doy_dens) + 0.025,
  #                                   label = paste("Median:", wt_med_dep_date)),
  #           size = 6, hjust = 0) +
  geom_text(data = dep_days_sum, aes(x = med_dep_doy + .5, y = 5, 
                                      label = paste("Median:", med_dep_date)), 
            size = 6, hjust = 0) +
  geom_bar(fill = "grey75", color = "black") +
  geom_text(data = dep_days_sum, aes(x = min(dep_range) - 1, y = 5, label = paste("n =", n)), 
            hjust = 0, size = 6) +
  scale_x_continuous("Date of departure", breaks = make_breaks(dep_range + c(-1, 1), by = 3), 
                     minor_breaks = make_breaks(dep_range + c(-1, 1), by = 1), limits = dep_range + c(-1, 1),
                     labels = label_days) +
  labs(y = "# Ruddy Turnstones") +
  theme_bw(base_size = 16)

# Approach 2
dep_all <- rutu_dep %>%
  mutate(date_rng = purrr::pmap(., function(min_depart, max_depart, ...) {
    seq.Date(min_depart, max_depart, by = 1)
  })) %>%
  tidyr::unnest(cols = c(date_rng)) %>%
  mutate(dt_wt = 1 / tdiff,
         doy = yday(date_rng))

# calculate weighted median departure date by banding location for display
dep_all_sum <- dep_all %>%
  summarize(n = n_distinct(motusTagID),
            wt_med_dep_doy = median_wtd(doy, dt_wt),
            wt_med_dep_wk = wt_med_dep_doy / 7 + 1,
            wt_med_dep_date = days_of_year[round(wt_med_dep_doy)])

dep_all_wt <- dep_all %>%
  mutate(dep_wk = week(date_rng)) %>%
  group_by(doy) %>%
  summarize(date_wt = sum(dt_wt))

dep_fig_wt_spr <- ggplot(dep_all_wt, aes(doy)) +
  geom_vline(xintercept = dep_all_sum$wt_med_dep_doy,
             color = "grey50", lwd = 4) +
  geom_text(data = dep_all_sum, aes(x = wt_med_dep_doy + .5, 
                                    y = max(dep_all_wt$date_wt) + 0.5,
                                    label = paste("Median:", wt_med_dep_date)),
            size = 6, hjust = 0) +
  geom_bar(aes(y = date_wt), stat = "identity", fill = "grey75", color = "black") +
  geom_text(data = dep_all_sum, aes(x = min(dep_range) - 1,
                                     y = max(dep_all_wt$date_wt) + 0.5, 
                                        label = paste("n =", n)), 
            hjust = 0, size = 6) +
  scale_x_continuous("Date of departure", breaks = make_breaks(dep_range + c(-1, 1), by = 3), 
                     minor_breaks = make_breaks(dep_range + c(-1, 1), by = 1), limits = dep_range + c(-1, 1),
                     labels = label_days) +
  labs(y = "# Ruddy Turnstones") +
  theme_bw(base_size = 16)

# Assemble
all <- dep_fig_precise | dep_fig_wt_spr
all + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 22))
ggsave("Output/rutu_departures.png", dpi = 600, height = 5, width = 10)
