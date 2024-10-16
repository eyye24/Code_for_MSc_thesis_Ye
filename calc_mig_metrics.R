#Purpose: Calculating migration metrics
#Ellen Ye
#MSc thesis
#14.10.2024 Konstanz, GER


library(geosphere)
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(lme4)
library(performance)

#subsampling for 1 hr gps 
#gps: all_gps_apr15_24.rds

gps <-  gps %>% 
  mutate(date = as.Date(timestamp))

gps_sub <- gps %>% 
  group_by(individual_local_identifier, date) %>% 
  mutate(ind_ID_day = paste0(individual_local_identifier, "_", yday(as.Date(timestamp))),
         dt_hr = round_date(timestamp, "hour")) %>%
  ungroup()

#filter data for individual migration periods
ind_id <- unique(as.character(gps$individual_local_identifier))[-32]

subset_ind_mig <- function(ind_id, start_date, end_date) {
  ind_data <- gps_sub %>% 
    filter(individual_local_identifier == ind_id,
           date >= start_date,
           (is.na(end_date) | date <= end_date))
  
  
  return(ind_data)
}


gps_mig_ind <- mig_dates_df %>%
  rowwise() %>%  #each row on its own
  mutate(data = list(subset_ind_mig(ind_id, migration_start, migration_end))) %>%  #subset_ind_mig function
  ungroup() %>%  # Ungroup after the operation
  select(ind_id, data) 

gps_mig <- as.data.frame(bind_rows(gps_mig_ind$data))

#check date ranges
date_summary <- gps_mig %>%
  group_by(individual_local_identifier) %>%
  summarise(
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    n_observations = n()
  ) %>%
  arrange(individual_local_identifier)

print(date_summary, n = Inf)


#############calculating distance, speed####################
gps_mig <- gps_mig %>%
  filter(location_lat != 0 & location_long != 0)

#filter out gappy or dead  inds
#dead or gappy
dead <-  c("D299_270", "D225_236", "D326_193", "D329_015", "D324_513")
gps_sub <- gps_sub %>% 
  filter(!(individual_local_identifier%in%dead))

gps_sub <- gps_sub %>%
  filter(!is.na(location_lat) & !is.na(location_long))

gps_mig <-  gps_sub %>%
  mutate(hour = floor_date(timestamp, "hour"))

gps_sub <- gps_mig %>% 
  arrange(individual_local_identifier, hour) %>%
  distinct(individual_local_identifier, hour, .keep_all = TRUE)

gps_dis <- gps_sub %>%
  mutate(
    lat_next = lead(location_lat), #shift for estimating pairwise distance
    long_next = lead(location_long),
    
    # Use distHaversine to calculate distance between each point and the next one (in meters)
    distance_m = distHaversine(cbind(location_long, location_lat), 
                               cbind(long_next, lat_next)),
    distance_km = distance_m / 1000
  )

gps_dis <- gps_dis %>% 
  filter(distance_km <= 1000)

gps_dis_sp <-  gps_dis %>% 
  group_by(individual_local_identifier) %>%
  mutate(
    time_diff_hours = as.numeric(difftime(timestamp, lag(timestamp), units = "hours")),
    speed_kmh = ifelse(time_diff_hours > 0, distance_km / time_diff_hours, NA)
  ) %>%
  ungroup()

mig_sum <- gps_dis_sp %>%
  mutate(date = as.Date(timestamp)) %>% 
  group_by(individual_local_identifier) %>%
  summarize(
    total_distance_km = sum(distance_km, na.rm = TRUE),
    avg_speed_kmh = mean(speed_kmh, na.rm = TRUE),
    total_duration_days = round(as.numeric(difftime(max(timestamp), min(timestamp), units = "days"))),
    overall_speed_kmh = total_distance_km / total_duration_days
  ) %>%
  arrange(desc(total_distance_km))


#calculate daily values

# calculate daily VeDBA
daily_vedba <- function(df) {
  df %>%
    mutate(date = as.Date(timestamp)) %>%
    group_by(individual_local_identifier, date) %>%
    summarize(
      daily_vedba = mean(VeDBA, na.rm = TRUE),
      .groups = 'drop'
    )
}

daily_vedba_df <- map_dfr(DBA_mig, daily_vedba) 



#add daily distance and total migration time
daily_sum_dis_sp <- gps_dis_sp %>%
  group_by(individual_local_identifier, date) %>%
  summarize(
    daily_distance_km = sum(distance_km, na.rm = TRUE),
    daily_speed_kmh = mean(speed_kmh, na.rm = TRUE),
    .groups = 'drop'
  )


daily_wobble_mig <- daily_vedba_df %>%
  left_join(daily_sum_dis_sp, by = c("individual_local_identifier", "date")) %>%
  left_join(wobble_daily_data, by = c("individual_local_identifier", "date"))


daily_wobble_mig <-  daily_wobble_mig %>% 
  drop_na(c(daily_distance_km, daily_vedba))

