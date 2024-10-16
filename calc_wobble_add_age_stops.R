library(tidyverse)
library(dplyr)
library(tidyr)

#Purpose: Calculating wobble and adding age and stopovers
#Ellen Ye
#MSc thesis
#14.10.2024 Konstanz, GER


#set directory to folder with classified bursts (Scacco et al. 2019)
setwd("/home/ellen/Documents/Honey buzzard project/Data/GPS/GPS segmentation/gps_seg_apr24")
rds_files <- list.files(pattern = "*.rds", full.names = TRUE)

calculate_wobble <- function(file) {
  
  ind <- readRDS(file)
  
  ind_gps_df <- ind %>% 
    dplyr::select(timestamp, individual_local_identifier, burst_id, altitude_diff, 
                  vert_speed, turn_angle, vert_speed_smooth, turn_angle_cum, flight_clust_sm3) 
 
  
  #group by individual and burst_id to calculate sd, mean 
  ind_wobble <- ind_gps_df %>%
    group_by(individual_local_identifier, burst_id) %>% 
    summarise(
      sd_vert_speed = sd(vert_speed, na.rm = TRUE),
      sd_turn_angle = sd(turn_angle, na.rm = TRUE),
      mean_vert_speed = mean(vert_speed, na.rm = TRUE),
      mean_turn_angle = mean(turn_angle, na.rm = TRUE),
      flight_class3 = first(flight_clust_sm3),
      date = as.Date(first(timestamp)),
      timestamp = first(timestamp),
      .groups = 'drop'
    ) 
  
  
  
  ind_wobble
}

all_ind_wobble <- map_dfr(rds_files, calculate_wobble) #instead of lapply to keep df format


#add migration dates and life stage based on the dates
all_ind_wobble_mig <- all_ind_wobble %>%
  left_join(mig_dates_df, by = c("individual_local_identifier" = "ind_id")) #mig_dates_df from 1.Data_extraction_wHMM.R


#post-fledging and migration
all_ind_wobble_life <- all_ind_wobble_mig %>%
  mutate(
    life_stage = case_when(
      # Dependency period: from tagging until first exploration
      date < first_exploration ~ "dependency_period",
      #exploration period: period between first exploration and migration start
      date >= first_exploration & date < migration_start ~ "exploration",
      date >= migration_start & date <= migration_end ~ "migration", 
      date >= migration_start & date <= is.na(migration_end) ~ "migration",
      TRUE ~ "wintering"
    )
  )

all_ind_wobble_life_mig <-  all_ind_wobble_life %>% 
  filter(life_stage != "wintering")


#calculate circling ratio per day
circling_ratio <- all_ind_wobble_life_mig %>%
  mutate(date = as.Date(timestamp)) %>%
  group_by(date) %>%
  summarize(
    circular_count = sum(flight_class3 == "circular_soaring", na.rm = TRUE),
    shallow_circular_count = sum(flight_class3 == "shallow_circular_soaring", na.rm = TRUE),
    total_count = sum(!is.na(flight_class3)),
    circling_ratio = (circular_count + shallow_circular_count) / total_count
  ) %>%
  ungroup()

all_ind_wobble_life <- all_ind_wobble_life %>%
  left_join(select(circling_ratio, date, circling_ratio), by = "date")



#filter for circling
all_ind_wobble_circ <- all_ind_wobble_life %>% 
  filter(flight_class3 %in% c("circular_soaring", "shallow_circular_soaring"))


#Adding estimated age with average wing growth rate 8.825 mm/day
#meta data: EHB_metadata.csv
meta <- meta %>% 
  mutate(est_age_deployment = wing_length_mm / 8.825)

all_ind_wobble_circ_age <- all_ind_wobble_circ %>%
              left_join(meta %>% 
                       dplyr::select(ring_ID, sex, deployment_dt_utc, wing_length_mm, est_age_deployment),
                       by = c("individual_local_identifier" = "ring_ID"))


all_ind_wobble_circ_age <- all_ind_wobble_circ_age %>%
  mutate(
    deployment_date = as.Date(deployment_dt_utc),
    est_age = est_age_deployment + as.numeric(date - deployment_date)
  )

all_ind_wobble_circ_age<- all_ind_wobble_circ_age %>%
        mutate(est_age = round(est_age))



#calculate stopover duration for each individual in each stage
#stopovers_df from script "Date_ectraction_wHMM"

calculate_stopover_duration <- function(stopovers, daily_data_circ) {
  range_of_interest <- range(daily_data_circ$date)
  
  # prep stopovers data
  stopovers_prep <- stopovers %>%
    mutate(start_date = as.Date(start),
           end_date = as.Date(end)) %>%
    filter(start_date >= range_of_interest[1] & start_date <= range_of_interest[2])
  
  #identify life stage of animal to identify exploration stops of migration stopovers
  get_life_stage <- function(date, individual_local_identifier) {
    stage <- daily_data_circ %>%
      filter(individual_local_identifier == !!individual_local_identifier, date == !!date) %>%
      pull(life_stage)
    
    if (length(stage) == 0) {
      return(NA_character_)  # Return NA if no matching row is found
    } else {
      return(stage[1])
    }
  }
  
  #calculating durations for each type of stopover
  result <- stopovers_prep %>%
    rowwise() %>%
    mutate(life_stage = get_life_stage(start_date, individual_local_identifier)) %>%
    group_by(individual_local_identifier) %>%
    summarise(
      pre_mig_stopover_duration = round(sum(duration[life_stage %in% c("dependence_period", "exploration")], na.rm = TRUE)),
      mig_stopover_duration = round(sum(duration[life_stage == "migration"], na.rm = TRUE)),
      total_stopover_duration = round(sum(duration, na.rm = TRUE))
    )
  
  return(result)
}


all_stopovers <- calculate_stopover_duration(stopovers_df, all_ind_wobble_circ_age)


all_ind_wobble_circ_compl <- all_ind_wobble_circ_age %>%
  left_join(all_stopovers %>%
              dplyr::select(individual_local_identifier,
                            pre_mig_stop_dur = pre_mig_stopover_duration,
                            mig_stop_dur = mig_stopover_duration,
                            total_stop_dur = total_stopover_duration),
            by = "individual_local_identifier")



#create daily wobble
daily_data <- all_ind_wobble_circ_compl %>%
  mutate(date = as.Date(date)) %>%
  group_by(individual_local_identifier,date, burst_id) %>%
  summarise(
    avg_sd_vert_speed = mean(sd_vert_speed, na.rm = TRUE),
    avg_sd_turn_angle = mean(sd_turn_angle, na.rm = TRUE),
    age = first(est_age),
    life_stage = first(life_stage),
    circling_ratio = first(circling_ratio),
    flight_type = first(flight_class3),
    pre_mig_stop_dur = first(pre_mig_stop_dur),
    mig_stop_dur = first(mig_stop_dur),
    total_stop_dur = first(total_stop_dur),
    .groups = "drop"
  )







