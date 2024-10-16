#Purpose: Calculating angular velocity and rolling mean quat
#Ellen Ye
#MSc thesis
#14.10.2024 Konstanz, GER

library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)  # for rolling window calculations
library(patchwork) #for adding plots
library(lubridate) #for handling dates
library(fuzzyjoin)


#Put quat values in the right format
preprocess_quat <-  function(rds_files, individuals) {
  results <- list()
  
  # Loop through each individual
  for (individual in individuals) {
    cat(paste("Processing individual:", individual, "\n")) #keep comments to keep track on process
    
    # Find the file for this individual
    file_for_individual <- rds_files[grep(individual, rds_files)]
    
    if (length(file_for_individual) == 0) {
      warning(paste("No file found for individual:", individual))
      next
    }}
  
  individual_data <- readRDS(file_for_individual)
  individual_data <- individual_data %>% 
    mutate(date = as.Date(timestamp)) %>%
    filter(date <= as.Date("2023-11-30"))
  
  # Preprocess the data (quat_split step)
  cat("Preprocessing data...\n")
  quat_split <- individual_data %>%
    mutate(
      roll_split = strsplit(as.character(roll_deg), " "),
      pitch_split = strsplit(as.character(pitch_deg), " "),
      yaw_split = strsplit(as.character(yaw_deg), " ")
    ) %>%
    unnest(roll_split, pitch_split, yaw_split) %>%
    group_by(across(all_of(names(individual_data)))) %>%
    mutate(
      roll_index = row_number(),
      roll_value = as.numeric(roll_split),
      pitch_value = as.numeric(pitch_split),
      yaw_value = as.numeric(yaw_split),
      delta_t = 1 / orientation_quaternions_sampling_frequency
    ) %>%
    arrange(timestamp, roll_index) 
  
  return(quat_split)
  
}

save_split_data <- function(quat_split, individual) {
  output_dir <- "/home/ellen/Documents/Honey buzzard project/Data/IMU/QUAT/quat_split"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  file_name <- file.path(output_dir, paste0(individual, "_quat_split.rds"))
  saveRDS(quat_split, file = file_name)
  cat("Saved split data for", individual, "to", file_name, "\n")
}

process_quat <- function(rds_files, individuals) {
  for (ind in individuals) {
    # Process and split data for the individual
    quat_split <- preprocess_quat(rds_files, ind)
    
    if (!is.null(quat_split)) {
      # Save the split data
      save_split_data(quat_split, ind)
    }
  }
  cat("Processing and saving completed for all individuals.\n")
}

process_quat(rds_files, individuals)


#############calculating angular velocity###############
setwd("/home/ellen/Documents/Honey buzzard project/Data/IMU/QUAT/quat_split")
split_rds_files <- list.files(pattern = "*.rds", full.names = TRUE) #here: wd, otherwise use path=

individuals <- c("D326_193", "D329_012", "D329_013", "D329_014", "D329_015", "D326_192")

calculate_angular_velocity <- function(rds_files, ind) {
  
  cat("Processing individual:", ind, "\n")
  
  # Read the RDS file for the current individual
  file_name <- paste0(ind, "_quat_split.rds")
  file_path <- file.path(getwd(), file_name)
  
  if (!file.exists(file_path)) {
    cat("File not found for individual:", ind, "\n")
    return(NULL)
  }
  
  quat_split <- readRDS(file_path)
  quat_split <- quat_split %>% 
    filter(flight_type_sm3 %in% c("circular_soaring", "shallow_circular_soaring"))
  
  # Standardize burst windows to 8s lengths
  quat_s <- quat_split %>%
    group_by(imu_burst_id) %>%
    mutate(
      time_seconds = (row_number() - 1) / 20  # divided by sampling rate to create real time interval
    ) %>%
    ungroup()
  
  quat_8s <- quat_s %>%
    group_by(imu_burst_id) %>%
    mutate(
      segment_id = floor(time_seconds / 8), #each burst should be 8 seconds
      new_burst_id = paste(imu_burst_id, segment_id, sep = "_") #give distinct burst id
    ) %>%
    ungroup() %>%
    group_by(new_burst_id) %>%
    mutate(
      time_seconds_adjusted = time_seconds - (segment_id * 8) #time_seconds is cumulative and time_seconds_adjusted is ranges 0-8
    ) %>%
    ungroup()
  
  #calculating angular velocity for each timestamp (second) within 8 second bursts
  cat("Calculating angular velocity...\n")
  angular_velocity <- quat_8s %>%
    arrange(new_burst_id, time_seconds_adjusted) %>%
    group_by(new_burst_id, time_seconds_adjusted) %>%
    summarize(
      roll = mean(roll_split, na.rm = TRUE),
      pitch = mean(pitch_split, na.rm = TRUE),
      yaw = mean(yaw_split, na.rm = TRUE),
      individual = first(individual_local_identifier),
      timestamp = first(timestamp),
      .groups = "drop"
    ) %>%
    group_by(new_burst_id) %>%
    mutate(
      roll_velocity = (roll - lag(roll)) / (time_seconds_adjusted - lag(time_seconds_adjusted)),
      pitch_velocity = (pitch - lag(pitch)) / (time_seconds_adjusted - lag(time_seconds_adjusted)),
      yaw_velocity = (yaw - lag(yaw)) / (time_seconds_adjusted - lag(time_seconds_adjusted)),
      combined_av = sqrt(yaw_velocity^2 + pitch_velocity^2 + roll_velocity^2),
      vav_r_p = sqrt(pitch_velocity^2 + roll_velocity^2)
    ) %>%
    ungroup()
  
  
  return(angular_velocity)
}


save_av_data <- function(data, individual) {
  output_dir <- "/home/ellen/Documents/Honey buzzard project/Data/IMU/QUAT/angular_velocity"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  file_name <- file.path(output_dir, paste0(individual, "_angular_velocity.rds"))
  saveRDS(data, file = file_name)
  cat("Saved file for", individual, "to", file_name, "\n")
}


process_quat_av <- function(rds_files_split, individuals) {
  for (ind in individuals) {
    # Process and split data for the individual
    quat_av <- calculate_angular_velocity(split_rds_files, ind)
    
    if (!is.null(quat_av)) {
      # Save the split data
      save_av_data(quat_av, ind)
    }
  }
  cat("Processing and saving completed for all individuals.\n")
}

process_quat_av(split_rds_files, individuals)

#create one df and join with dates
setwd("/home/ellen/Documents/Honey buzzard project/Data/IMU/QUAT/angular_velocity")
rds_files <- list.files(pattern = "*.rds", full.names = TRUE)

av_all_df <- do.call(rbind, lapply(rds_files, readRDS))

av_all_df_dates <- av_all_df %>%
  left_join(mig_dates_df, by = c("individual" = "ind_id"))

av_all_df_life <- av_all_df_dates %>%
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


############Rolling mean and variance of quat ##################
#filter for D326_192
burst_data <- quat_split %>%
  filter(imu_burst_id == "2015")  # Replace 666 with your desired burst ID

# calculate delta_yaw and spread out points within each second
burst_data_spread <- burst_data %>%
  arrange(timestamp) %>%
  mutate(
    time_seconds = as.numeric(difftime(timestamp, min(timestamp), units = "secs")),
    delta_yaw = yaw_split - lag(yaw_split),
    delta_roll = roll_split - lag(roll_split),
    whole_second = floor(time_seconds)
  ) %>%
  group_by(whole_second) %>%
  mutate(
    spread_factor = seq(0, 0.999, length.out = n()),
    time_seconds_spread = whole_second + spread_factor
  ) %>%
  ungroup()

# calculate the rolling mean of delta_yaw
burst_data_spread <- burst_data_spread %>%
  arrange(time_seconds_spread) %>%
  mutate(rolling_mean_delta_yaw = zoo::rollmean(delta_yaw, k=10, fill=NA, align="right"),
         rolling_mean_delta_roll = zoo::rollmean(delta_roll, k=10, fill=NA, align="right"))


#include variance
burst_data_spread <- burst_data_spread %>%
  arrange(time_seconds) %>%
  mutate(
    delta_yaw = yaw_split - lag(yaw_split),
    delta_roll = roll_split - lag(roll_split),
    delta_pitch = pitch_split - lag(pitch_split),
    whole_second = floor(time_seconds)
  ) %>%
  group_by(whole_second) %>%
  mutate(
    spread_factor = seq(0, 0.999, length.out = n()),
    time_seconds_spread = whole_second + spread_factor
  ) %>%
  ungroup()

# calculate the rolling mean and variance of delta_yaw
window_size <- 10 #adjust if needed 
burst_data_spread <- burst_data_spread %>%
  arrange(time_seconds_spread) %>%
  mutate(
    rolling_mean_delta_yaw = zoo::rollmean(delta_yaw, k=window_size, fill=NA, align="right"),
    rolling_var_delta_yaw = zoo::rollapply(delta_yaw, width=window_size, FUN=var, fill=NA, align="right"),
    rolling_mean_delta_roll = zoo::rollmean(delta_roll, k=window_size, fill=NA, align="right"),
    rolling_var_delta_roll = zoo::rollapply(delta_roll, width=window_size, FUN=var, fill=NA, align="right"),
    rolling_mean_delta_pitch = zoo::rollmean(delta_pitch, k=window_size, fill=NA, align="right"),
    rolling_var_delta_pitch = zoo::rollapply(delta_pitch, width=window_size, FUN=var, fill=NA, align="right")
  )

#plotting
plot_mean_r <- ggplot(burst_data_spread, aes(x=time_seconds_spread, y=rolling_mean_delta_roll)) +
  geom_line(color="midnightblue") +
  geom_point(size = 1, alpha = 0.5, color="midnightblue") +
  theme_minimal() +
  labs(title = expression(paste("Rolling mean of ", Delta,  "roll for imu_burst_id ", 3055)),
       x = "Time (seconds)",
       y = expression(paste("Rolling mean of ", Delta,  "roll"))) +
  scale_x_continuous(breaks = seq(0, max(burst_data_spread$time_seconds_spread), by = 60))

plot_var_r <- ggplot(burst_data_spread, aes(x=time_seconds_spread, y=rolling_var_delta_roll)) +
  geom_line(color="cornflowerblue") +
  geom_point(size = 1, alpha = 0.5, color="cornflowerblue") +
  theme_minimal() +
  labs(title = expression(paste("Rolling variance of ", Delta, "roll for imu_burst_id ", 3055)),
       x = "Time (seconds)",
       y = expression(paste("Rolling variance of ", Delta, " roll"))) +
  scale_x_continuous(breaks = seq(0, max(burst_data_spread$time_seconds_spread), by = 60))

plot_mean_y <- ggplot(burst_data_spread, aes(x=time_seconds_spread, y=rolling_mean_delta_yaw)) +
  geom_line(color="gold2") +
  geom_point(size = 1, alpha = 0.5, color="gold2") +
  theme_minimal() +
  labs(title = expression(paste("Rolling mean of ", Delta,  "yaw for imu_burst_id ", 3055)),
       x = "Time (seconds)",
       y = expression(paste("Rolling mean of ", Delta,  "yaw"))) +
  scale_x_continuous(breaks = seq(0, max(burst_data_spread$time_seconds_spread), by = 60))

plot_var_y <- ggplot(burst_data_spread, aes(x=time_seconds_spread, y=rolling_var_delta_yaw)) +
  geom_line(color="goldenrod1") +
  geom_point(size = 1, alpha = 0.5, color="goldenrod1") +
  theme_minimal() +
  labs(title = expression(paste("Rolling variance of ", Delta, "yaw for imu_burst_id ", 3055)),
       x = "Time (seconds)",
       y = expression(paste("Rolling variance of ", Delta, " yaw"))) +
  scale_x_continuous(breaks = seq(0, max(burst_data_spread$time_seconds_spread), by = 60))

plot_mean_p <- ggplot(burst_data_spread, aes(x=time_seconds_spread, y=rolling_mean_delta_pitch)) +
  geom_line(color="red") +
  geom_point(size = 1, alpha = 0.5, color="red") +
  theme_minimal() +
  labs(title = expression(paste("Rolling mean of ", Delta,  "pitch for imu_burst_id ", 3055)),
       x = "Time (seconds)",
       y = expression(paste("Rolling mean of ", Delta,  "pitch"))) +
  scale_x_continuous(breaks = seq(0, max(burst_data_spread$time_seconds_spread), by = 60))

plot_var_p <- ggplot(burst_data_spread, aes(x=time_seconds_spread, y=rolling_var_delta_pitch)) +
  geom_line(color="firebrick1") +
  geom_point(size = 1, alpha = 0.5, color="firebrick1") +
  theme_minimal() +
  labs(title = expression(paste("Rolling variance of ", Delta, "pitch for imu_burst_id ", 3055)),
       x = "Time (seconds)",
       y = expression(paste("Rolling variance of ", Delta, " pitch"))) +
  scale_x_continuous(breaks = seq(0, max(burst_data_spread$time_seconds_spread), by = 60))

# Combine plots
combined_plot_quat <- plot_mean_y/plot_var_r/plot_var_p/plot_vs_av/acc_plot