#Purpose: Calculating DBA 
#Ellen Ye
#MSc thesis
#14.10.2024 Konstanz, GER


library(dplyr)
library(tidyr)
library(purrr)

setwd("/home/ellen/Documents/Honey buzzard project/Data/IMU/ACC/matched_gps_acc")
csv_files <- list.files(pattern = "*.csv", full.names = TRUE) #here: wd, otherwise use path=



###prepare data for later filtering ### 
extract_ring_id <- function(file_name) {
  parts <- strsplit(basename(file_name), "_")[[1]]
  paste(parts[1], parts[2], sep = "_")
}

assign_life_stage <- function(files) {
  
  data <- read.csv(files, stringsAsFactors = FALSE)
  ring_id <- extract_ring_id(files)
  
  id_data <- data[data$individual_local_identifier == ring_id,]
  id <- unique(id_data$individual_local_identifier)
  
  deployment_t <-  meta[meta$ring_ID == id,]
  mig_start <- as.Date(mig_dates[[id]]$migration_start)
  mig_end <- as.Date(mig_dates[[id]]$migration_end)
  
  if (is.na(mig_end) || is.na(mig_start) || length(deployment_t) == 0 || is.na(deployment_t)) {
    print(paste("Skipping", ring_id))
    return(NULL)
  }
  
  
  ind_data <- data %>%
    mutate(
      date = as.Date(timestamp),
      life_stage = ifelse(  #ifelse(test, yes, no)
        timestamp < as.Date(deployment_t$deployment_dt_utc), "pre_deployment",
        ifelse(timestamp >= as.Date(deployment_t$deployment_dt_utc) &
                 timestamp < as.Date(mig_dates[[deployment_t$ring_ID]]$migration_start), "pre_migration",
               ifelse(timestamp >= as.Date(mig_dates[[deployment_t$ring_ID]]$migration_start) &
                        timestamp < as.Date(mig_dates[[deployment_t$ring_ID]]$migration_end), "migration",
                      ifelse(timestamp > as.Date(mig_dates[[deployment_t$ring_ID]]$migration_end), "post_migration",
                             "unknown"
                      )
               )
        )
      )
    )
  
  return(ind_data)
  
}

acc_life <- lapply(csv_files, assign_life_stage)
acc_life <- Filter(Negate(is.null), acc_life) #remove empty lists (base)

names(acc_life) <- sapply(acc_life, function(df) unique(df$individual_local_identifier)[1])


###################### Calculate DBA ###################### 
#function to calculate VeDBA and ODBA
DBA <- function(data) {
  MyDBA <- function(acc_string, sampling_rate = 24, axes = 3) {
    ACC_values <- as.numeric(unlist(strsplit(acc_string, " ")))
    complete_samples <- floor(length(ACC_values) / (sampling_rate * axes))
    ACC1 <- matrix(ACC_values[1:(complete_samples * sampling_rate * axes)], 
                   ncol = axes, byrow = TRUE)
    VeDBA <- mean(sqrt(rowSums((ACC1 - colMeans(ACC1))^2)))
    ODBA <- mean(rowSums(abs(ACC1 - colMeans(ACC1))))
    return(c(VeDBA = VeDBA, ODBA = ODBA))
  }
  
  # Apply MyDBA to each row of the data frame
  DBA_results <- t(sapply(data$eobs_acceleration_g, MyDBA))
  
  # Add results to the original data frame
  data$VeDBA <- DBA_results[, "VeDBA"]
  data$ODBA <- DBA_results[, "ODBA"]
  
  return(data)
}

# Apply DBA function to each data frame in ind_data_ls
all_DBA <- lapply(acc_life, function(df) {
  tryCatch({
    DBA(df)
  }, error = function(e) {
    warning(paste("Error processing data frame:", "\nError:", e$message))
    return(NULL)
  })
})

names(all_DBA) <- names(acc_life) #add existing names of list elements

all_DBA <- lapply(all_DBA, function(df) {
  df %>% mutate(date = as.Date(timestamp))
})

DBA_mig <- bind_rows(DBA_mig, .id = "individual_local_identifier")


DBA_mig <- all_DBA %>%
  filter(life_stage == "migration")


#plotting quat comparison
match_timestamps <- function(acc_df, reference_df, time_column = "timestamp", tolerance = 1) {
  
  acc_df[[time_column]] <- as.POSIXct(acc_df[[time_column]], tz = "UTC")
  reference_df[[time_column]] <- as.POSIXct(reference_df[[time_column]], tz = "UTC")
  
  #find time rrange
  start_time <- min(reference_df[[time_column]]) - seconds(tolerance)
  end_time <- max(reference_df[[time_column]]) + seconds(tolerance)
  
  # Subset acc_df to the relevant time range
  acc_subset <- acc_df %>% 
    filter(!!sym(time_column) >= start_time & !!sym(time_column) <= end_time)
  
  # Function to find matches within tolerance
  find_matches <- function(x, y, tol) {
    abs(difftime(x, y, units = "secs")) <= tol
  }
  
  
  matched_list <- list()
  
  #loop through timestamp
  for (i in seq_along(acc_subset[[time_column]])) {
    acc_time <- acc_subset[[time_column]][i]
    
    # Find all matches in reference_df within the tolerance
    ref_matches <- reference_df %>% 
      filter(find_matches(!!sym(time_column), acc_time, tolerance))
    
    #Append rows
    if (nrow(ref_matches) > 0) {
      #don't forget renaming before combining
      ref_matches <- ref_matches %>%
        rename_with(~paste0("ref_", .), all_of(time_column))
      
      #combining rows from acc_subset and reference matches
      matched_rows <- cbind(acc_subset[i, ], ref_matches)
      matched_list[[length(matched_list) + 1]] <- matched_rows
    }
  }
  
  matched_data <- do.call(rbind, matched_list)
  
  return(matched_data)
}



acc_burst <- match_timestamps(acc_df, qb_spread)

# create precise time 
acc_burst <- acc_burst %>%
  arrange(timestamp) %>%
  mutate(
    start_time = min(timestamp),
    elapsed_time = as.numeric(difftime(timestamp, start_time, units = "secs")),
    row_within_second = row_number() - 1,
    precise_time = elapsed_time + (row_within_second / (n() / n_distinct(timestamp)))
  )

# scaling
time_scale_factor <- max(acc_burst$time_seconds_spread) / max(acc_burst$precise_time)

#scale time so it shows 240 s range (instead of 400))
acc_burst <- acc_burst %>%
  mutate(scaled_precise_time = precise_time * time_scale_factor)


#reshape for plotting
acc_long <- acc_burst %>%
  dplyr::select(scaled_precise_time, surge_g, sway_g, heave_g, time_seconds_spread) %>%
  pivot_longer(cols = c(surge_g, sway_g, heave_g), 
               names_to = "acceleration_type", 
               values_to = "acceleration")

#plotting
acc_plot <- ggplot(acc_long, aes(x = scaled_precise_time, y = acceleration, color = acceleration_type)) +
  geom_line(size = 0.5) +
  labs(
    x = "Time (seconds)",
    y = "Acceleration (g)",
    color = "Axes") +
  theme_classic() +
  scale_color_manual(values = c("surge_g" = "#D95319", "sway_g" = "#0072BD", "heave_g" = "#EDB120"),
                     labels = c("Heave", "Sway", "Surge")) +
  theme(legend.position = "bottom", 
        text = element_text(size = 16),
        axis.title = element_text(size = 20),  # Axis titles
        axis.text = element_text(size = 18),   
        plot.subtitle = element_text(size = 20),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 17),
        plot.margin = margin(r = 20))
