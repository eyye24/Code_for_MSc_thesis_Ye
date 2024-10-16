#Purpose: Extracting migration, exploration and stopover dates from HMM results 
#Ellen Ye
#MSc thesis
#14.10.2024 Konstanz, GER

library(tidyverse)
library(dplyr)
library(tidyr)


##########################Selecting migration time with HMM results##########################
#Prepare data for migration date extraction, group by day
dt_transform <- function(data) {
  data$date <- as.Date(data$timestamp)
  
  sum_state <-  data %>%
    mutate(date = as.Date(timestamp)) %>%
    group_by(date, state) %>%
    summarise(count = n(),
              location_lat = mean(location_lat, na.rm = TRUE)) %>% 
    ungroup()
  
  return(sum_state)
}
#Parameters used in mig_time function: objects keeping track of counts for conditions
state_3_count <- 0
consecutive_days <- 0
migration_start <- NA

# Find migration start
mig_time <- function(data) {
  data = dt_transform(data)
  
  #condition 1: select the date for state 3 occurrences
  potential_starts <- data %>% 
    filter(state == 3) %>%
    arrange(date)
  
  #condition 2: check if there's at least one state 3 in the following 3 days
  check_following_days <- function(start_date) {
    following_days <- seq(start_date + 1, by = "day", length.out = 3) #length.out: desired length of the sequence; by: increment of the sequence
    any(data$date %in% following_days & data$state == 3)
  }
  
  #find first date that has both conditions: 
  migration_start <- NA
  for (i in 1:nrow(potential_starts)) {
    if (check_following_days(potential_starts$date[i])) { #if true it is -> migration_start
      migration_start <- potential_starts$date[i]
      break
    }
  }
  
  #if no suitable date found, return NA to avoid errors
  if (is.na(migration_start)) {
    return(list(
      migration_start = NA,
      migration_end = NA
    ))
  }
  
  
  
  #add column with days since migration to select for migration start
  data <- data %>%
    mutate(days_since_mig = as.integer(date - migration_start))
  
  data_pmig <- data %>% filter(date > as.Date(migration_start))
  
  #create variables for counting
  last_state_3_day <- 0
  migration_end <- NA
  
  #migration_end
  for (i in 1:nrow(data_pmig)) {
    current_day <- data_pmig$days_since_mig[i] 
    lat <- data_pmig$location_lat[i]
    
    if (data_pmig$state[i] == 3) {
      last_state_3_day <- current_day #if state is 3, counter gets assigned current day
    } else if (current_day - last_state_3_day >= 3 && !is.na(lat) && lat <= 16) {  #check if the difference between the last migration day and the current date surpass 3 days
      migration_end <- data_pmig$date[i]
      break
    }
  }
  
  return(list(
    migration_start = migration_start,
    migration_end = migration_end 
  ))
}

mig_dates <- lapply(results, mig_time)

mig_dates_df <- do.call(rbind, lapply(names(mig_dates), function(id) {
  data.frame(
    ind_id = id,
    migration_start = mig_dates[[id]]$migration_start,
    migration_end = mig_dates[[id]]$migration_end,
    stringsAsFactors = FALSE
  )
}))


##########################Selecting exploration time with HMM results##########################
#function selects the first occurrence of a state
first_occ <- function(df) {
  df %>%
    mutate(date = as.Date(timestamp)) %>% 
    filter(state == "2") %>%
    slice_head(n = 1) %>% #select for first row
    pull(date)
}

expl_date <- lapply(results, first_occ)

#add to df
expl_vector <- unlist(expl_date)
expl_vector <- as.Date(expl_vector)

mig_dates_df$first_exploration <- expl_vector[match(mig_dates_df$ind_id, names(expl_vector))]


##########################Selecting stopovers##########################
id_stopovers <- function(data, lat_threshold = 0.1, min_duration = 2) { # 1 deg lat is approx 111 km; duration in days
  stopovers <- list()
  in_stopover <- FALSE
  stopover_start <- NULL
  last_lat <- NULL
  
  for (i in 1:(nrow(data) - 1)) {
    current_state <- data$state[i]
    next_state <- data$state[i + 1]
    current_timestamp <- data$timestamp[i]
    current_lat <- data$location_lat[i]
    
    if (is.na(current_state) || is.na(next_state) || is.na(current_lat)) {
      print(paste("NA values detected at row", i))
      next
    }
    
    #Condition 1 & 2: set potential start of a potential stopover (selecting starting point)
    if (!in_stopover && (current_state == 2 || current_state == 3) && next_state == 1) {
      stopover_start <- current_timestamp
      last_lat <- current_lat
      in_stopover <- TRUE
    }
    
    #condition 3: continue stopover if in state 1 or 2 with small latitude change
    if (in_stopover && (current_state == 1 || current_state == 2) && abs(current_lat - last_lat) < lat_threshold) {
      last_lat <- current_lat
    }
    
    #end of stopover: state 3 or significant latitude change
    if (in_stopover && ((current_state == 3) || (abs(current_lat - last_lat) >= lat_threshold))) {
      stopover_end <- current_timestamp
      duration <- as.numeric(difftime(stopover_end, stopover_start, units = "days"))
      
      #only adding stopover if duration is greater than minimum and there's a change in latitude
      if (duration > min_duration && abs(current_lat - last_lat) > 0) {
        stopovers <- c(stopovers, list(list(
          start = stopover_start,
          end = stopover_end,
          duration = duration,
          start_lat = last_lat,
          end_lat = current_lat
        )))
      }
      
      in_stopover <- FALSE
      stopover_start <- NULL
      last_lat <- NULL
    }
  }
  
  #handle case where stopover continues until the end of the data
  if (in_stopover) {
    stopover_end <- data$timestamp[nrow(data)]
    duration <- as.numeric(difftime(stopover_end, stopover_start, units = "days"))
    end_lat <- data$location_lat[nrow(data)]
    
    # Only add final stopover if duration is greater than minimum and there's a change in latitude
    if (duration > min_duration && abs(end_lat - last_lat) > 0) {
      stopovers <- c(stopovers, list(list(
        start = stopover_start,
        end = stopover_end,
        duration = duration,
        start_lat = last_lat,
        end_lat = end_lat
      )))
    }
  }
  
  return(stopovers)
}

#for the entire list
stopovers <- lapply(names(results), function(id) {
  individual_data <- results[[id]]
  stopovers <- id_stopovers(individual_data, lat_threshold = 0.1, min_duration = 1)  #set latitude threshold, and minimum stay (day)
  
  #if stopovers found, convert to df and ID
  if (length(stopovers) > 0) {
    stopover_df <- do.call(rbind, lapply(stopovers, as.data.frame))
    stopover_df$ind_id <- id
    return(stopover_df)
  } else {
    return(NULL)
  }
})

names(stopovers) <- names(results)
stopovers_df <- bind_rows(stopovers)


