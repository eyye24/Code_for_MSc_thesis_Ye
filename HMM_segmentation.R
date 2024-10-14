#Purpose: HMM behavioural segementation. Here, modeling for 3 states
#Ellen Ye
#MSc thesis
#14.10.2024 Konstanz, GER

library(tidyverse)
library(dplyr)
library(tidyr)
library(lubridate)
library(anytime)
library(momentuHMM)
library(CircStats)
# library(circular)
# library(zoo)


##########################Running the HMM model##########################
#create vector with individuals
ind_id = as.character(unique(gps$individual_local_identifier))[-32]

gps <- gps %>% 
  drop_na(location_lat)

results <- lapply(ind_id, function(id){
  
  ind_df <- subset(gps, individual_local_identifier == id)
  ind_df <- ind_df %>% 
    filter(!is.na(location_lat))
  print(paste("Processing individual:", id))
  
  ind_df <- ind_df[!duplicated(ind_df$timestamp), ] #deal with duplicated timestamps
  
  ind_1hr = ind_df %>% 
    mutate(dt_1hr = round_date(timestamp, "1 hour")) %>% 
    group_by(dt_1hr) %>% 
    slice(1) %>%
    ungroup()
  
  ind_1hr = ind_1hr %>%
    dplyr::select(timestamp, location_long, location_lat) %>% as.data.frame()
  
  #prepare the data for momentuHMM
  prep_1hr = prepData(ind_1hr, type = "LL", 
                      coordNames = c("location_long", "location_lat")) 
  #type: LL (longitude/latitude), UTM (easting/northing), covNames = c("other variables of interest", "")
  
  # prep_1hr <- prep_1hr[complete.cases(prep_1hr[, c("step", "angle")]), ] #deal with nas
  
  set.seed(66) 
  
  stepMean0 <- c(
    m1 = mean(subset(prep_1hr, step <= 0.5)$step),
    m2 = mean(subset(prep_1hr, step > 0.5 & step <= 22)$step),
    m3 = mean(subset(prep_1hr, step > 22)$step)
  )
  
  stepSD0 <- c(
    sd1 = sd(subset(prep_1hr, step <= 0.5)$step),
    sd2 = sd(subset(prep_1hr, step > 0.5 & step <= 22)$step),
    sd3 = sd(subset(prep_1hr, step > 22)$step)
  )
  
  angleMean0 <- c(
    m1 = mean(subset(prep_1hr, angle > -0.5 & angle < 0.5)$angle),  #stationary
    m2 = mean(subset(prep_1hr, angle <= -0.5 | angle >= 0.5)$angle),  #undirected movement
    m3 = mean(subset(prep_1hr, angle > -0.5 & angle < 0.5)$angle)   #directed movement
  )
  
  angleCon0 <- c(
    rho1 = est.rho(subset(prep_1hr, angle > -0.5 & angle < 0.5)$angle),  
    rho2 = est.rho(subset(prep_1hr, angle <= -0.5 | angle >= 0.5)$angle),  
    rho3 = est.rho(subset(prep_1hr, angle > -0.5 & angle < 0.5)$angle)   
  )
  
  whichzero <- which(prep_1hr$step == 0) #check if there is step length 0, add zero mass (between 0 and 1) if so
  if(length(whichzero) > 0) {
    zero_mass <- length(whichzero) / nrow(prep_1hr) #Proportion of steps of length zero in the data set
    stepZM0 <- c(z1 = zero_mass, z2 = zero_mass, z3 = zero_mass)
    Par0 <- list(step = c(stepMean0, stepSD0, stepZM0), angle = c(angleCon0))
  } else {
    Par0 <- list(step = c(stepMean0, stepSD0), angle = c(angleCon0))
  } #select initial parameters(initial state-dependent probability distribution)
  #gamma distribution cannot handle step length 0, additional parameter needed    (https://cran.r-project.org/web/packages/moveHMM/vignettes/moveHMM-starting-values.pdf)
  
  
  print(paste("Saved initial parameters for individual:", id))
  print(c(stepMean0, stepSD0, angleCon0))
  ## ----------------------------------------------------------------------------------------------
  stateNames = c("stationary","exploration","migration")
  #for stationary: < 0.5 km
  #for exploration: < 22 km
  #for migration: > 22 km
  
  print(paste("Set state names for individual:", id))
  ## ----------------------------------------------------------------------------------------------
  dist = list(step = "gamma", angle = "vm") #choose the distribution
  
  print(paste("Set HMM parameters for individual:", id))
  print(Par0)
  
  ind_fit = fitHMM(data = prep_1hr, nbStates = 3, dist = dist, 
                   Par0 = Par0, stateNames = stateNames)
  
  print(ind_fit)
  
  ## ----------------------------------------------------------------------------------------------
  #reconstructs the most probable states sequence, using the Viterbi algorithm for a given model
  hmm_states = viterbi(ind_fit) 
  #str(hmm_states)
  
  
  ## ----------------------------------------------------------------------------------------------
  # add states to data frame
  ind_1hr$state = hmm_states
  
  return(ind_1hr)
})

#add IDs to results
names(results) <- ind_id




