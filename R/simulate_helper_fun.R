#This script simulates an outbreak from when the first case of each lineage was 
#introduced to when the lineage was detected to give an estimate of how many cases were 
#circulating on Romblon in this time. 

# #R0 parameter values from Townsend et al., 2013
# R0 = 1.20
# dispersion = 1.33
# #SI parameter values from Mancy et al., 2022
# SI_meanlog <- 2.85
# SI_sdlog <-  0.966

# Function to generate secondary cases
infect <- function(IDs, today, R0, dispersion, SI_meanlog, SI_sdlog){
  offspring <- rnbinom(length(IDs), mu = R0, size = dispersion)
  offspringTotal <- sum(offspring)
  if(offspringTotal > 0){
    for(j in 1:offspringTotal){
      caseID <- 1:offspringTotal
      parentID <- rep(IDs, offspring)
      infD <- rep(today, offspringTotal)
      transD <- today + rlnorm(offspringTotal, meanlog = SI_meanlog, sdlog = SI_sdlog)
    }
  } else { 
    parentID <- caseID <- infD <- transD <- c()
  }
  return(data.frame(caseID, parentID, infD, transD))
}
# TEST
# test <- infect(c(5,12,13), 25, R0, dispersion, SI_meanlog, SI_sdlog); nrow(test)

# Function to simulate outbreaks
outbreak_sim <- function(R0, dispersion, SI_meanlog, SI_sdlog, lastday) {
  
  # Set up vectors to track events
  outbreak_df <- data.frame(caseID = 1, parentID = 0, infD = 0, transD = 1) # set up outbreak dataframe
  today <- 0
  
  # loop thru every day
  for(i in 1:lastday){ 
    today <- today + 1 
    cases <- outbreak_df[which(round(outbreak_df$transD) == today),]
    
    # If there are cases today generate new infections
    if(nrow(cases)>0){ 
      newCases <- infect(cases$caseID, today, R0, dispersion, SI_meanlog, SI_sdlog)
      # If there are new cases the newCases$ID columnn needs updating
      if(nrow(newCases > 0)){ newCases$caseID <-  max(outbreak_df$caseID) + (1:nrow(newCases)) } 
    } else { newCases <- data.frame(caseID = c(), parentID = c(), infD = c(), transD = c()) } # Make sure to replace newCases with an empty dataframe!

    outbreak_df <- rbind(outbreak_df, newCases) # update outbreak dataframe
  }
  return(outbreak_df)
}
# TEST
# test <- outbreak_sim(R0, dispersion, SI_meanlog, SI_sdlog, 50); test

#______________________________________________________________________
# Set up N outbreak simulations that exceed the last day detection!  

# N <- 10 # 10 successful outbreaks
# introductions <- 1 # for tracking introductions
# outbreak_size <- c() # for tracking successful outbreak sizes
# detectionD <- 100
# 
# # Plot outbreaks that persist too
# par(mfrow = c(3,4))
# 
# while(length(outbreak_size) < N) { # simulate outbreaks until the correct number are reached!
#   outbreak <- outbreak_sim(R0, dispersion, SI_meanlog, SI_sdlog, detectionD) # simulate
#   outbreak_detected <- length(which(outbreak$transD >= detectionD)) # find out if detected
#   print(outbreak_detected) # print!
#   
#   if(outbreak_detected > 0) { # if outbreak persistent record outbrea size and plot outbreak time series
#     outbreak_size <- c(outbreak_size, length(which(outbreak$infD < detectionD)))
#     hist(outbreak$infD, breaks = seq(0, detectionD + 40, 30), ylim = c(0,100), xlim = c(0, detectionD))
#   }
#   introductions <- introductions + 1 # increment introductions!
# }

#______________________________________________________________________
# Function to simulate N outbreaks until a given detection date!  
undetected <- function(N, R0, dispersion, SI_meanlog, SI_sdlog, detectionD) {
  introductions <- 1 # for tracking introductions
  outbreak_size <- c() # for tracking successful outbreak sizes
  
  while(length(outbreak_size) < N) { # simulate outbreaks until the correct number are reached!
    outbreak <- outbreak_sim(R0, dispersion, SI_meanlog, SI_sdlog, detectionD) # simulate
    outbreak_detected <- length(which(outbreak$transD >= detectionD)) # find out if detected
    print(outbreak_detected) # print!
    
    if(outbreak_detected > 0) { # if outbreak persistent record outbrea size and plot outbreak time series
      outbreak_size <- c(outbreak_size, length(which(outbreak$infD < detectionD)))
    }
    introductions <- introductions + 1 # increment introductions!
  }
  outbreak_info <- list(introductions = introductions, 
                        outbreak_size = outbreak_size-1) # remove the detected case!
  return(outbreak_info)
}

# set.seed(14)
# undetected(10, R0, dispersion, SI_meanlog, SI_sdlog, 100)
  
# # cluster 1 (inc. lineages 1 and 5) ------------------------------------------
# timeDiff1_5 = as.numeric(difftime(as.Date("2022-09-30"), #detection date of first case in lineage
#                                   as.Date("2021-07-28"), #tMRCA of lineage
#                                   units = "days"))
# undetect1_5 <- undetected(1000, R0, dispersion, SI_meanlog, SI_sdlog, timeDiff1_5)
# undetect1_5$introductions # 4676
# quantile(undetect1_5$outbreak_size, c(0.025, 0.5, 0.975)) # 31, 493, 3616
# mean(undetect1_5$outbreak_size) # 830
# 
# # lineage 3 -----------------------------------------------
# timeDiff3 = as.numeric(difftime(as.Date("2022-11-14"), as.Date("2022-06-26"),
#                                 units = "days"))
# undetect3 <- undetected(1000, R0, dispersion, SI_meanlog, SI_sdlog, timeDiff3)
# undetect3$introductions
# quantile(undetect3$outbreak_size, c(0.025, 0.5, 0.975)) #34, 3007.975
# mean(undetect3$outbreak_size) #716.124
# 
# 
# ##if cluster 1 is actually two separate introductions:--------------------
# #lineage 1 alone -----------------------------------------------
# timeDiff1 = as.numeric(difftime(as.Date("2022-09-30"), #detection date of first case in lineage
#                                 as.Date("2022-09-24"), #tMRCA of lineage
#                                 units = "days"))
# undetect1 <- undetected(1000, R0, dispersion, SI_meanlog, SI_sdlog, timeDiff1)
# undetect1$introductions
# quantile(undetect1$outbreak_size, c(0.025, 0.5, 0.975)) #34, 3007.975
# mean(undetect1$outbreak_size) #716.124
# 
# #lineage 5 alone -----------------------------------------------
# timeDiff5 = as.numeric(difftime(as.Date("2022-12-25"), as.Date("2022-08-26"),
#                                 units = "days"))
# undetect5 <- undetected(1000, R0, dispersion, SI_meanlog, SI_sdlog, timeDiff5)
# undetect5$introductions
# quantile(undetect5$outbreak_size, c(0.025, 0.5, 0.975)) #34, 3007.975
# mean(undetect5$outbreak_size) #716.124
