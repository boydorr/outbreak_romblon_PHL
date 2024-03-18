#This script simulates an outbreak from when the first case of each lineage was 
#introduced to when the lineage was detected to give an estimate of how many cases were 
#circulating on Romblon in this time. 

set.seed(743)

#R0 parameter values from Townsend et al., 2013
R0 = 1.20
dispersion = 1.33
#SI parameter values from Mancy et al., 2022
SI_meanlog <- 2.85
SI_sdlog <-  0.966

#cluster 1 (inc. lineages 1 and 5) ------------------------------------------
timeDiff1_5 = as.numeric(difftime(as.Date("2022-09-30"), #detection date of first case in lineage
                                as.Date("2021-07-28"), #tMRCA of lineage
                                units = "days"))

prevCases1_5 = c()
while(length(prevCases1_5) < 1000){
  outbreak = data.frame(caseID = 1,
                        parentID = 0,
                        infD = 0,
                        transD = 0)
  
  last_gen = c(1)
  while(all(outbreak$transD[outbreak$caseID %in% last_gen] > timeDiff1_5) == F){ 
    for(i in last_gen){
      offspring = rnbinom(1, mu = R0, size = dispersion)
      if(offspring > 0){
        for(j in 1:offspring){
          secondary = data.frame(caseID = nrow(outbreak) + 1,
                                 parentID = i,
                                 infD = outbreak$transD[outbreak$caseID == i],
                                 transD = outbreak$transD[outbreak$caseID == i]+as.numeric(rlnorm(1, meanlog = SI_meanlog, sdlog = SI_sdlog)))
          
          outbreak = rbind(outbreak, secondary)
        }
      }
    }
    last_gen = outbreak$caseID[outbreak$caseID > max(last_gen)]
    
  }
  if(length(last_gen) > 0){
    prevCases1_5 = append(prevCases1_5, nrow(outbreak[outbreak$transD < timeDiff1_5,]))
  }
}

mean(prevCases1_5) #716.124
quantile(prevCases1_5, c(0.025, 0.975)) #34, 3007.975
median(prevCases1_5) #431

#lineage 3 -----------------------------------------------

timeDiff3 = as.numeric(difftime(as.Date("2022-11-14"), as.Date("2022-06-26"),
                                units = "days"))

prevCases3 = c()
while(length(prevCases3) < 1000){
  outbreak = data.frame(caseID = 1,
                        parentID = 0,
                        infD = 0,
                        transD = 0)
  
  last_gen = c(1)
  while(all(outbreak$transD[outbreak$caseID %in% last_gen] > timeDiff3) == F){
    for(i in last_gen){
      offspring = rnbinom(1, mu = R0, size = dispersion)
      if(offspring > 0){
        for(j in 1:offspring){
          secondary = data.frame(caseID = nrow(outbreak) + 1,
                                 parentID = i,
                                 infD = outbreak$transD[outbreak$caseID == i],
                                 transD = outbreak$transD[outbreak$caseID == i]+as.numeric(rlnorm(1, meanlog = SI_meanlog, sdlog = SI_sdlog)))
          
          outbreak = rbind(outbreak, secondary)
        }
      }
    }
    last_gen = outbreak$caseID[outbreak$caseID > max(last_gen)]
    
  }
  if(length(last_gen) > 0){
    prevCases3 = append(prevCases3, nrow(outbreak[outbreak$transD < timeDiff3,]))
  }
}

mean(prevCases3) #40.8
quantile(prevCases3, c(0.025, 0.975)) #2, 146.025
median(prevCases3) #29

##if cluster 1 is actually two separate introductions:--------------------
#lineage 1 alone -----------------------------------------------

timeDiff1 = as.numeric(difftime(as.Date("2022-09-30"), #detection date of first case in lineage
                                as.Date("2022-09-24"), #tMRCA of lineage
                                units = "days"))

prevCases1 = c()
while(length(prevCases1) < 1000){
  outbreak = data.frame(caseID = 1,
                        parentID = 0,
                        infD = 0,
                        transD = 0)
  
  last_gen = c(1)
  #make sure all branches have reached the observation date or gone extinct before ending the sim!
  while(all(outbreak$transD[outbreak$caseID %in% last_gen] > timeDiff1) == F){ 
    for(i in last_gen){
      offspring = rnbinom(1, mu = R0, size = dispersion)
      if(offspring > 0){
        for(j in 1:offspring){
          secondary = data.frame(caseID = nrow(outbreak) + 1,
                                 parentID = i,
                                 infD = outbreak$transD[outbreak$caseID == i],
                                 transD = outbreak$transD[outbreak$caseID == i]+as.numeric(rlnorm(1, meanlog = SI_meanlog, sdlog = SI_sdlog)))
          
          outbreak = rbind(outbreak, secondary)
        }
      }
    }
    last_gen = outbreak$caseID[outbreak$caseID > max(last_gen)]
    
  }
  if(length(last_gen) > 0){
    prevCases1 = append(prevCases1, nrow(outbreak[outbreak$transD < timeDiff1,]))
  }
}

mean(prevCases1) #1.2
quantile(prevCases1, c(0.025, 0.975)) #1, 3
median(prevCases1) #1

#lineage 5 alone -----------------------------------------------

timeDiff5 = as.numeric(difftime(as.Date("2022-12-25"), as.Date("2022-08-26"),
                                units = "days"))

prevCases5 = c()
while(length(prevCases5) < 1000){
  outbreak = data.frame(caseID = 1,
                        parentID = 0,
                        infD = 0,
                        transD = 0)
  
  last_gen = c(1)
  while(all(outbreak$transD[outbreak$caseID %in% last_gen] > timeDiff5) == F){
    for(i in last_gen){
      offspring = rnbinom(1, mu = R0, size = dispersion)
      if(offspring > 0){
        for(j in 1:offspring){
          secondary = data.frame(caseID = nrow(outbreak) + 1,
                                 parentID = i,
                                 infD = outbreak$transD[outbreak$caseID == i],
                                 transD = outbreak$transD[outbreak$caseID == i]+as.numeric(rlnorm(1, meanlog = SI_meanlog, sdlog = SI_sdlog)))
          
          outbreak = rbind(outbreak, secondary)
        }
      }
    }
    last_gen = outbreak$caseID[outbreak$caseID > max(last_gen)]
    
  }
  if(length(last_gen) > 0){
    prevCases5 = append(prevCases5, nrow(outbreak[outbreak$transD < timeDiff5,]))
  }
}

mean(prevCases5) #30.5
quantile(prevCases5, c(0.025, 0.975)) #1, 114.1
median(prevCases5) #20.5

# CHECKS
source("simulate_helper_fun.R")

pop <- 164012
gr <- 1.023
yrs <- 2015:2024
pop2015_2024 <- pop * gr^(0:length(yrs))
hdr <- 6
dogs2015_2024 <- pop2015_2024/hdr
dogs <- dogs2015_2024[which(yrs == 2021)]
dog_inc <- 0.01
max_cases <- dogs * (timeDiff1_5/ 365) * dog_inc


# cluster 1 (inc. lineages 1 and 5) ------------------------------------------
undetect1_5 <- undetected(1000, R0, dispersion, SI_meanlog, SI_sdlog, timeDiff1_5)
undetect1_5$introductions # 4676
plausible <- which(undetect1_5$outbreak_size < max_cases)
quantile(undetect1_5$outbreak_size[plausible], c(0.025, 0.5, 0.975)) # 31, 493, 3616
mean(undetect1_5$outbreak_size[plausible]) # 830

# lineage 3 -----------------------------------------------
undetect3 <- undetected(1000, R0, dispersion, SI_meanlog, SI_sdlog, timeDiff3)
undetect3$introductions
quantile(undetect3$outbreak_size, c(0.025, 0.5, 0.975)) #34, 3007.975
mean(undetect3$outbreak_size) #716.124


##if cluster 1 is actually two separate introductions:--------------------
#lineage 1 alone -----------------------------------------------
undetect1 <- undetected(1000, R0, dispersion, SI_meanlog, SI_sdlog, timeDiff1)
undetect1$introductions
quantile(undetect1$outbreak_size, c(0.025, 0.5, 0.975)) #34, 3007.975
mean(undetect1$outbreak_size) #716.124

#lineage 5 alone -----------------------------------------------
undetect5 <- undetected(1000, R0, dispersion, SI_meanlog, SI_sdlog, timeDiff5)
undetect5$introductions
quantile(undetect5$outbreak_size, c(0.025, 0.5, 0.975)) #34, 3007.975
mean(undetect5$outbreak_size) #716.124

