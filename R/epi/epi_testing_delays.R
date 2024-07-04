library(here)
library(ggthemes)
library(gsheet)
library(tidyverse)
library(janitor)
library(zoo)

dat_outbreak <- read_csv("data/clean/dat_outbreak_all.csv") ## includes negative cases

start_date <- as.Date("2022-09-30")
end_date <- as.Date("2023-09-30")

dat_outbreak <- dat_outbreak %>% filter(province == "Romblon") %>% filter(date >= start_date & date <= end_date)

## Time between bite incidence and RDT
## RDTs have been performed by RADDL and SPEEDIER (MAO and PVO however no date has been provided)
dat_outbreak_delay <- dat_outbreak %>% 
  filter(!is.na(tests_raddl_rdt) | !is.na(tests_speedier)) %>% 
  filter(tests_raddl_rdt != "Unknown")

## create one date for the laboratory date
dat_outbreak_delay <- dat_outbreak_delay %>% 
  mutate(date_lab_confirmation = as.Date(ifelse(is.na(date_tests_raddl), date_tests_speedier, date_tests_raddl))) %>% 
  mutate(date_lab_confirmation = as.Date(ifelse(is.na(date_lab_confirmation), date_sample_collection, date_lab_confirmation)))

## Only look at those where there was an exposure (biting incident)
ani_td <- dat_outbreak_delay %>% filter(exposures >= 1)

## Remove those which do not have a biting incident date or lab confirmation date
ani_td <- ani_td %>% filter(!is.na(date_biting_incident)) 
ani_td <- ani_td %>% filter(!is.na(date_lab_confirmation))

## Calculate time difference
ani_td <- ani_td %>% 
  mutate(time_delay = date_lab_confirmation - date_biting_incident) %>% 
  mutate(time_delay = as.numeric(time_delay))

## mean and CI
mean(ani_td$time_delay)
median(ani_td$time_delay)

l_model <- lm(time_delay ~ 1, ani_td)
confint(l_model, level=0.95)


## Time between biting incident and dog death
##  only interested in confirmed and probable cases
ani_td_death <- dat_outbreak %>% filter(exposures >= 1)

ani_td_death <- ani_td_death %>% filter(case_status != "dFAT Negative")

ani_td_death <- ani_td_death %>% filter(!is.na(date_biting_incident)) 

ani_td_death <- ani_td_death %>% 
  mutate(date_died = ifelse(is.na(date_animal_died),date_sample_collection,date_animal_died))

ani_td_death <- ani_td_death %>% filter(!is.na(date_died))

ani_td_death <- ani_td_death %>% 
  mutate(time_delay = date_animal_died - date_biting_incident) %>% 
  mutate(time_delay = as.numeric(time_delay))

mean(ani_td_death$time_delay)
median(ani_td_death$time_delay)

l_model <- lm(time_delay ~ 1, ani_td_death)
confint(l_model, level=0.95)

## Calculate the average number of exposures
dat_outbreak_exp <- dat_outbreak %>% 
  filter(case_status != "dFAT Negative") %>% 
  mutate(exposures = as.numeric(exposures))

mean(dat_outbreak_exp$exposures, na.rm = T)



## Time delay between biting incident and PEP
##  only interested in confirmed and probable cases
ani_td_pep <- dat_outbreak %>% filter(exposures >= 1)


ani_td_pep <- ani_td_pep %>% filter(!is.na(date_biting_incident)) 

ani_td_pep <- ani_td_pep %>% 
  mutate(pep_date = ifelse(is.na(pep_date), pep_date_init, pep_date)) %>% 
  filter(!is.na(pep_date)) %>% 
  mutate(pep_date = as.Date(pep_date))

ani_td_pep <- ani_td_pep %>% 
  mutate(time_delay =  pep_date - date_biting_incident) %>% 
  mutate(time_delay = as.numeric(time_delay))

mean(ani_td_pep$time_delay)
median(ani_td_pep$time_delay)

l_model <- lm(time_delay ~ 1, ani_td_pep)
confint(l_model, level=0.95)

#% of rabies-negative dogs had a history of rabies vaccination.

dat_outbreak_vac = dat_outbreak %>% 
  filter(case_status == "dFAT Negative")

dat_outbreak_vac <- dat_outbreak_vac %>% 
  mutate(vacstatus_biting_animal = ifelse(is.na(vacstatus_biting_animal), "Unknown",vacstatus_biting_animal))

table(dat_outbreak_vac$vacstatus_biting_animal)
nrow(dat_outbreak_vac)

(nrow(dat_outbreak_vac[dat_outbreak_vac$vacstatus_biting_animal == "Vaccinated",])/nrow(dat_outbreak_vac))*100

