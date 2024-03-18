library(gsheet)
library(tidyverse)
library(janitor)
library(zoo)
source("rabid_creds.r")

dat_outbreak <- gsheet2tbl(dat_outbreak_link)
dat_human <- gsheet2tbl(dat_human_link)

write_csv(dat_outbreak, "data/raw/dat_animal.csv")
write_csv(dat_human, "data/raw/dat_human.csv")

dat_human <- dat_human %>%
  clean_names(.,"snake")  ## Using janitor package to clean names and change to lower case

dat_outbreak <- dat_outbreak %>% row_to_names(row_number = 1) %>%
  clean_names(.,"snake")  ## Using janitor package to clean names and change to lower case

table(dat_outbreak$case_status, dat_outbreak$province)

dat_outbreak <- dat_outbreak %>%
  select(case_no,
         samples_sequenced_seq_id,
         barangay,
         municipality,
         province,
         case_status,
         date_sample_collection,
         date_biting_incident,
         date_lab_submission_form,
         date_animal_died,
         biting_animal,
         is_owned_biting_animal,
         date_tests_raddl,
         date_tests_ritm,
         date_tests_pvo,
         date_tests_mao,
         date_tests_speedier,
         tests_speedier,
         tests_raddl_rdt,
         tests_pvo,
         tests_mao,
         exposures)

dat_outbreak <- dat_outbreak %>%
  mutate(across(contains('date'), ymd))

dat_outbreak <- dat_outbreak %>%
  mutate(date = if_else(is.na(date_biting_incident),date_sample_collection,date_biting_incident)) %>%
  mutate(date = if_else(is.na(date),date_lab_submission_form,date)) %>%
  mutate(date = if_else(is.na(date),date_tests_raddl,date)) %>%
  mutate(date = as_date(date)) %>%
  mutate(month_year = as.yearmon(date))

## Only interested in confirmed and probable cases 
dat_outbreak_all <- dat_outbreak %>%
  filter(case_status != "Duplicate") %>% 
  drop_na(date) ## No cases dopped, all have a date


## Only interested in confirmed and probable cases 
dat_outbreak <- dat_outbreak %>%
  filter(case_status == "Confirmed" | case_status == "Probable") %>%
  filter(case_status != "Duplicate") %>% 
  drop_na(date) ## No cases dopped, all have a date

# dat_outbreak <- dat_outbreak %>% 
#   mutate(biting_animal = ifelse(biting_animal == "DOG","Dog",NA))

## create one date column 
## if available use when_bitten
## else use date_stated symptoms

dat_human <- dat_human %>%
  mutate(when_bitten = ymd(when_bitten)) %>%
  # mutate(date_started_symptoms = parse_date_time2(date_started_symptoms, 
  #                                                 orders = c("mdY","mdy", "Ydm","Ymd"))) %>%
  mutate(date_started_symptoms = ymd(date_started_symptoms)) %>%
  mutate(date = if_else(is.na(when_bitten), date_started_symptoms,when_bitten)) 

dat_human <- dat_human %>% 
  mutate(date_died = ymd(date_died))
  # mutate(date_died = parse_date_time2(date_died, 
  #                                     orders = c("mdY","mdy", "Ydm","Ymd")))

dat_human <- dat_human %>%
  mutate(month_year = as.yearmon(date)) %>%
  mutate(year =  lubridate::year(date)) %>% 
  select(date, date_died, when_bitten, date_started_symptoms, province,municipality,barangay,month_year, case_status)

write_csv(dat_human, "data/clean/dat_human.csv")
write_csv(dat_outbreak, "data/clean/dat_animal.csv")
write_csv(dat_outbreak_all, "data/clean/dat_outbreak_all.csv")


