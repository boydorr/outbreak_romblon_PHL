library(dplyr)
library(tidyr)
library(here)
library(readr)
library(sf)
library(ggplot2)
library(ggthemes)
library(lubridate)

## Additional processing of animal cases for transmission tree analyses

# Kept separately from main data file so as not to create confusion: this script
# mainly generates proxies / assumes values for data properties that in reality
# are unavailable.

dat_outbreak <- read_csv(here("data/clean/dat_animal.csv"))

# Further case data processing ----

dat_outbreak <- dat_outbreak %>%
  filter(province == "Romblon") %>%
  mutate(
  # Add numeric id
         id = as.numeric(substring(case_no, 3)),
  # Exclude case as progenitor if species other than dog (currently all are dogs)
         exclude_progen = if_else(biting_animal %in% "Dog", FALSE, TRUE),
  # Clean up coding for is_owned_biting_animal, set to FALSE if NA/missing
         owned = if_else(is_owned_biting_animal == "Yes", TRUE, FALSE, 
                         missing = FALSE),
  # Add a proxy date for symptoms_started:
  # If biting or sampling date is available, use 2 days before that date;
  # else use 7 days before lab form submission or raddl tests.
         date_symptoms_started_proxy = 
    if_else((is.na(date_biting_incident) & is.na(date_sample_collection)), 
            date - 7, date - 2),
  # Add (assumed) uncertainty around dates:
  # If date is based on biting or sampling date, use 2 days, else use 7 days.
         days_uncertain = 
    if_else((is.na(date_biting_incident) & is.na(date_sample_collection)), 7, 2),
  )

# Filter by date: paper is Sep 2022-Sep 2023 ----

dat_outbreak <- dat_outbreak %>%
  filter(date < ymd("2023-10-01"))

# Add and convert location coordinates: centroids ----

PHL_centroids <- 
  read_csv(here("data/gis/PHL_centroids.csv")) %>% 
  filter(grepl("Romblon",Loc_ID)) %>% 
  st_as_sf(coords=c("Longitude","Latitude"), crs = 4326) %>%
  st_transform(crs = 32651) %>%
  mutate(utm_easting_centroid = do.call(rbind, geometry)[,1],
         utm_northing_centroid = do.call(rbind, geometry)[,2]) %>%
  dplyr::select(Loc_ID, utm_easting_centroid, utm_northing_centroid)

dat_outbreak <-
  dat_outbreak %>%
  mutate(locality = paste(province,municipality,barangay,sep="-")) %>%
  left_join(PHL_centroids,
            by = join_by(locality == Loc_ID))
dat_outbreak[c('geometry')]<- NULL

# Filter columns used in analysis ----

dat_outbreak <-
  dat_outbreak %>%
  dplyr::select(case_no,
                samples_sequenced_seq_id,
                id,
                date_symptoms_started_proxy,
                days_uncertain,
                exclude_progen,
                owned,
                locality,
                utm_easting_centroid,
                utm_northing_centroid)

# Add genetic data ----

lineages <- read_csv(here("output/phylogenetic_distances/0.0004clust.csv")) %>%
  `colnames<-`(c("sample_id", "lineage"))

dat_outbreak %>%
  left_join(lineages, by = join_by(samples_sequenced_seq_id == sample_id)) %>%
  tidyr::replace_na(list(lineage = 0)) -> dat_outbreak

# Give dog R45, traced as biter of the sequenced human case H-23-011Sk_12, the lineage of this human case
dat_outbreak[dat_outbreak$id == 45,"lineage"] <- 
  lineages[lineages$sample_id == "H-23-011Sk_12", "lineage"]

write_csv(dat_outbreak, "data/clean/dat_animal_for_Romblon_tm_trees.csv")

