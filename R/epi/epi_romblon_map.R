library(sf)
library(here)
library(tidyverse)
library(janitor)
library(zoo)
library(ggthemes)
library(ggspatial)

## Code used to create epi map of cases (fig 2c)
dat_outbreak <- read_csv("data/clean/dat_animal.csv")
dat_human <- read_csv("data/clean/dat_human.csv")

start_date <- as.Date("2022-09-30")
end_date <- as.Date("2023-09-30")


## Only interested in confirmed and probable cases in Romblon from Sep 2022 to Sep 2023
dat_outbreak <- dat_outbreak %>% filter(province == "Romblon") %>% filter(date >= start_date & date <= end_date)
dat_human <- dat_human %>% filter(province == "Romblon") %>% filter(date >= start_date & date <= end_date)

## Reading in mapping data
province <- st_read(dsn = here("data", "gis"), layer = "PHLsmallTEST_fixed")
municipality <- st_read(dsn = here("data", "gis"), layer = "PHL_municipality")
province_cropped <- province %>% filter(NAME_1 != "Romblon")
village <- st_read(dsn = here("data", "gis"), layer = "PHL_village")
village <- village %>% filter(NAME_1 == "Romblon")
village_loc <- read.csv(here("data", "gis", "PHL_village_centroids.csv"))
village_loc <- village_loc[grep("Romblon-", village_loc$Loc_ID),]



dat_human <- dat_human %>%
  mutate(date_bitten = ymd(when_bitten),
         date = date_died) ## plotting the date died rather than date bitten

### Combine animal and human case data within one data frame for epicurve
## Creating a "Human Death" case status so that the human case can be included on the same epi curve as the animal cases
dat_human <- dat_human %>%
  select(date, province,municipality,barangay,month_year) %>% 
  mutate(case_status = "Human Death") 

dat_outbreak <- dat_outbreak %>%
  select(barangay,municipality,province,case_status,date,month_year)

## Combining the data frames
dat_outbreak <- rbind(dat_outbreak,dat_human) %>%
  mutate(case_status_fct = 
           factor(case_status, 
                  levels = c("Human Death","Confirmed","Probable")))

levels(dat_outbreak$case_status_fct) <- 
  c("Human Death", 'Confirmed (Animal)',  "Probable (Animal)")


## Combining the data frames
dat_outbreak <- dat_outbreak %>% 
  mutate(case_type = ifelse(case_status == "Human Death", "Human", "Animal")) %>% 
  mutate(case_type = factor(case_type, 
                            levels = c("Human","Animal")))

levels(dat_outbreak$case_type) <-  c("Human ", "Animal")

## process GIS data
village_loc <- village_loc[grep("Romblon-", village_loc$Loc_ID),]
dat_outbreak$loc_id <- paste0(dat_outbreak$province, "-", dat_outbreak$municipality, "-", dat_outbreak$barangay)
dat_outbreak$Lat <- village_loc$Latitude[match(dat_outbreak$loc_id, village_loc$Loc_ID)]
dat_outbreak$Lon <- village_loc$Longitude[match(dat_outbreak$loc_id, village_loc$Loc_ID)]

dat_outbreak <- dat_outbreak %>% 
  left_join(village, by = c("loc_id" = "Loc_ID"))


dat_outbreak <- dat_outbreak %>%
  mutate(recent_case = ifelse(date >= (Sys.Date() - 30),"Yes","No"))

dat_outbreak <- dat_outbreak %>%
  mutate(case_type = 
           ifelse(case_status == "Confirmed" | case_status == "Probable", "Animal","Human"))

### subset for Romblon
municipality_df_romblon <- filter(municipality, NAME_1 == "Romblon")

# partially transparent points by setting `alpha = 0.5`

# fillcolors <- c("#014D93", "#c1121f", "#ff9f1c")
fillcolors <- c("#ae2012", "#001219")


dat_outbreak_grp <- dat_outbreak %>% 
  group_by(barangay, case_type,loc_id) %>% 
  summarise(n = n(),
            Lat = first(Lat),
            Lon = first(Lon))

dat_outbreak_grp <- dat_outbreak_grp %>% 
  mutate(case_type = factor(case_type, 
                            levels = c("Human","Animal"))) 

dat_outbreak_grp <- dat_outbreak_grp %>% 
  mutate(n = ifelse(case_type == "Human", 1.5, n))

dat_outbreak_grp <- dat_outbreak_grp %>% 
  mutate(case_type = ifelse(case_type == "Human", "Human death", "Animal case"))

dat_outbreak_grp <- dat_outbreak_grp %>% 
  mutate(Lat = ifelse(case_type == "Human death", Lat + 0.001, Lat)) %>% 
  mutate(Lon = ifelse(case_type == "Human death", Lon + 0.002, Lon))


map_plot <- ggplot() + 
  geom_sf(data = village, 
          alpha = 0.1,
          colour = "grey80",
          linewidth = 0.3) +
  geom_sf(data = municipality_df_romblon,
          col = "grey50",
          alpha = 0,
          linewidth = 0.6) +
  coord_sf(xlim = c(121.8,122.4), ylim = c(12,12.75), expand = FALSE) +
  geom_point(data = dat_outbreak_grp, aes(x = Lon,
                                          y = Lat,
                                          col = case_type,
                                          shape = case_type,
                                          alpha = case_type,
                                          size = n))+
  scale_color_manual(values=c("#ae2012", "#001219"), name = "") +
  scale_shape_manual(values= c(16, 17), name = "") +
  scale_alpha_manual(values = c(0.9,1)) +
  scale_size_continuous(range = c(1.5,5)) +
  theme_map(base_size = 12) + 
  theme(legend.position = c(-0.2, 0.5),
        legend.box="vertical",
        legend.title = element_text(size=12,face = "plain"),
        legend.text = element_text(size=11,face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 4, shape = c(16,17)), order = 1),
         size = guide_legend(title = "Number of animal cases", override.aes = list(col = "#ae2012", order = 2)),
         alpha = "none", 
         shape = "none") +
  annotation_scale(location = "br", width_hint = 0.25)
map_plot

png(filename="output/figures/fig2c.png", width = 550, height = 450)
map_plot
dev.off()


pdf(file = "output/figures/fig2c.pdf",width = 8, height = 7)
map_plot
dev.off()

