library(sf)
library(here)
library(ggthemes)
library(ggspatial)
library(dplyr)

## Reading in mapping data
province <- st_read(dsn = here("gis_data", "PHL"), layer = "PHLsmallTEST_fixed")
municipality <- st_read(dsn = here("gis_data", "PHL"), layer = "PHL_municipality")
municipality_loc <- read.csv(here("gis_data", "PHL", "PHL_municipality_centroids.csv"))

## Reading in metadata
dat_cluster3 <- read.csv("processed_data/romblon_subtrees/cluster3.csv")
dat_cluster3 <- dat_cluster3 %>% filter(isTip == "TRUE" & Province == "Romblon")

## process GIS data
dat_cluster3$loc_id <- paste0(dat_cluster3$Province, "-", dat_cluster3$Municipality)
dat_cluster3$Lat <- municipality_loc$Latitude[match(dat_cluster3$loc_id, municipality_loc$Loc_ID)]
dat_cluster3$Lon <- municipality_loc$Longitude[match(dat_cluster3$loc_id, municipality_loc$Loc_ID)]

cluster_grp <- dat_cluster3 %>%
  group_by(cluster) %>%
  summarise(n=n()) 
cluster_grp

dat_cluster3 <- dat_cluster3 %>% 
  mutate(cluster = factor(cluster, 
                          levels = c("1")))

romblon <- province[province$NAME_1 == "Romblon",]
cluster3.map <- ggplot() + 
  geom_sf(data = romblon, #aes(fill = Province, geometry = provigeometry), 
          alpha = 0,
          linewidth = 0.4,
          colour = "#ADB5BD")+
  geom_sf(data=romblon , fill="darkred", alpha=0.5)+
  # geom_point(data = dat_cluster, aes(x = Lon,
  #                                    y = Lat,
  #                                    fill = factor(lineage_chain)),size=3, pch=21) +
  geom_jitter(data = dat_cluster3, aes(x = Lon,
                                      y = Lat, fill = factor(lineage_chain)),
              width = 0.02, height = 0.02,size=3, pch=21) +
  scale_fill_manual(values=c("#EBAC23"),
                    labels = c('3'),na.translate=FALSE, 
                    name = "Lineage")+
  guides(fill="none")+
  theme_map(); #this is to differently format your map, feel free to remove :) 
cluster3.map


