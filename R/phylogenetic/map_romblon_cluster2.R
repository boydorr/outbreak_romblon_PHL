library(sf)
library(here)
library(ggthemes)
library(dplyr)
library(ggspatial)

## Reading in mapping data
province <- st_read(dsn = here("gis_data", "PHL"), layer = "PHLsmallTEST_fixed")
municipality <- st_read(dsn = here("gis_data", "PHL"), layer = "PHL_municipality")
municipality_loc <- read.csv(here("gis_data", "PHL", "PHL_municipality_centroids.csv"))

## Reading in metadata
dat_cluster2 <- read.csv("processed_data/romblon_subtrees/cluster2_cases.csv")
# remove the label filter if want to plot related tc 3 cases
dat_cluster2 <- dat_cluster2 %>% filter(isTip == "TRUE" & Province == "Romblon" & label=="H-23-011Sk_12")
dat_cluster2$cluster <- "2"

## process GIS data
dat_cluster2$loc_id <- paste0(dat_cluster2$Province, "-", dat_cluster2$Location)
dat_cluster2$Lat <- municipality_loc$Latitude[match(dat_cluster2$loc_id, municipality_loc$Loc_ID)]
dat_cluster2$Lon <- municipality_loc$Longitude[match(dat_cluster2$loc_id, municipality_loc$Loc_ID)]

cluster_grp2 <- dat_cluster2 %>%
  group_by(cluster) %>%
  summarise(n=n()) 
cluster_grp2

dat_cluster2 <- dat_cluster2 %>% 
  mutate(cluster = factor(cluster, 
                          levels = c("1")))

romblon <- province[province$NAME_1 == "Romblon",]
cluster2.map <- ggplot() + 
  geom_sf(data = romblon, #aes(fill = Province, geometry = provigeometry), 
          alpha = 0,
          linewidth = 0.4,
          colour = "#ADB5BD")+
  scale_fill_manual(values=c("#B80058","#EBAC23"),
                    labels = c('2','3'),na.translate=FALSE, 
                    name = "Lineage")+
  guides(fill="none")+
  geom_sf(data=romblon , fill="darkred", alpha=0.5)+

  geom_jitter(data = dat_cluster2, aes(x = Lon,
                                       y = Lat, fill = as.factor(lineage_chain)),
              width = 0.02, height = 0.02,size=3, pch=21) +
  guides(fill="none")+
  theme_map(); #this is to differently format your map, feel free to remove :) 
cluster2.map


