


library(sf)
library(here)
library(ggthemes)
library(dplyr)

## Reading in mapping data
province <- st_read(dsn = here("data/gis"), layer = "PHLsmallTEST_fixed")
municipality <- st_read(dsn = here("data/gis"), layer = "PHL_municipality")
municipality_loc <- read.csv(here("data/gis", "PHL_municipality_centroids.csv"))
#province_cropped <- province %>% filter(NAME_1 == "Romblon, Batangas, Bulacan, Laguna, Metropolitan Manila, Quezon")

## Reading in metadata
dat_cluster <- read.csv("output/phylogenetic/romblon_subtrees/cluster1.csv")
dat_cluster <- dat_cluster %>% filter(isTip == "TRUE" & Province == "Romblon")


## process GIS data
#village_loc <- village_loc[grep("Romblon-", village_loc$Loc_ID),]
dat_cluster$loc_id <- paste0(dat_cluster$Province, "-", dat_cluster$Municipality)
dat_cluster$Lat <- municipality_loc$Latitude[match(dat_cluster$loc_id, municipality_loc$Loc_ID)]
dat_cluster$Lon <- municipality_loc$Longitude[match(dat_cluster$loc_id, municipality_loc$Loc_ID)]

dat_cluster$provigeometry <- province$geometry[match(dat_cluster$Province,province$NAME_1)]

cluster_grp <- dat_cluster %>%
  group_by(cluster) %>%
  summarise(n=n()) 
cluster_grp

dat_cluster <- dat_cluster %>% 
  mutate(cluster = factor(cluster, 
                          levels = c("1")))
romblon <- province[province$NAME_1 == "Romblon",]
cluster1.map <- ggplot() + 
  geom_sf(data = romblon, #aes(fill = Province, geometry = provigeometry), 
          alpha = 0,
          linewidth = 0.4,
          colour = "#ADB5BD")+
  geom_sf(data=romblon , fill="darkred", alpha=0.5)+
  coord_sf(xlim = c(121.8,122.4), ylim = c(12,12.75), expand = FALSE) +
  # geom_point(data = dat_cluster, aes(x = Lon,
  #                                    y = Lat,
  #                                    fill = factor(lineage_chain)),size=3, pch=21) +
  geom_jitter(data = dat_cluster, aes(x = Lon,
                                      y = Lat, fill = factor(lineage_chain)),
              width = 0.02, height = 0.02,size=3, pch=21) +
  scale_fill_manual(values=c("#008CF9", "#006E00", "#00BBAD"),
                     labels = c('1','4','5'),na.translate=FALSE, 
                     name = "Lineage")+
  guides(fill="none")+
  theme_map() #this is to differently format your map, feel free to remove :) 

cluster1.map


