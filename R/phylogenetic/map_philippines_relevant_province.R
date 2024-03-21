
library(sf)
library(here)
library(ggthemes)
library(ggspatial)

## Reading in mapping data
province <- st_read(dsn = here("gis_data", "PHL"), layer = "PHLsmallTEST_fixed")
municipality <- st_read(dsn = here("gis_data", "PHL"), layer = "PHL_municipality")
municipality_loc <- read.csv(here("gis_data", "PHL", "PHL_municipality_centroids.csv"))


# scale_fill_manual(values=c("deeppink","blue4","lightsalmon","darkred"),
#                   labels = sort(unique(cluster3$data$Province)) )
#"Batangas" "Laguna"   "Quezon"   "Romblon" 
batangas <- province[province$NAME_1 == "Batangas",]
laguna <- province[province$NAME_1 == "Laguna",]
quezon <- province[province$NAME_1 == "Quezon",]
romblon <- province[province$NAME_1 == "Romblon",]
bulacan <- province[province$NAME_1 == "Bulacan",]
metro <- province[province$NAME_1 == "Metropolitan Manila",]
romblon <- province[province$NAME_1 == "Romblon",]

philippines <- ggplot() + 
  geom_sf(data = province, #aes(fill = Province, geometry = provigeometry), 
          alpha = 1,
          linewidth = 0.2,
          colour = "#ADB5BD") +
  geom_sf(data=batangas , fill="deeppink",col="black")+
  geom_sf(data=laguna , fill="blue4",col="black")+
  geom_sf(data=quezon , fill="lightsalmon",col="black")+
  geom_sf(data=romblon , fill="darkred",col="black")+
  geom_sf(data=bulacan , fill="purple",col="black")+
  geom_sf(data=metro , fill="pink",col="black")+
coord_sf(xlim = c(120,123), ylim = c(12,16), expand = FALSE) +
  
  theme_map() +#this is to differently format your map, feel free to remove :) 
  theme(panel.background = element_rect(fill = "white"),
        panel.border= element_rect(color = "black", fill=NA),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) +
  annotation_scale(style="ticks");
philippines





