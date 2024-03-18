library(tidyverse)
library(sf)
library(ggthemes)
library(here)
library(ggspatial)

## Code used to create fig1c highlighting the region around Romblon and the ports

## Reading in mapping data
province <- st_read(dsn = here("data", "gis"), layer = "PHLsmallTEST_fixed")
province_romblon <- province %>% filter(NAME_1 == "Romblon")

## GIS location of Manila, RADDL 4b and RITM
manila <- list("Longitude" =  14.5995, "Latitude" =  120.9842)
raddl <- list("Longitude" = 13.273705324860375, "Latitude" =  121.23846544455016)
ritm <- list("Longitude" = 14.410539319358085, "Latitude" =  121.03819309623576)
  
ports = read_csv("data/clean/ports.csv")

ports_romblon = ports %>% filter(Name != "airport")
airport = ports %>% filter(Name == "airport")


map = ggplot() + 
  geom_sf(data = province,
          fill = "#DDE0E3",
          linewidth = 0.5) +
  geom_sf(data = province_romblon,
          linewidth = 0.5,
          col = "black",
          fill = "white") +
  geom_point(aes(y = airport$Longitude,
                 x = airport$Latitude),
             col = "goldenrod1"
  ) +
  geom_point(data = ports_romblon, 
             aes(y = Longitude,
                 x = Latitude),
                 colour = "dodgerblue3"
  ) +
  geom_point(aes(y = manila$Longitude,
                 x = manila$Latitude),
             fill = "black",
             shape = 17,
             size = 2
  ) +
  
  geom_point(aes(y = raddl$Longitude,
                 x = raddl$Latitude),
             col = "forestgreen",
             shape = 18,
             size = 2.5
  ) +
  
  geom_point(aes(y = ritm$Longitude,
                 x = ritm$Latitude),
             col = "forestgreen",
             shape = 18,
             size = 2.5
  ) +
  # geom_sf(data = municipality_romblon, linewidth = 0.2, alpha = 0) +
  coord_sf(xlim = c(120.5,123.4), ylim = c(11,14.8), expand = FALSE) +
  theme_map() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.background = element_rect(fill = "aliceblue")) +
  annotation_scale(location = "bl", width_hint = 0.25)

map

png(filename="output/figures/fig1c.png", width = 550, height = 450)
map
dev.off()  

pdf(file="output/figures/fig1c.pdf", width = 5.729, height = 4.6875)
map
dev.off()  
