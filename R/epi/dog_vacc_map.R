library(tidyverse)
library(readxl)
library(sf)
library(here)

dat = read_excel("data/raw/MDV_2022_2023.xlsx")

dat_sum = dat %>% 
  group_by(Municipality) %>% 
  summarise(n = sum(Vaccinated_Dogs))


municipality <- st_read(dsn = here("data", "gis"), layer = "PHL_municipality")
municipality_romblon <- municipality %>% filter(NAME_1 == "Romblon")
municipality_df_romblon_cropped <- municipality_romblon %>% 
  filter(NAME_2 != "Corcuera" & NAME_2 != "Banton" & NAME_2 != "Concepcion" & NAME_2 != "San Fernando" & NAME_2 != "Cajidiocan" & NAME_2 != "Magdiwang")

municipality_df_romblon_cropped <- municipality_df_romblon_cropped %>%
  left_join(dat_sum, by = c("NAME_2" = "Municipality")) 
  # mutate(n = ifelse(is.na(n),0,n))

centroids = read_csv("data/gis/PHL_centroids.csv")
centroids = centroids %>% filter(Type == "Municipality")

dat_points = municipality_df_romblon_cropped %>% 
  mutate(Loc_ID = paste0(NAME_1,"-",NAME_2)) %>% 
  left_join(centroids, by = "Loc_ID")


map_vac = ggplot() + 
  geom_sf(data = municipality_df_romblon_cropped,
          alpha = 0.8, 
          colour = "dimgrey",
          linewidth = 0.5) +
  geom_point(data = dat_points,
             aes(x = Longitude, 
                 y = Latitude, 
                 size = n),
             col = "darkblue",
             alpha = 0.5) +
  scale_size(range = c(1,14), breaks = c(500,1000,1500,2000,2500,3000,4000),name = "Dog vaccination") +
  coord_sf(xlim = c(121.8,122.4), ylim = c(12,12.75), expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank()) 
map_vac

png(filename="output/figures/dog_vac_map.png", width = 550, height = 450)
map_vac
dev.off() 
