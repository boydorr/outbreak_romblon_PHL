library(dplyr)
library(tidyr)
library(here)
library(readr)
library(sf)
library(raster)
library(data.table)

## Preparation of population size data.frame to use for simulating case locations

# Load Romblon island shape file
PHL_shp <- st_read(here("data/gis/PHL_village.shp"))
Romblon_shp <- PHL_shp[PHL_shp$NAME_1 == "Romblon",]

# # WorldPop population estimates at 100m x 100m resolution, for 2020
# # Downloaded from https://hub.worldpop.org/geodata/summary?id=6316
# PHL_pop_unconstrained <- raster(here("data/gis/phl_ppp_2020.tif"))
# # Using top-down unconstrained model, as in the constrained model some Barangays 
# # (including some associated with cases) have only cells with value NA,
# # being "unpopulated". 
# # See https://www.worldpop.org/methods/top_down_constrained_vs_unconstrained/ 
# # for model info. 
# # Crop pop estimates to Romblon to make manageable size for Github repo
# Romblon_pop_unconstrained <-
#   crop(PHL_pop_unconstrained, Romblon_shp, here("data/gis/Romblon_ppp_2020.tif"))
Romblon_pop_unconstrained <- raster(here("data/gis/Romblon_ppp_2020.tif"))

# The above tif file is by WorldPop and is licensed under CC BY 4.0. 
# For more details, see the `LICENSE` file. 
# For the full license text, please refer to: https://creativecommons.org/licenses/by/4.0/


# Extract population estimates for all 100m x 100m cells in Romblon's barangays
# (This counts a cell as being in a barangay if its central point is located 
# within the barangay boundaries)
Romblon_pop_unconstrained %>%
  raster::extract(., Romblon_shp, cellnumbers = TRUE, df = TRUE) %>%
  tibble() %>%
  mutate(loc_ID = Romblon_shp$Loc_ID[ID],
         cell = cell,
         population = Romblon_ppp_2020,
         .keep = "none") -> Romblon_pop_by_cell

# Extract cell coordinates & convert to UTM
Romblon_pop_by_cell %<>%
  xyFromCell(object = Romblon_pop_unconstrained, cell = .$cell, 
             spatial = FALSE) %>%
  as.data.frame() %>%
  st_as_sf(coords=c("x","y"), crs = 4326) %>%
  st_transform(crs = 32651) %>%
  mutate(utm_easting = do.call(rbind, geometry)[,1],
         utm_northing = do.call(rbind, geometry)[,2]) %>%
  as_tibble() %>%
  dplyr::select(-geometry) %>%
  cbind(Romblon_pop_by_cell, .)

fwrite(Romblon_pop_by_cell, here("output/population_sizes/Romblon_pop_by_cell.csv"))
