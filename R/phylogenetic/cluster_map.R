#'---------------------------------------------------------
#'title: Cluster map 
#'author: Criselda T. Bautista
#'date: 06/03/2024
#'---------------------------------------------------------
# Scripts to consolidate variables and fix names for cleaning data
rm(list=ls())  # This clears anything currently in the environment
setwd("~/Documents/GitHub/Romblon_Outbreak")

#############################################
#  INSTALL PACKAGES                         #
############################################# 
library(tidyverse)
library(lubridate)
library(incidence)
library(mapview)
library(sp) 
library(sf) 
library(rgdal)
library(ggplot2)
library(rworldmap)
library(cleangeo)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(mapproj)
library(sf)
library(here)

#############################################
#              LOAD FUNCTIONS               #
#############################################
source("romblon_clean/scripts/R/shp_to_df.R")
source("romblon_clean/scripts/R/create_discrete_scale_for_GLUE.R")
# These scripts were created by Rachel Steenson 31/01/2020

#############################################
#  IMPORT DATA                              #
#############################################
# Load shapefiles and create df (for plots)
#province <- readOGR("data", "PHLsmallTEST_fixed", stringsAsFactors=FALSE)
province <- read_sf("romblon_clean/gis_data/PHL/PHLsmallTEST_fixed.shp")
prov <- read.csv("romblon_clean/gis_data/PHL/PHL_province_centroids.csv")
muni <- read.csv("romblon_clean/gis_data/PHL/PHL_municipality_centroids.csv")
municipality <- read_sf("romblon_clean/gis_data/PHL/PHL_municipality.shp")
#muni <- readOGR("data/PHL_municipality.shp")
#case_muni <- unique(mdo_meta$municipality)

## Reading in metadata
cluster <- read.csv("romblon_clean/processed_data/romblon_subtrees/cluster1_3_cases.csv")
cluster <- cluster %>% filter(isTip == "TRUE")

plot(province)

#############################################
#           PROCESS THE DATA                #
#############################################

cluster$loc_id <- paste0(cluster$Province, "-", cluster$Municipality)
cluster$Lat <- muni$Latitude[match(cluster$loc_id, muni$Loc_ID)]
cluster$Lon <- muni$Longitude[match(cluster$loc_id, muni$Loc_ID)]

cluster$provgeometry <- province$geometry[match(cluster$Province,province$NAME_1)]

cluster$munigeometry <- municipality$geometry[match(cluster$Municipality,municipality$NAME_2)]

# Make a quick table of each cluster and the number of samples associated with it
cluster_grp <- cluster %>%
  group_by(cluster) %>%
  summarise(n=n()) 
cluster_grp # Here you can see all the messy province names still!

cluster <- cluster %>% 
  mutate(cluster = factor(cluster, 
                            levels = c("1","3")))

basemap <- ggplot()+
  geom_sf(data=province, fill="white", color="black") # TO PLOT THIS, THE SHAPEFILE NEEDED IMPORTING USING read_sf (I changed on line 38)
basemap

png(paste0("figures/basemap.png"), width=25,height=12,units="cm",res=300) #png makes an image 
basemap 
dev.off() ; basemap # done making an image 

cluster_prov <- basemap+
  geom_point(data = cluster, aes(x = Lon, y = Lat, color  = cluster ))
cluster_prov

basemap_municipality <- ggplot()+
  geom_sf(data=municipality, fill="white", color="black") # TO PLOT THIS, THE SHAPEFILE NEEDED IMPORTING USING read_sf (I changed on line 38)
basemap_municipality

png(paste0("figures/basemap_muni.png"), width=25,height=12,units="cm",res=300) #png makes an image 
basemap_municipality 
dev.off() ; basemap_municipality # done making an image 

cluster_muni <- basemap_municipality+
  geom_sf(data = cluster, aes(fill = cluster, geometry = munigeometry)) + 
  scale_y_continuous(breaks = 34:36)
cluster_muni

png(paste0("figures/cluster_muni.png"), width=25,height=12,units="cm",res=300) #png makes an image 
cluster_muni 
dev.off() ; cluster_muni # done making an image 

