######################################################################################
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  OUTPUT PROCESSING TO CREATE MAPS   ############################
#############  filibert.heim@posteo.de and steffen.oppel@gmail.com      ##############
######################################################################################

##### Montserrat forest bird surveys data preparation ####
##### written in September 2024 by Filibert Heim  #### 

##### 1: load packages and make other preparations ---------- ####

# load packages
rm(list=ls())
library(RODBC)
library(tmap)
library(tidyverse)
library(data.table)
library(dtplyr)
library(dplyr)
library(basemaps)  ### loads basemap in EPSG:3857 and therefore all other sf objects must be st_transform(3857)
library(sf)
library(terra)
library("rnaturalearth")
library("rnaturalearthdata")
library(gridExtra)
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter


# set working directory
getwd()
#setwd('C:/Users/filib/Documents/Praktika/Sempach/Montserrat') # for Filibert
setwd('C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/UKOT/Montserrat/Analysis/Population_status_assessment/AnnualMonitoring/Montserrat') # for Steffen
# connect with database to load data from query and afterwards close connection again 
db <- odbcConnectAccess2007('data/Montserrat_Birds_2024.accdb') # change name of the db to the actual one 
tblLoc <- sqlFetch(db, 'tblLocation') # coordinates for the points
odbcClose(db)
load("data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData")



#### 2: COLLATE DATA FROM INDIVIDUAL MODEL OUTPUT FILES    -----------------------

#setwd("C:\\STEFFEN\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat\\output")
setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat\\output")

#setwd("C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat\\Montserrat\\output")
nsites<-length(unique(siteCov$Point))
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")
allout<-list.files(pattern="trend_estimates2024_withN.csv")
annestimates<-tibble()
trendout<-tibble()
for (f in allout){
  x<-fread(f)
  annestimates<-annestimates %>% bind_rows(x %>% filter(substr(parameter,1,1)=="N"))
}


### COMBINE INFO WITH COORDINATES ###

points<- tblLoc %>%
  filter(Point %in% unique(siteCov$Point)) %>%
  #filter(FOREST=="CH") %>%
  select(Point,Eastings,Northings,Altitude,Habitat_Code) %>%
  arrange(Point) %>%
  mutate(order=seq_along(Point))


out<-annestimates %>% separate_wider_delim(.,cols=parameter,delim=", ", names=c("Point","Year"))
out$Point<-str_replace_all(out$Point,pattern="[^[:alnum:]]", replacement="")
out$Year<-str_replace_all(out$Year,pattern="[^[:alnum:]]", replacement="")
out$order<-str_replace_all(out$Point,pattern="N", replacement="")

mapdata<-out %>%
  mutate(Year=as.numeric(Year)+2010) %>%
  mutate(order=as.numeric(order)) %>%
  select(-Point) %>%
  left_join(points, by="order") %>%
  select(species,Year,Point, Eastings,Northings,Altitude,Habitat_Code,mean, median, lcl,ucl)

fwrite(mapdata,"Annual_estimates2024_mapdata.csv")







############################################################################################################
####   PRODUCE FIGURES FOR POPULATION TREND AND DETECTION PROBABILITY          #############################
############################################################################################################

mapdata_sf<-mapdata %>%
  st_as_sf(coords = c("Eastings", "Northings"), crs=2004) %>%
  st_transform(4326)                        ## Montserrat is EPSG 4604 or 2004


# create bounding box
bbox <- st_sfc(st_point(c(-62.25,16.68)), st_point(c(-62.1, 17.3)), crs = 4326) %>% st_bbox()

# extract tiles - check the provider argument, there are many options
basemap <- maptiles::get_tiles(x = bbox, 
                               zoom = 12,
                               crop = TRUE, provider = "OpenTopoMap")



tmap_mode("view")
tm_shape(basemap)+
  tm_rgb()+
  tm_shape(mapdata_sf)  +
  tm_symbols(col = "red", size = 0.2)



