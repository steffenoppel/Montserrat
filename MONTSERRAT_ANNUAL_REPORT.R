#.....................................................................................
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  CREATION OF ANNUAL REPORT    ##########################################
#############  steffen.oppel@gmail.com      ##########################################
#.....................................................................................

### based on multi-year trend model of Kery et al. 2009
### REQUIRES ANALYIS TO BE COMPLETED IN  "MONTSERRAT_BIRD_MONITORING_Nimble.R"


#.....................................................................................
# 1.  Load required packages       ---------------------
#.....................................................................................

rm(list=ls())
library(tidyverse)
library(data.table)
library(dplyr)
library(dtplyr)
library(lubridate)
library(ggplot2)
library(knitr)
library(rmarkdown)
library(tmap)
library(basemaps)  ### loads basemap in EPSG:3857 and therefore all other sf objects must be st_transform(3857)
library(sf)
library(terra)
library("rnaturalearth")
library("rnaturalearthdata")
library(gridExtra)
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

YEAR<- if_else(month(Sys.time())<6,year(Sys.time())-1,year(Sys.time()))


#.....................................................................................
# 2. Set your working directory to find and load model outputs ---------
#.....................................................................................

#setwd("C:\\STEFFEN\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")
setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")
#setwd("C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat\\Montserrat")



#.....................................................................................
## 2.1 load the pre-prepared modelling INPUT dataset					--------
#.....................................................................................

load(sprintf("data/MONTSERRAT_ANNUAL_DATA_INPUT%s.RData",YEAR))
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")



#.....................................................................................
## 2.2 load the pre-prepared modelling OUTPUT estimates					----------
#.....................................................................................



annestimates<-read.table(sprintf("output/Annual_estimates%s.csv",YEAR), header=T, sep=",")
trendout<-read.table(sprintf("output/Trend_estimates%s.csv",YEAR), header=T, sep=",")
mapdata<-fread(sprintf("output/Annual_estimates%s_mapdata.csv",YEAR))




#.....................................................................................
## 2.3 manipulate the loaded data					----------
#.....................................................................................

## add full species names
annestimates$fullspec<-fullnames[match(annestimates$species, SPECIES)]
trendout$fullspec<-fullnames[match(trendout$species, SPECIES)]
mapdata$fullspec<-fullnames[match(mapdata$species, SPECIES)]

## add colour
trendout<-trendout %>%
  mutate(col=ifelse(lcl<0,ifelse(ucl<0,"darkred","black"),ifelse(ucl>0,"forestgreen","black"))) %>%
  mutate(col=ifelse(species=="CAEL","darkred",col))


## add numeric year  
annestimates$Year<-str_replace_all(annestimates$parameter,pattern="[^[:alnum:]]", replacement="")
annestimates$Year<-str_replace_all(annestimates$Year,pattern="totalN", replacement="")
annestimates$Year<-as.numeric(annestimates$Year)+2010

## define dimensions of arrays
nsites<-length(unique(siteCov$Point))
nyears<-length(unique(countdata$year))


## create simple feature

mapdata_sf<-mapdata %>%
  st_as_sf(coords = c("Eastings", "Northings"), crs=2004) %>%
  st_transform(4326)                        ## Montserrat is EPSG 4604 or 2004



#.....................................................................................
# 3. CREATE BASIC PLOTS FROM OUTPUT ---------
#.....................................................................................


#.....................................................................................
## 3.1 PLOT FOR ABUNDANCE TREND -------					
#.....................................................................................

annestimates %>% 
  filter(Year!=2020) %>%
  mutate(col = as.factor(trendout$col[match(species,trendout$species)])) %>%
  
ggplot()+
	geom_line(aes(x=Year, y=mean,col=col), linewidth=1)+
	facet_wrap(~fullspec, ncol=2, scales="free_y")+
geom_point(aes(x=Year, y=mean,col=col), size=2)+
#geom_ribbon(data=annestimates,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
geom_errorbar(aes(x=Year, ymin=lcl,ymax=ucl,col=col), width=.1) +

## remove the legend
theme(legend.position="none")+
guides(scale="none",fill=FALSE)+
theme(legend.title = element_blank())+
theme(legend.text = element_blank())+

## format axis ticks
scale_x_continuous(name="Year", breaks=seq(2011,YEAR,2), labels=as.character(seq(2011,YEAR,2)))+
#scale_y_continuous(name="Number of Birds at 67 Sampling Points", breaks=seq(0,4000,500), labels=as.character(seq(0,4000,500)))+
  ylab(sprintf("Number of birds at %i sampling points",nsites)) +
  scale_color_manual(values = c("black","darkred", "forestgreen"))+
## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=18, color="black"),
        axis.title=element_text(size=18),
        strip.text.x=element_text(size=18, color="black"),
	  axis.title.y=element_text(margin=margin(0,20,0,0)),
        strip.background=element_rect(fill="white", colour="black"))





#.....................................................................................
## 3.2 PROVIDE BASIC SUMMARY REPORT FOR ANNUAL NEWS ARTICLE  -------					
#.....................................................................................

surveys<-obsCov %>% filter(year==YEAR)
birds<-countdata %>% filter(year==YEAR) %>% select(-Time,-Rain,-Wind,-day,-time,-activity,-Date,-Time,-VisitID) %>%
	gather(key=Species, value=N, -year,-Point,-Count) %>%
	mutate(Point=as.integer(as.character(Point)))
summary<-surveys %>% select(year, Point, Count, Date) %>%
  left_join(birds, by=c('year','Point','Count')) %>%
  group_by(Count, Species) %>%
  summarise(N=sum(N, na.rm=T))

totals<-summary %>% group_by(Count) %>%
  summarise(N=sum(N), n_spec=length(unique(Species)))

table1<-summary %>% spread(key=Count, value=N, fill = 0) %>%
  filter(!is.na(Species)) %>%
  mutate(Species=species$Species[match(Species,species$SpeciesCode)]) %>%
  arrange(desc(`1`))

table2<-trendout %>% filter(parameter %in% c('trend','trend2')) %>% # trend2 is the quadratic trend estimate
  mutate(dir=ifelse(lcl<0,ifelse(ucl<0,"decrease","stable"),ifelse(ucl>0,"increase","stable"))) %>%  
  mutate(dir=ifelse(parameter=="trend", dir,
                    ifelse(lcl<0,ifelse(ucl<0,"concave","flat"),ifelse(ucl>0,"convex","flat")))) %>%  ## ADJUST THIS FOR trend2 separately: https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faqhow-do-i-interpret-the-sign-of-the-quadratic-term-in-a-polynomial-regression/
  mutate(conf=paste(round(mean, 2)," (",round(lcl,2)," - ",round(ucl,2),")", sep="")) %>%
  mutate(BayesP = round(BayesP,2)) %>%
  select(fullspec,dir,conf,BayesP)


#.....................................................................................
## 3.3 PROVIDE BASIC MAPS OF RECENT ABUNDANCE  -------					
#.....................................................................................

abundance_sf<-mapdata_sf %>%
  filter(Year==YEAR) %>%
  group_by(Point, geometry) %>%
  summarise(N=sum(median))

tmap_mode("view")
tm_basemap(c(StreetMap = "OpenStreetMap",
             TopoMap = "OpenTopoMap")) +
  # tm_shape(basemap)+
  #   tm_rgb()+
  tm_shape(abundance_sf)  +
  tm_symbols(col="N", palette = c("lightgreen", "darkred"),
             size=0.8, alpha=0.7,
             title.col = "Total N")







#.....................................................................................
# 4. CREATE OUTPUT REPORT WITH MARKDOWN ---------
#.....................................................................................


### create HTML report for overall summary report
 Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")

rmarkdown::render('Annual_abundance_report_modelled_with_maps.Rmd',
                  output_file = sprintf("Montserrat_ForestBird_AnnualSummary%s.html",YEAR),
                  output_dir = './output')






