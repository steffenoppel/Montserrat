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

#setwd("S:/ConSci/DptShare/SteffenOppel/RSPB/Montserrat/Analysis/Population_status_assessment/AnnualMonitoring")
setwd("C:\\STEFFEN\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring")

#annestimates<-read.table("Annual_estimates2017.csv", header=T, sep=",")
#trendout<-read.table("Trend_estimates2017.csv", header=T, sep=",")


annestimates$fullspec<-fullnames[match(annestimates$species, SPECIES)]
trendout$fullspec<-fullnames[match(trendout$species, SPECIES)]



################ PLOT FOR ABUNDANCE TREND ####################
trendout<-trendout %>%
  mutate(col=ifelse(lcl<0,ifelse(ucl<0,"darkred","black"),ifelse(ucl>0,"forestgreen","black"))) %>%
  mutate(col=ifelse(species=="CAEL","darkred",col))
annestimates %>% 
  mutate(Year=rep(seq(2011,2024), length(allout))) %>%
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
  scale_x_continuous(name="Year", breaks=seq(2011,2023,2), labels=as.character(seq(2011,2023,2)))+
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

ggsave("Montserrat_ForestBird_Trends_2024.pdf", width=13, height=16)






######################################################################################
#############  PROVIDE BASIC SUMMARY REPORT FOR ANNUAL NEWS ARTICLE        ##################
######################################################################################


surveys2024<-obsCov %>% filter(year==2024)
birds2024<-countdata %>% filter(year==2024) %>% select(-Time,-Rain,-Wind,-day,-time,-activity,-Date,-Time,-VisitID) %>%
  gather(key=Species, value=N, -year,-Point,-Count) %>%
  mutate(Point=as.integer(as.character(Point)))
summary2024<-surveys2024 %>% select(year, Point, Count, Date) %>%
  left_join(birds2024, by=c('year','Point','Count')) %>%
  group_by(Count, Species) %>%
  summarise(N=sum(N, na.rm=T))

totals2024<-summary2024 %>% group_by(Count) %>%
  summarise(N=sum(N), n_spec=length(unique(Species)))

table1<-summary2024 %>% spread(key=Count, value=N, fill = 0) %>%
  filter(!is.na(Species)) %>%
  mutate(Species=species$Species[match(Species,species$SpeciesCode)]) %>%
  arrange(desc(`1`))

table2<-trendout %>%
  mutate(dir=ifelse(lcl<0,ifelse(ucl<0,"decrease","stable"),ifelse(ucl>0,"increase","stable"))) %>%
  mutate(dir=ifelse(species=="CAEL","(decrease)",dir)) %>%
  mutate(conf=paste(round(mean,3)," (",round(lcl,3)," - ",round(ucl,3),")", sep="")) %>%
  select(fullspec,dir,conf,BayesP,Rhat) %>%
  arrange(desc(BayesP))



######################################################################################
#############  SIMPLE plot for the key species     ########################
######################################################################################

### ALTERNATIVE IF MODELS DO NOT CONVERGE

# ggplot()+
#   geom_line(data=summary, aes(x=Year, y=mean), linewidth=1)+
#   facet_wrap(~Species, ncol=2, scales="free_y")+
#   geom_point(data=summary, aes(x=Year, y=mean), size=1.5,col='black')+
#   #geom_ribbon(data=annestimates,aes(x=Year, ymin=lower95CI,ymax=upper95CI),alpha=0.2)+
#   geom_errorbar(data=summary,aes(x=Year, ymin=(mean-0.5*sd),ymax=(mean+0.5*sd)), colour="black", width=.1) +
#
#   ## remove the legend
#   theme(legend.position="none")+
#   guides(fill=FALSE)+
#   theme(legend.title = element_blank())+
#   theme(legend.text = element_blank())+
#
#   ## format axis ticks
#   scale_x_continuous(name="Year", breaks=seq(2011,2023,1), labels=as.character(seq(2011,2023,1)))+
#   #scale_y_continuous(name="Number of Birds at 67 Sampling Points", breaks=seq(0,4000,500), labels=as.character(seq(0,4000,500)))+
#   ylab("Number of Birds per Sampling Point") +
#
#   ## beautification of the axes
#   theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text=element_text(size=18, color="black"),
#         axis.title=element_text(size=18),
#         strip.text.x=element_text(size=18, color="black"),
#         axis.title.y=element_text(margin=margin(0,20,0,0)),
#         strip.background=element_rect(fill="white", colour="black"))




######## THIS WILL ONLY WORK ONCE THE ANNUAL MODELS HAVE BEEN RUN ##

# annestimates<-read.table("Annual_estimates2023.csv", header=T, sep=",")
# fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler","Antillean Crested Hummingbird","Purple-throated Carib","Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher",
# 			"Scaly-naped Pigeon","Caribbean Elaenia", "Bananaquit")
# annestimates$fullspec<-fullnames[match(annestimates$species, SPECIES)]

### create HTML report for overall summary report
Sys.setenv(RSTUDIO_PANDOC="C:/Users/Inge Oppel/AppData/Local/Pandoc")
# Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")
#
rmarkdown::render('C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Annual_abundance_report.Rmd',
                  output_file = "Montserrat_ForestBird_AnnualSummary2024.html",
                  output_dir = 'C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring')

rmarkdown::render('C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat\\Annual_abundance_report_modelled.Rmd',
                  output_file = "Montserrat_ForestBird_AnnualSummary2024.html",
                  output_dir = 'C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat')







