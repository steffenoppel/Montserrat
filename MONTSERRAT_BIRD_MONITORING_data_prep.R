######################################################################################
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  DATA PREPARATION FOR MULTI-YEAR ANALYSIS   ############################
#############  steffen.oppel@gmail.com      ##########################################
######################################################################################

### generic data extraction and preparation for analysis
### adapted from spOccupancy script, but removed modelling and only retained data preparation

YEAR<-2023									## ENTER THE YEAR OF THE MOST RECENT MONITORING
SPECIES<-c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA')				## ENTER THE SPECIES FOR WHICH YOU WANT ANALYSES TO BE DONE


######################################################################################
#############  Load required packages       ##########################################
######################################################################################

library(RODBC)
library(tidyverse)
library(sf)
library(lubridate)
sf::sf_use_s2(FALSE)
filter<-dplyr::filter
select<-dplyr::select
rename<-dplyr::rename





######################################################################################
#############  Import data from the queries in the Access Database        ############
######################################################################################
#setwd("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring")
#load("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_RODBC_DATA_INPUT.RData")
#setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring")
#load("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_RODBC_DATA_INPUT.RData")

setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data")
MTOR_count <- odbcConnectAccess2007('Montserrat_Birds_2023.accdb')
birds <- sqlQuery(MTOR_count, "SELECT * FROM Total_count_allpoints_by_species_by_year_vertical")
obsCov <- sqlQuery(MTOR_count, "SELECT * FROM obsCovariates_simple")
siteCov <- sqlQuery(MTOR_count, "SELECT * FROM SiteCov")
species <- sqlQuery(MTOR_count, "SELECT * FROM lkSpecies")
odbcClose(MTOR_count)

nyears<-YEAR-2010


######################################################################################
#############  READ SHAPEFILES         #############
######################################################################################

setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\GIS")
#setwd("C:\\STEFFEN\\RSPB\\Montserrat\\GIS")


elevation_map <- sf::read_sf("dem_shp.shp") 
elevation_map %>% 
  ggplot() +
  geom_sf()



######################################################################################
#############  REMOVE POINTS THAT ARE NOT INDEPENDENT         ########################
######################################################################################
## these points are too close to other points in the forest and may lead to double-counts
## this is arguably not critical for trend monitoring, but was for total population size extrapolation
 
removal<-c(99,76)																			## removes only two points that cause error in JAGS

siteCov<-siteCov[!(siteCov$Point %in% removal),]
obsCov<-obsCov[!(obsCov$Point %in% removal),]
birds<-birds[!(birds$Point %in% removal),]





######################################################################################
#############  RELABEL HABITAT VARIABLE NAMES AND FILL IN ZEROs          #############
######################################################################################


names(siteCov)[4:9]<-c('hab','dbh','dist','treeheight','alt','slope')
siteCov[is.na(siteCov)]<-0
siteCov$ridge<-ifelse(siteCov$Location=="Ridge",1,0)				## convert location to a two-level factor

### scale the X and Y coordinates
siteCov<-siteCov %>%
  mutate(x=Eastings-min(Eastings), y=Northings-min(Northings)) %>%
  dplyr::select(Point,x,y,hab,dbh,dist,treeheight,alt,slope,ridge,Eastings,Northings)


head(siteCov)



######################################################################################
#############  CREATE NEW NUMERIC VARIABLES FOR RAIN AND WIND AND DAY    #############
######################################################################################

obsCov$rain<- 1
obsCov$rain[is.na(obsCov$Rain)] <- 0

obsCov$wind<-as.ordered(match(obsCov$Wind,levels(obsCov$Wind)))
obsCov$wind<-ifelse(obsCov$wind<3,0,1)					### convert wind into two-level factor yes and no
obsCov$Wind<-NULL
names(obsCov)[c(7,8)]<- c("obs", "skill")
head(obsCov)


### CREATE A NUMERIC VARIABLE FOR THE DAY
obsCov$Year<-as.numeric(format(strptime(obsCov$Date, format="%Y"),format="%Y"))
obsCov$Day<-1
for (y in 2011:YEAR){
startdate<-as.Date(sprintf("%i-03-01 00:00",y))  				#set a startdate, here 1 March (must be smaller than first survey day)
obsCov$Day[obsCov$Year==y]<-as.numeric(as.Date(obsCov$Date[obsCov$Year==y])-startdate) 	#calculate the difference between the start date and the survey date
}

obsCov$Date<-NULL

starttime<-min(as.numeric(obsCov$Time-60)) 				#sets the first survey time-1 min as the min time (here 5:49 am)
obsCov$time<-(as.numeric(obsCov$Time)- starttime)/60  		#calculate the difference between 5:49 am and the time of the survey in minutes (default is in seconds)
obsCov$Time<-NULL

head(obsCov)


################ CREATE A BIRD ACTIVITY COVARIATE (ALL SPECIES) FOR EACH POINT COUNT ##########
#### NEEDS TO BE SCALED FOR EACH POINT TO OVERCOME SPATIAL DIFFERENCES #########

head(birds)
pointmax<- birds %>% group_by(Year,Count,Point) %>%
	summarise(activity=sum(SumOfNumber1, na.rm=T)) %>%
	group_by(Point) %>%
	summarise(max_act=max(activity))

activity<- birds %>% group_by(Year,Count,Point) %>%
	summarise(activity=sum(SumOfNumber1, na.rm=T)) %>%
	inner_join(pointmax) %>%
	mutate(ACT=activity/max_act) %>%
	select(Year,Count,Point,ACT)
names(activity)[1]<-'year'

obsCov<-merge(obsCov,activity, by=c('year','Point','Count'), all.x=T)
head(obsCov)


######################################################################################
#############  ADD ZERO COUNTS FOR TARGET SPECIES FROM CUT DOWN POINT 15 #############
######################################################################################
add<-data.frame(Year=2015, Species=rep(SPECIES,each=3), Point='15', Count=rep(1:3,length(SPECIES)), SumOfNumber1=0)
birds<-rbind(birds,add)

head(obsCov)
add<-data.frame(year=2015, Point='15', Count=rep(1:3), Rain=NA, obs=2, skill=1, rain=0, wind=0, Year=2015, Day=37, time=150, ACT=0.2)
obsCov<-rbind(obsCov,add)



######################################################################################
#############  CREATE BLANK ARRAYS FOR INPUT DATA          ###########################
######################################################################################

COUNTDATA<-data.frame(Year=rep(seq(2011,YEAR,1),length(SPECIES)*3*length(unique(siteCov$Point))),
			Species=rep(SPECIES, each=length(seq(2011:YEAR))*3*length(unique(siteCov$Point))),
			Point=rep(rep(unique(siteCov$Point),each=length(seq(2011:YEAR))), length(SPECIES)*3),
			Count=rep(rep(1:3, each=length(seq(2011:YEAR))*length(unique(siteCov$Point))),length(seq(SPECIES))))

SURVEYDATA<-data.frame(Year=rep(seq(2011,YEAR,1),3*length(unique(siteCov$Point))),
			Point=rep(rep(unique(siteCov$Point),each=length(seq(2011:YEAR))), 3),
			Count=rep(1:3, each=length(seq(2011:YEAR))*length(unique(siteCov$Point))))




######################################################################################
#############  FILL IN THE ACTUAL DATA FOR THE FULL ARRAY     ########################
######################################################################################

COUNTDATA<-merge(COUNTDATA, birds, by=c("Species", "Year", "Point", "Count"), all.x=T)
names(COUNTDATA)[5]<-"N"
COUNTDATA$N[is.na(COUNTDATA$N)]<-0										### sets all NA to 0, which includes counts that did not take place
SURVEYDATA<-merge(SURVEYDATA, obsCov, by=c("Year", "Point", "Count"), all.x=T)

SURVEYDATA<-SURVEYDATA %>%
  arrange(Point,Year,Count)




###############################################################################
############## CREATE OCCURRENCE COVARIATE DATA INPUT LIST   ######################
###############################################################################
nsites<-length(unique(siteCov$Point))
# occ.covs is a list of variables included in the occurrence portion of the model
# Each list element is a different occurrence covariate, which can be site level or site/primary timer period level
# Site-level covariates are specified as a vector of length J while
# site/primary time period level covariates are specified as a matrix with rows corresponding to sites and columns correspond to primary time periods.
occ.covs<-list()
siteCov<-siteCov %>% arrange(Point)
occ.covs$hab<-siteCov$hab
occ.covs$alt<-scale(siteCov$alt)
occ.covs$treeheight<-scale(siteCov$treeheight)
occ.covs$trend<-matrix(data=scale(1:nyears), ncol=nyears, nrow=nsites, byrow=T)

coords<-as.matrix(siteCov %>% select(x,y))


###############################################################################
############## CREATE DETECTION COVARIATE DATA INPUT LIST  ###############
###############################################################################
# det.covs is a list of variables included in the detection portion of the model
# each list element corresponding to an individual variable
# In addition to site-level and/or site/primary time period-level, detection covariates can also be observational-level
# Observation-level covariates are specified as a three-dimensional array with first dimension corresponding to sites, second dimension corresponding to primary time period, third dimension corresponding to replicate.
## SORT THE TABLE SO IT HAS THE SAME ORDER AS THE BIRD DATA

obsCov<-obsCov[order(obsCov$Point,obsCov$Year, obsCov$Count, decreasing=F),] 
head(obsCov)


### STANDARDIZE COVARIATES FOR WINBUGS

meant<-mean(SURVEYDATA$time, na.rm = TRUE)
sdt<-sd(SURVEYDATA$time, na.rm = TRUE)
SURVEYDATA$time<-(SURVEYDATA$time-meant)/sdt

meant<-mean(SURVEYDATA$Day, na.rm = TRUE)
sdt<-sd(SURVEYDATA$Day, na.rm = TRUE)
SURVEYDATA$Day<-(SURVEYDATA$Day-meant)/sdt


### only needs standardisation if measured in mm, not as 0/1 variable
#meant<-mean(SURVEYDATA$rain, na.rm = TRUE)
#sdt<-sd(SURVEYDATA$rain, na.rm = TRUE)
#SURVEYDATA$rain<-(SURVEYDATA$rain-meant)/sdt
nreps<-length(unique(SURVEYDATA$Count))
nyears<-length(unique(SURVEYDATA$Year))


### create array for each covariate

Day<-array(NA, dim=c(nsites,nyears, nreps))
#wind<-array(NA, dim=c(nsites,nyears, nreps))
time<-array(NA, dim=c(nsites,nyears, nreps))
ACT<-array(NA, dim=c(nsites,nyears, nreps))				## REPLACED ON 2 MAY WITH RAINFALL AMOUNT

### fill in array for each covariate
for (c in 1:nreps){
  obsC<-SURVEYDATA %>% dplyr::filter(Count==c) %>% #filter(Year==2023, Point = 45)
    mutate(PointCount=paste(Point,Count,sep="_")) %>%
    ungroup() %>%
    select(PointCount,Year,ACT,Day,time)
  
  x<-obsC %>%
    select(PointCount,Year,ACT) %>% 
    spread(key=Year, value='ACT') %>%
    arrange(PointCount)
  x2<-as.matrix(x[,2:(nyears+1)], dimnames=NULL)
  ACT[,,c]<-x2
  
  x<-obsC %>%
    select(PointCount,Year,Day) %>%
    spread(key=Year, value='Day') %>%
    arrange(PointCount)
  x2<-as.matrix(x[,2:(nyears+1)], dimnames=NULL)
  Day[,,c]<-x2
  
  x<-obsC %>%
    select(PointCount,Year,time) %>%
    spread(key=Year, value='time') %>%
    arrange(PointCount)
  x2<-as.matrix(x[,2:(nyears+1)], dimnames=NULL)
  time[,,c]<-x2
}


###############################################################################
####   REPLACE ALL NA IN COVARIATES otherwise "undefined node" error    #######
###############################################################################

for (d in 1:nyears){							### replace missing dates with mean for each survey round in each year
  ACT[is.na(ACT[,d,1]),d,1]<-mean(ACT[,d,1], na.rm=T)
  ACT[is.na(ACT[,d,2]),d,2]<-mean(ACT[,d,2], na.rm=T)
  ACT[is.na(ACT[,d,3]),d,3]<-mean(ACT[,d,3], na.rm=T)
  Day[is.na(Day[,d,1]),d,1]<-mean(Day[,d,1], na.rm=T)
  Day[is.na(Day[,d,2]),d,2]<-mean(Day[,d,2], na.rm=T)
  Day[is.na(Day[,d,3]),d,3]<-mean(Day[,d,3], na.rm=T)
  time[is.na(time[,d,1]),d,1]<-mean(time[,d,1], na.rm=T)
  time[is.na(time[,d,2]),d,2]<-mean(time[,d,2], na.rm=T)
  time[is.na(time[,d,3]),d,3]<-mean(time[,d,3], na.rm=T)
}


det.covs<-list()
det.covs$Day<-Day
det.covs$time<-time
det.covs$ACT<-ACT



######################################################################################
#############  START THE LOOP OVER EVERY SPECIES          ############################
######################################################################################

for (s in SPECIES){
  
  ######################################################################################
  #############  TAKE SUBSET OF DATA FOR SELECTED SPECIES AND SORT THE TABLES    ###################
  ######################################################################################
  
  bird_s<-COUNTDATA %>% filter(Species==s) %>%
    arrange(Point,Year,Count)
  

  
  ### RE-INTRODUCE THE NAs for COUNTS THAT DID NOT TAKE PLACE #####
  bird_s$N[is.na(SURVEYDATA$time)]<-NA
  
  
  ###############################################################################
  ############## CREATE BIRD DATA INPUT MATRIX   ################################
  ###############################################################################
  ## y is a three-dimensional array with first dimension equal to the number of sites (J)
  ## second dimension equal to the maximum number of primary time periods (i.e., years or seasons)
  ## third dimension equal to the maximum number of replicates at a given site (3 in Montserrat)
  
  ## note that for single season models you will only need two dimensions!
  

  nyears<-length(unique(bird_s$Year))
  nreps<-length(unique(bird_s$Count))
  
  
  
  ### create array to be filled with data
  BIRD.y<-array(NA, dim=c(nsites,nyears,nreps))

  #### CAST THE MOLTEN DATA FRAME INTO MATRIX WITH 1 COLUMN PER YEAR and fill in array
  for (c in 1:nreps){
    x<- bird_s %>% filter(Count==c) %>%
      mutate(PointCount=paste(Point,Count,sep="_")) %>%
      mutate(N=ifelse(N>0,1,0)) %>%
      mutate(N=ifelse(is.na(N),0,N)) %>%
      spread(key=Year, value=N) %>%
      arrange(PointCount) %>%
      dplyr::select(-Species,-Point,-Count)
    BIRD.y[,,c]<-as.matrix(x[,2:(nyears+1)], dimnames=NULL)
    #Nst[,y]<-apply(dis[,2:4], MARGIN=1, FUN=max, na.rm=T)+1
  }
  
  BIRD.y[,nyears-(YEAR-2020),]<-NA   ### COVID-19 year (2020) with no surveys
  

  


##### END THE LOOP ACROSS SPECIES ####
}


