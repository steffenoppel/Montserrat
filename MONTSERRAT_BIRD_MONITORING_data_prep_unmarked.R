######################################################################################
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  DATA PREPARATION FOR MULTI-YEAR ANALYSIS   ############################
#############  steffen.oppel@gmail.com      ##########################################
######################################################################################
### data preparation to facilitate analysis in package unmarked
### commented out the main data extraction process, load saved workspace and simply reformat data for package unmarked

YEAR<-2023									## ENTER THE YEAR OF THE MOST RECENT MONITORING
SPECIES<-c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA')				## ENTER THE SPECIES FOR WHICH YOU WANT ANALYSES TO BE DONE


######################################################################################
#############  Load required packages       ##########################################
######################################################################################

library(RODBC)
library(reshape)
library(tidyverse)
library(lubridate)
filter<-dplyr::filter
select<-dplyr::select
rename<-dplyr::rename
library(unmarked)
library(here)

# ######################################################################################
# #############  Import data from the queries in the Access Database        ############
# ######################################################################################
# #setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Montserrat\\Raw_Data")
# #setwd("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data")
# #setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data")
# #system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data\\RODBC_annual_count_import.r")), wait = TRUE, invisible = FALSE)
# #system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data\\RODBC_annual_count_import.r")), wait = TRUE, invisible = FALSE)
# 
# 
# #setwd("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring")
# #load("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_RODBC_DATA_INPUT.RData")
# #setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data")
# MTOR_count <- odbcConnectAccess2007(here('Montserrat_Birds_2023.accdb'))
# birds <- sqlQuery(MTOR_count, "SELECT * FROM Total_count_allpoints_by_species_by_year_vertical")
# obsCov <- sqlQuery(MTOR_count, "SELECT * FROM obsCovariates_simple")
# siteCov <- sqlQuery(MTOR_count, "SELECT * FROM SiteCov")
# species <- sqlQuery(MTOR_count, "SELECT * FROM lkSpecies")
# odbcClose(MTOR_count)
# nyears<-YEAR-2010
# 
# 
# 
# 
# ######################################################################################
# #############  PROVIDE BASIC SUMMARY FOR ANNUAL NEWS ARTICLE        ##################
# ######################################################################################
# 
# 
# 
# surveys2023<-obsCov %>% filter(year==2023) %>% rename(Year=year)
# birds2023<-birds %>% filter(Year==2023) %>% rename(N=SumOfNumber1) %>% mutate(Point=as.integer(as.character(Point))) %>% filter(!Species %in% c('UNK','NA'))
# summary2023<-surveys2023 %>% select(Year, Point, Count, Date) %>%
#   left_join(birds2023, by=c('Year','Point','Count')) %>%
#   group_by(Count, Species) %>%
#   summarise(N=sum(N, na.rm=T))
#   
# totals2023<-summary2023 %>% group_by(Count) %>%
#   summarise(N=sum(N), n_spec=length(unique(Species)))
# 
# table1<-summary2023 %>% spread(key=Count, value=N, fill = 0) %>%
#   filter(!is.na(Species)) %>%
#   mutate(Species=species$Species[match(Species,species$SpeciesCode)]) %>%
#   arrange(desc(`1`))
# 
# 
# 
# 
# ######################################################################################
# #############  REMOVE POINTS THAT ARE NOT INDEPENDENT         ########################
# ######################################################################################
# ## these points are too close to other points in the forest and may lead to double-counts
# ## this is arguably not critical for trend monitoring, but was for total population size extrapolation
#  
# #removal<-c(6,8,109,46,103,92,122,125,72,69,70,90,121,115,124,126,1,2,16,34,10,19,45,50,73,75,81,113, 99,76)		## exhaustive list
# #removal<-c(6,8,109,46,103,92,122,125,72,69,70,90,121,115,124,126,1,2,16,34,10,19,45,50,73,75,81,113, 176:300)		## exhaustive list for ANDY (includes 76,99, 174, and 175)
# removal<-c(99,76)																			## removes only two points that cause error in JAGS
# 
# siteCov<-siteCov[!(siteCov$Point %in% removal),]
# obsCov<-obsCov[!(obsCov$Point %in% removal),]
# birds<-birds[!(birds$Point %in% removal),]
# 
# 
# 
# 
# 
# ######################################################################################
# #############  RELABEL HABITAT VARIABLE NAMES AND FILL IN ZEROs          #############
# ######################################################################################
# 
# 
# names(siteCov)[2:7]<-c('hab','dbh','dist','treeheight','alt','slope')
# siteCov[is.na(siteCov)]<-0
# siteCov$ridge<-ifelse(siteCov$Location=="Ridge",1,0)				## convert location to a two-level factor
# head(siteCov)
# 
# 
# 
# ######################################################################################
# #############  CREATE NEW NUMERIC VARIABLES FOR RAIN AND WIND AND DAY    #############
# ######################################################################################
# 
# obsCov$rain<- 1
# obsCov$rain[is.na(obsCov$Rain)] <- 0
# 
# obsCov$wind<-as.ordered(match(obsCov$Wind,levels(obsCov$Wind)))
# obsCov$wind<-ifelse(obsCov$wind<3,0,1)					### convert wind into two-level factor yes and no
# obsCov$Wind<-NULL
# names(obsCov)[c(7,8)]<- c("obs", "skill")
# head(obsCov)
# 
# 
# ### CREATE A NUMERIC VARIABLE FOR THE DAY
# obsCov$Year<-as.numeric(format(strptime(obsCov$Date, format="%Y"),format="%Y"))
# obsCov$Day<-1
# for (y in 2011:YEAR){
# startdate<-as.Date(sprintf("%i-03-01 00:00",y))  				#set a startdate, here 1 March (must be smaller than first survey day)
# obsCov$Day[obsCov$Year==y]<-as.numeric(as.Date(obsCov$Date[obsCov$Year==y])-startdate) 	#calculate the difference between the start date and the survey date
# }
# 
# obsCov$Date<-NULL
# 
# starttime<-min(as.numeric(obsCov$Time-60)) 				#sets the first survey time-1 min as the min time (here 5:49 am)
# obsCov$time<-(as.numeric(obsCov$Time)- starttime)/60  		#calculate the difference between 5:49 am and the time of the survey in minutes (default is in seconds)
# obsCov$Time<-NULL
# 
# head(obsCov)
# 
# 
# ################ CREATE A BIRD ACTIVITY COVARIATE (ALL SPECIES) FOR EACH POINT COUNT ##########
# #### NEEDS TO BE SCALED FOR EACH POINT TO OVERCOME SPATIAL DIFFERENCES #########
# 
# head(birds)
# pointmax<- birds %>% group_by(Year,Count,Point) %>%
# 	summarise(activity=sum(SumOfNumber1, na.rm=T)) %>%
# 	group_by(Point) %>%
# 	summarise(max_act=max(activity))
# 
# activity<- birds %>% group_by(Year,Count,Point) %>%
# 	summarise(activity=sum(SumOfNumber1, na.rm=T)) %>%
# 	inner_join(pointmax) %>%
# 	mutate(ACT=activity/max_act) %>%
# 	select(Year,Count,Point,ACT)
# names(activity)[1]<-'year'
# 
# obsCov<-merge(obsCov,activity, by=c('year','Point','Count'), all.x=T)
# head(obsCov)
# 
# 
# ######################################################################################
# #############  ADD ZERO COUNTS FOR TARGET SPECIES FROM CUT DOWN POINT 15 #############
# ######################################################################################
# add<-data.frame(Year=2015, Species=rep(SPECIES,each=3), Point='15', Count=rep(1:3,length(SPECIES)), SumOfNumber1=0)
# birds<-rbind(birds,add)
# 
# head(obsCov)
# add<-data.frame(year=2015, Point='15', Count=rep(1:3), Rain=NA, obs=2, skill=1, rain=0, wind=0, Year=2015, Day=37, time=150, ACT=0.2)
# obsCov<-rbind(obsCov,add)
# 
# 
# 
# ######################################################################################
# #############  CREATE BLANK ARRAYS FOR INPUT DATA          ###########################
# ######################################################################################
# 
# COUNTDATA<-data.frame(Year=rep(seq(2011,YEAR,1),length(SPECIES)*3*length(unique(siteCov$Point))),
# 			Species=rep(SPECIES, each=length(seq(2011:YEAR))*3*length(unique(siteCov$Point))),
# 			Point=rep(rep(unique(siteCov$Point),each=length(seq(2011:YEAR))), length(SPECIES)*3),
# 			Count=rep(rep(1:3, each=length(seq(2011:YEAR))*length(unique(siteCov$Point))),length(seq(SPECIES))))
# 
# SURVEYDATA<-data.frame(Year=rep(seq(2011,YEAR,1),3*length(unique(siteCov$Point))),
# 			Point=rep(rep(unique(siteCov$Point),each=length(seq(2011:YEAR))), 3),
# 			Count=rep(1:3, each=length(seq(2011:YEAR))*length(unique(siteCov$Point))))
# 
# 
# 
# 
# ######################################################################################
# #############  FILL IN THE ACTUAL DATA FOR THE FULL ARRAY     ########################
# ######################################################################################
# 
# COUNTDATA<-merge(COUNTDATA, birds, by=c("Species", "Year", "Point", "Count"), all.x=T)
# names(COUNTDATA)[5]<-"N"
# COUNTDATA$N[is.na(COUNTDATA$N)]<-0										### sets all NA to 0, which includes counts that did not take place
# SURVEYDATA<-merge(SURVEYDATA, obsCov, by=c("Year", "Point", "Count"), all.x=T)
# 
# 
# ######################################################################################
# #############  SAVE WORKSPACE TO RUN ANALYSIS IN CLOUD OR ON 64-bit SYSTEM    ########
# ######################################################################################
# #save.image("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_ANNUAL_DATA_INPUT2023.RData")
# save.image("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_ANNUAL_DATA_INPUT2023.RData")
load(here("MONTSERRAT_ANNUAL_DATA_INPUT2023.RData"))
ls()





#################################################################################################
### CREATE AN UNMARKED DATA FRAME FOR FUNCTIONS occu   ###########
#################################################################################################

### CODE ONLY PROVIDED AS EXAMPLE FOR ONE SPECIES AND YEAR
#####  CREATE INPUT DATA FROM 2023 ####
occdata<-COUNTDATA %>% dplyr::filter(Species=="MTOR") %>%
  filter(Year==2023) %>%
  mutate(occ=ifelse(N>0,1,0)) %>%   # convert counts to detection/non-detection
  dplyr::select(-Species,-Year,-N) %>%
  spread(key=Count, value=occ) %>%
  arrange(Point)

siteCovOccu<-siteCov %>% mutate(Point=as.numeric(Point)) %>%
  arrange(Point)

obsCovOccu<-COUNTDATA %>% dplyr::filter(Species=="MTOR") %>%
  filter(Year==2023) %>%
  left_join((obsCov %>% mutate(Point=as.numeric(Point))), by=c("Point","Count","Year")) %>%
  select(-N,-Species,-Rain,-skill) %>%
  arrange(Point,Count)
dim(occdata)
dim(siteCovOccu)
dim(obsCovOccu)/3

OCCU_UMF <- unmarkedFrameOccu(occdata[,2:4],
                              siteCovs = siteCovOccu,
                              obsCovs = obsCovOccu)		# create data frame for occupancy models
summary(OCCU_UMF)

siteCovs(OCCU_UMF)[c(2,3,5:17,19,21:38)] <- scale(siteCovs(OCCU_UMF)[c(2,3,5:17,19,21:38)])     	### standardize numeric covariate to avoid numerical fitting problems
obsCovs(OCCU_UMF)[c(5,6,8,9,10)] <- scale(obsCovs(OCCU_UMF)[c(5,6,8,9,10)])     					### standardize numeric covariate to avoid numerical fitting problems

test.model <- occu(~ ACT + time + I(time^2) + Day ~ Altitude, OCCU_UMF, se=T)
summary(test.model)



#################################################################################################
### CREATE AN UNMARKED DATA FRAME FOR FUNCTIONS colext   ###########
#################################################################################################
### CODE ONLY PROVIDED AS EXAMPLE FOR ONE SPECIES

numPrimary=YEAR-2010
occdataCOLEXT<-COUNTDATA %>% dplyr::filter(Species=="MTOR") %>%
  mutate(occ=ifelse(N>0,1,0)) %>%   # convert counts to detection/non-detection
  mutate(Season=paste(Year,Count,sep="_")) %>%
  dplyr::select(-Species,-N,-Year,-Count) %>%
  spread(key=Season, value=occ) %>%
  arrange(Point)

siteCovColext<-siteCov %>% mutate(Point=as.numeric(Point)) %>%
  arrange(Point)

obsCovColext<-COUNTDATA %>% dplyr::filter(Species=="MTOR") %>%
  left_join((obsCov %>% mutate(Point=as.numeric(Point))), by=c("Point","Count","Year")) %>%
  select(-N,-Species,-Rain,-skill) %>%
  arrange(Point,Year,Count)
dim(occdataCOLEXT)
dim(siteCovColext)
dim(obsCovColext)/(3*numPrimary)

COLEXT_UMF <- unmarkedMultFrame(occdataCOLEXT[,2:40],
                                numPrimary=numPrimary,
                              siteCovs = siteCovColext,
                              obsCovs = obsCovColext)		# create data frame for occupancy models
summary(COLEXT_UMF)

siteCovs(COLEXT_UMF)[c(2,3,5:17,19,21:38)] <- scale(siteCovs(COLEXT_UMF)[c(2,3,5:17,19,21:38)])     	### standardize numeric covariate to avoid numerical fitting problems
obsCovs(COLEXT_UMF)[c(5,6,8,9,10)] <- scale(obsCovs(COLEXT_UMF)[c(5,6,8,9,10)])     					### standardize numeric covariate to avoid numerical fitting problems


test.model<-colext(psiformula= ~treeheight, gammaformula =  ~ Altitude, epsilonformula = ~ Altitude,    pformula = ~ ACT + time + Day, COLEXT_UMF, se=TRUE)
