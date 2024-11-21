#### script for automating the yearly analysis of centre hills forest bird monitoring on Montserrat ####
#### script written by Filibert Heim, filibert.heim@posteo.de, but compomonents of this script most importantly the N-mix model in Nimble were written by Steffen Oppel 

# load packages
library(tidyverse)
library(RODBC)
library(data.table)
library(lubridate)
library(reticulate)
library(splitstackshape)
library(MCMCvis)
library(nimble)
library(basicMCMCplots) # for trace plots called chainsPlot
library(parallel)
library(foreach)
library(doParallel)
library(dtplyr)
filter<-dplyr::filter
select<-dplyr::select
rename<- dplyr::rename

# for working local I have to set my working directory, I guess later this is done in the yaml script and can be removed
#setwd('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/')

################################################################################
#### LOAD DATA FROM DATA BASE AND FROM ESRIs SURVEY123 #########################
################################################################################

# get data from the MS Access data base - this Database has to be manually updated which might have been stopped in 2024
db <- odbcConnectAccess2007('Montserrat_Birds_2024.accdb') # connect to db, this might require Windows and a certain bit-version 
tblVisit <- sqlFetch(db, 'tblVisit') # observation level covariates
tblBirdData <- sqlFetch(db, 'tblBirdData') # bird count data
Point_hab <- sqlFetch(db, 'point_habitat_data')
Point_tree <- sqlFetch(db, 'point_tree_data')
tblSpecies <- sqlFetch(db, 'lkSpecies')
tblObserverAtSite <- sqlFetch(db, "tblObserverAtSite")
odbcClose(db) # close connection to data base

# check files
head(tblVisit)
head(tblBirdData)
head(tblSpecies)

# unzip the zip file in which ESRI's Survey123 App provides the data 
unzip(zipfile = 'Montserrat_Forest_Bird_Survey.zip')

# load data from Survey123 .csv files and remove all data which already exists in data base 
import <- fread("Montserrat_Forest_Bird_Survey_0.csv") # read in data from Survey123 App 
import <- import %>% filter(!(GlobalID %in% tblVisit$GlobalID)) %>%
  filter(Editor!="edwards.alice_Ext")
head(import)

################################################################################
#### manipulate VISIT data #####################################################
################################################################################

# tidy up bird data from database
tblVisit <- tblVisit %>% select(-Redundant, -NumObservers, -CensusID, -Cloud) %>% 
  mutate(Date = as.Date(Date))
head(tblVisit)

# manipulate from Survey123 imported data to match the data structure (tblVisit) from the data base AND remove all data before 2025
impVisit <- import %>% select(2:6,8,9,`Wind:`,`Observer name(s)/initials:`) %>%
  mutate(Date=as.Date(mdy_hms(`Date:`))) %>%
  mutate(day=day(Date), month=month(Date), year=year(Date), season=2, jday=yday(Date)) %>%
  mutate(VisitID=seq(max(tblVisit$VisitID)+1,max(tblVisit$VisitID)+nrow(import),1)) %>%
  mutate(Time = as.POSIXct(paste0(date(Date), ' ', `Start time:`), format = "%Y-%m-%d %H:%M")) %>% 
  rename(Point=`Point:`,Wind=`Wind:`,Observers=`Observer name(s)/initials:`, Route=`Route:`,Count= `Count:`,Rain=`Rain:`) %>%
  select(VisitID,GlobalID,Count,Date,Time,Rain,Wind,Observers,Point,Route,day,month,year,season,jday) %>%
  #filter(year == 2023) # this is just for DEMO purposes 
  filter(year > 2024) # filter for data after 2024 because all other data is already in the data base and any errors are corrected!
head(impVisit)

################################################################################
#### manipulate BIRD data ######################################################
################################################################################

# tidy up bird data from database
tblBirdData <- tblBirdData %>% select(VisitID, Species, Number)
head(tblBirdData)

# manipulate from Survey123 imported data to match the data structure (tblBirdData) from the data base // AND remove all data before 2025
impBirdData <- import %>% select(2,11:53) %>%
  gather(key=Species,value=Number,-GlobalID) %>%
  filter(!is.na(Number)) %>%
  mutate(VisitID=impVisit$VisitID[match(GlobalID,impVisit$GlobalID)]) %>%
  mutate(SpeciesCode=tblSpecies$SpeciesCode[match(tolower(Species),tolower(tblSpecies$Species))]) %>%
  select(VisitID,SpeciesCode,Number) %>%
  rename(Species=SpeciesCode) %>%
  filter(!is.na(VisitID)) 
head(impBirdData)


################################################################################
#### actually bring data together from imported tables from Survey123 and from the MS Access database
################################################################################

# combine bird data 
birds <- rbind(tblBirdData, impBirdData)
##### THINK ABOUT A PROPER CHECK THAT EVERYTHING WORKED OUT

# combine visit data 
str(impVisit$Time)
str(tblVisit$Time)
tblVisit_test <- rbind(impVisit, tblVisit)

# remove unneeded data from global environment (needed: birds, tblVisit, Point_hab, Point_tree, tblSpecies)
# rm(list = c('impBirdData', 'import', 'impVisit', 'tblObserverAtSite', 'tblBirdData', 'db'))
needed_objects <- c('birds', 'tblVisit', 'Point_hab', 'tblSpecies', 'Point_tree')
rm(list = setdiff(ls(), needed_objects)) # remove all unneeded objects 

################################################################################
### data preparation ###########################################################
################################################################################

# create year, choose species for analysis and select points that might cause isses for removing them 
YEAR <- year(Sys.Date()) # set the most recent year 
SPECIES <- c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA') # set species for which we actually want to conduct analysis 
removal<-c(99,76)	# these are the two Points that cause error in JAGS, but here 


#### 3: prepare siteCovs ####

# first manipulate/aggregate the tree data measured for 4 trees near each point by averaging to one value for each point
head(Point_tree)
tree <- Point_tree %>% rename(Point = Point_ID, treeheight = Tree_height) %>% 
  group_by(Point) %>% 
  mutate(DBH = mean(DBH), # take the mean for each point
         treeheight = mean(treeheight)) %>% 
  ungroup() %>% 
  select(Point, DBH, treeheight) %>% # remove all unneeded variables 
  distinct() # remove duplicate rows

# join the two siteCov data frames for easier handling and afterwards remove all unneeded variables/columns
siteCov <- Point_hab %>% 
  rename(Point = Point_ID) %>%
  left_join(tree, by="Point") %>%
  select(Point, Elevation, Location, DBH, treeheight,Canopy_cover) # select the habitat variables that are relevant

siteCov <- siteCov %>% 
  mutate(Point = as.numeric(Point), 
         Location = factor(Location)) # location as factor with 4 levels flat-open, midslope, ridge and valley

##### 3: prepare obsCov #### 

# remove point visits that are very obviously not needed (another survey protocol before 2011)
head(tblVisit)
tblVisit <- tblVisit %>% 
  filter(year %in% 2011:YEAR) %>% # remove all surveys which were done in another sampling scheme before 2011
  filter(Point %in% siteCov$Point) # remove all points that are not in the center hills and/or we don't have siteCovs
dim(tblVisit)



# tidy everything a bit up, remove unneeded columns and give the remaining a new order
head(tblVisit)
obsCov <- tblVisit %>% 
  select(-Route, -season, -GlobalID, -day, -month, -jday) %>% # remove unneeded columns
  select(VisitID, year, Point, Count, everything()) %>% # reorder columns 
  arrange(year, Point, Count)

# store the existing ones in the proper format
str(obsCov)
obsCov <- obsCov %>% 
  mutate(Point = as.numeric(Point),
         Rain = if_else(is.na(Rain), true = 0, false = 1), # assigns 0 for no rain and 1 for rain
         Wind = recode(Wind, 'calm' = 1, 'light' = 2, 'moderate' = 3, 'strong' = 4) %>% 
           factor(levels = c(1, 2, 3, 4), ordered = TRUE)) 

# calculate missing variables
obsCov <- obsCov %>% 
  mutate(march1st = as.POSIXct(paste0(year, '-03-01', format = '%Y-%m-%d')), # calculate difference to the first of March (therefore create a dummy variable)
         day = round(as.numeric(difftime(Date, march1st, units = 'days')), digits = 0), # calculate diff between dummy variable and Date
         time = as.numeric(difftime(Time, min(Time), units = 'mins'))) %>% # calculate time difference to the min time 
  select(-march1st) # delete dummy variable

# check if all surveys were done 
obsCov %>% group_by(year, Point) %>% summarise(number_counts = length(Count)) %>% filter(number_counts != 3) %>% ungroup()

# there are cases (23) where not all surveys were done in each year and point - include missing surveys by joining a df that includes all theoretically performed surveys to join with actual obsCov
# create df with all surveys that should have been done each year for each point
obsCov_theor <- data.frame(year = rep(2011:YEAR, each = length(siteCov$Point)), 
                           Point = rep(siteCov$Point, each = 3, times = length(2011:YEAR)), 
                           Count = rep(1:3, times = length(2011:YEAR)*length(siteCov$Point)))

# calculate the number of obsCovs / surveys that should theoretically be done 
length(siteCov$Point)*3*length(2011:YEAR) # this is the calculation where 2020 is still present (but no surveys were performed because of COVID)

# join theoretical df obsCov_theor with actual performed survey obsCov to get information of missing surveys
obsCov <- obsCov_theor %>% 
  left_join(obsCov, by = c('year', 'Point', 'Count'))


## find non-existing surveys and replace covariates with mean (irrelevant for modelling)
obsCov$Rain[is.na(obsCov$Rain)]<-0
obsCov$Wind[is.na(obsCov$Wind)]<-2   ## check the most common valuetable(obsCov$Wind)
obsCov$day[is.na(obsCov$day)]<-median(obsCov$day, na.rm=T)
obsCov$time[is.na(obsCov$time)]<-median(obsCov$time, na.rm=T)



#### 4: prepare bird countdata ####

##### digression to calculate activity covariate for obsCov BEFORE ALL OTHER SPECIES ARE ELIMINATED #####
ACT<-birds %>% 
  group_by(VisitID) %>% 
  summarise(activity = sum(Number,na.rm=T)) %>%
  filter(!is.na(VisitID))
dim(ACT)
ACT %>% filter(is.na(activity)) # this just checks if there are NAs, if 0 - everything is okay!

# first, check for duplicate entries in the data base bird count data for each species - here we will only care for all data between 2011 and the most recent year and SPECIES
# birds %>% 
#   filter(VisitID %in% obsCov$VisitID, Species %in% SPECIES) %>% # remove all observations from before 2011 which were done with another sampling scheme
#   group_by(VisitID, Species) %>%
#   summarise(n=length(Number)) %>%
#   filter(n>1) %>%
#   filter(!is.na(VisitID)) %>%
#   ungroup()
# okay, there are several duplicates to deal with - NEED TO BE SUMMED UP DUE TO HISTORIC DATA ENTRY (distance sampling)

# sum the numbers up, because a distance sampling framework was applied until 2011 and the duplicate entries refer to different distance bands, for all other cases: we cant solve it so we treat them the same 
## create a matrix with 0 counts
birdmat <- birds %>%
  filter(VisitID %in% obsCov$VisitID, Species %in% SPECIES) %>% # remove all observations from before 2011 which were done with another sampling scheme
  group_by(VisitID, Species) %>%
  summarise(Number = sum(Number)) %>% filter(!is.na(VisitID)) %>% 
  ungroup() %>%
  spread(key=Species, value=Number, fill=0)

# create df where for each species obsCov are prepared to fill in birddata using a join
# obsCov_spec <- bind_rows(replicate(length(SPECIES), obsCov, simplify = FALSE)) %>%
#   mutate(Species = rep(SPECIES, each = nrow(obsCov))) 
# head(obsCov_spec)

# bring both obsCov and bird data together 
# countdata <- birds %>% 
#   filter(VisitID %in% obsCov$VisitID, Species %in% SPECIES) %>% # remove all observations from before 2011 which were done with another sampling scheme
#   select(VisitID, Species, Number) %>% # select only the relevant columns
#   right_join(y = obsCov_spec, by = c('VisitID', 'Species')) %>% 
#   arrange(Species, year, Point, Count) %>%
#   mutate(Number = ifelse(is.na(Number), 0, Number), Number = ifelse(is.na(VisitID), NA, Number)) %>%  # fill in 0 for counts where no birds of the species were seen
#   select(VisitID, Species, year, Point, Count, Number)

countdata<-obsCov %>% 
  left_join(ACT, by="VisitID") %>%
  left_join(birdmat, by="VisitID")
dim(countdata)
dim(obsCov)


# calculate the theoretically number of observations that should have been done since 2011 at the selected points
length(siteCov$Point)*3*length(2011:YEAR)


#### 5: check prepared data ####

# siteCovs
head(siteCov)
dim(siteCov) # there are 88 points and 4 siteCovs as variables 

# obsCov
head(obsCov)
dim(obsCov) # there are number of points times number of replictes/counts each year per point times number of survey years observations - check (next line)
nrow(obsCov) == length(siteCov$Point)*3*length(2011:YEAR) # must be TRUE
obsCov %>% group_by(year, Point) %>% summarise(number_counts = length(Count)) %>% filter(number_counts != 3) %>% ungroup() %>% nrow() # must be 0

# countdata 
head(countdata)
nrow(countdata) == length(siteCov$Point)*3*length(2011:YEAR)# should be TRUE
unique(countdata$Species) # should contain all Species you want to study

countdata %>% filter(is.na(VisitID))

# fill in the objects you want to save from the environment for further analysis
needed_objects <- c('countdata', 'obsCov', 'siteCov', 'tblSpecies', 'YEAR', 'SPECIES')
rm(list = setdiff(ls(), needed_objects)) # remove all unneeded objects 

# save workspace to wd
#setwd() # set the file you want to store the data if it differs from the current wd
#save.image(paste0('data/MONTSERRAT_ANNUAL_DATA_INPUT', YEAR, '.RData'))
#save.image(paste0('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/Montserrat_git/MONTSERRAT_ANNUAL_DATA_INPUT', YEAR, '.RData'))

################################################################################
#### PREPARE DATA FOR ANALYSIS IN NIMBLE #######################################
################################################################################

SPECIES # print species codes 
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler", # save full species names in the same order as SPECIES 
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")

### define dimensions of arrays
nsites<-length(unique(siteCov$Point))
nyears<-length(unique(countdata$year))


### removed on 8 Oct 2024 because it did not explain anything
# ### summarise rainfall from Jan to March, productivity from PREVIOUS year will affect count in current year
# rain<-fread("data/MontserratRain2005_2023.csv",fill=TRUE) %>%
#   dplyr::filter(Variable=="RainMM") %>%
#   dplyr::filter(YEAR %in% seq(2010,2024,1)) %>%
#   dplyr::select(-Variable,-Total) %>%
#   gather(key="Month", value="mm",-YEAR) %>%
#   dplyr::filter(Month %in% c('JAN','FEB','MAR')) %>%
#   group_by(YEAR) %>%
#   summarise(rain=sum(mm)) %>%
#   mutate(rain=scale(rain)[,1])

###############################################################################
############## CREATE SITE COVARIATE DATA INPUT MATRIX   ######################
###############################################################################

siteCov<-siteCov %>% arrange(Point) %>%
  mutate(ridge=ifelse(Location=="Ridge",1,0)) %>%
  dplyr::select(Point,treeheight,Elevation,Canopy_cover,ridge) %>%
  mutate(tree=scale(treeheight)[,1], elev=scale(Elevation)[,1],canopy=scale(Canopy_cover)[,1])





###############################################################################
############## CREATE OBSERVATION COVARIATE DATA INPUT MATRIX   ###############
###############################################################################
SURVEYDATA<-countdata %>%
  arrange(Point,year,Count) %>%
  mutate(activity=ifelse(is.na(activity),mean(activity, na.rm=T),activity)) %>%
  mutate(time=scale(time),day=scale(day),activity=scale(activity))

### only needs standardisation if measured in mm, not as 0/1 variable
#meant<-mean(SURVEYDATA$rain, na.rm = TRUE)
#sdt<-sd(SURVEYDATA$rain, na.rm = TRUE)
#SURVEYDATA$rain<-(SURVEYDATA$rain-meant)/sdt


### create array for each covariate

wind<-array(NA, dim=c(nsites,3,nyears))
rain<-array(NA, dim=c(nsites,3,nyears))
time<-array(NA, dim=c(nsites,3,nyears))
day<-array(NA, dim=c(nsites,3,nyears))
ACT<-array(NA, dim=c(nsites,3,nyears))				## REPLACED ON 2 MAY WITH RAINFALL AMOUNT


### fill in array for each covariate
for (y in 2011:YEAR){
  obsC<-subset(SURVEYDATA, year==y)
  y<-match(y,c(2011:YEAR))						## translates the year (2011, 2012, etc.) into consecutive number (1,2,...) for array dimensions
  x<-obsC %>% dplyr::select(Point, Count, time) %>% tidyr::spread(key=Count, value=time) %>% dplyr::arrange(Point)
  time[,,y]<-as.matrix(x[,2:4])
  
  x<-obsC %>% dplyr::select(Point, Count, day) %>% tidyr::spread(key=Count, value=day) %>% dplyr::arrange(Point)
  day[,,y]<-as.matrix(x[,2:4])
  
  x<-obsC %>% dplyr::select(Point, Count, Wind) %>%
    mutate(Wind=ifelse(Wind<3,0,1)) %>%
    tidyr::spread(key=Count, value=Wind) %>% dplyr::arrange(Point)
  wind[,,y]<-as.matrix(x[,2:4])
  
  x<-obsC %>% dplyr::select(Point, Count, activity) %>% tidyr::spread(key=Count, value=activity) %>% dplyr::arrange(Point)
  ACT[,,y]<-as.matrix(x[,2:4])
  
  x<-obsC %>% dplyr::select(Point, Count, Rain) %>% tidyr::spread(key=Count, value=Rain) %>% dplyr::arrange(Point)
  rain[,,y]<-as.matrix(x[,2:4])
  
}





###############################################################################
####   REPLACE ALL NA IN COVARIATES otherwise "undefined node" error    #######
###############################################################################

for (d in 1:nyears){							### replace missing dates with mean for each survey round in each year
  ACT[is.na(ACT[,1,d]),1,d]<-mean(ACT[,1,d], na.rm=T)
  ACT[is.na(ACT[,2,d]),2,d]<-mean(ACT[,2,d], na.rm=T)
  ACT[is.na(ACT[,3,d]),3,d]<-mean(ACT[,3,d], na.rm=T)
  time[is.na(time[,1,d]),1,d]<-mean(time[,1,d], na.rm=T)
  time[is.na(time[,2,d]),2,d]<-mean(time[,2,d], na.rm=T)
  time[is.na(time[,3,d]),3,d]<-mean(time[,3,d], na.rm=T)
  day[is.na(day[,1,d]),1,d]<-mean(day[,1,d], na.rm=T)
  day[is.na(day[,2,d]),2,d]<-mean(day[,2,d], na.rm=T)
  day[is.na(day[,3,d]),3,d]<-mean(day[,3,d], na.rm=T)
  wind[is.na(wind[,1,d]),1,d]<-0
  wind[is.na(wind[,2,d]),2,d]<-0
  wind[is.na(wind[,3,d]),3,d]<-0
  rain[is.na(rain[,1,d]),1,d]<-0
  rain[is.na(rain[,2,d]),2,d]<-0
  rain[is.na(rain[,3,d]),3,d]<-0
}






######################################################################################
#############  WRITE THE NIMBLE MODEL AND SET INITS  ############################
######################################################################################

# Specify model in NIMBLE format
trend.model<-nimbleCode({
  
  ####  Priors ########
  loglam.int~dnorm(0,sd=1)          ##  mean abundance prior for intercept
  logitp.int~dnorm(0,sd=1)          ##  mean abundance prior for intercept
  trend~dnorm(0,sd=1)         ##  trend prior
  trend2~dnorm(0,sd=1)        ##  quadratic trend prior
  beta.elev~dnorm(0,sd=1)
  beta.canopy~dnorm(0,sd=1)
  beta.treeheight~dnorm(0,sd=1)
  bwind~dnorm(-1,sd=1)   ## wind can only have negative effect on detection
  brain~dnorm(-1,sd=1)   ## rain can only have negative effect on detection
  btime~dnorm(0,sd=1)
  b2time~dnorm(0,sd=1)
  bday~dnorm(0,sd=1)
  bridge~dnorm(0,sd=1)
  bact~dnorm(0,sd=1)
  
  ## SITE RANDOM EFFECT ##
  for(i in 1:nsite){
    lam.site[i]~dnorm(0,sd=sigma.site)    ## site-specific random effect
  }
  sigma.site~T(dt(0,1,2), 0, 10)
  
  ## YEAR RANDOM EFFECT FOR ANNUALLY VARYING DETECTION PROBABILITY ##
  for(year in 1:nyear){
    logit.p0[year]~dnorm(0,sd=sigma.year.p0)## detection probability
    lam.year[year]~dnorm(0,sd=sigma.year)    ## year-specific random effect with hierarchical centering from Kery email 5 June 2018
  }
  sigma.year~T(dt(0,1,2), 0, 10)
  sigma.year.p0~T(dt(0,1,2), 0, 10)
  
  ######### State and observation models ##############
  for(year in 1:nyear){
    for(i in 1:nsite){
      log(lambda[i,year])<- loglam.int+
        trend*primocc[year]+
        trend2*pow(primocc[year],2)+
        beta.elev*elev[i]+
        beta.treeheight*treeheight[i]+
        beta.canopy*canopy[i]+
        lam.site[i]+
        lam.year[year]
      N[i,year]~dpois(lambda[i,year]) ### T(dpois(lambda[i,year]),0,20) truncate N per point and year to a maximum of 20 to avoid ridiculously high values
      
      for(t in 1:nrep){
        M[i,t,year]~dbin(p[i,t,year],N[i,year])
        logit(p[i,t,year])<-logitp.int +
          logit.p0[year] +   ## random annual detection effect
          btime*time[i,t,year]+
          b2time * pow(time[i,t,year], 2) +
          bday*day[i,t,year]+
          bridge*ridge[i]+
          bwind*wind[i,t,year]+
          brain*rain[i,t,year]+
          bact*ACT[i,t,year]
        
      }
    }
    
    ### DERIVED PARAMETER FOR EACH YEAR ###
    totalN[year]<-sum(N[1:nsite,year])
    anndet[year]<-mean(p[1:nsite,1:nrep,year])
  }
  
  # Computation of fit statistic (Bayesian p-value)
  # Fit statistic for observed data
  # Also, generate replicate data and compute fit stats for them
  for(year in 1:nyear){
    for(i in 1:nsite){
      for(t in 1:nrep){
        
        # Actual data
        eval[i,t,year] <-N[i,year]*p[i,t,year] # Expected value
        sd.resi[i,t,year]<-sqrt(eval[i,t,year]*(1-p[i,t,year])) +0.5
        E[i,t,year]<-(M[i,t,year]-eval[i,t,year])/ sd.resi[i,t,year]
        E2[i,t,year] <- pow(E[i,t,year],2)
        
        # Replicate data sets
        M.new[i,t,year]~dbin(p[i,t,year],N[i,year])
        E.new[i,t,year]<-(M.new[i,t,year]-eval[i,t,year])/sd.resi[i,t,year]
        E2.new[i,t,year] <- pow(E.new[i,t,year], 2)
      }
    }
  }
  
  ### NIMBLE CANNOT SUM OVER 3 DIMENSIONS, hence the JAGS code does not work
  # fit <- sum(E2[1:nsite,1:nrep,1:nyear])# Sum up squared residuals for actual data set
  # fit.new <- sum(E2.new[1:nsite,1:nrep,1:nyear]) # Sum up for replicate data sets
  
  ## alternative solution from https://groups.google.com/g/nimble-users/c/fI8fXBpgIAE
  for(year in 1:nyear){
    for(i in 1:nsite){
      fsum[i,year]<- sum(E2[i,1:nrep,year])
      nsum[i,year] <- sum(E2.new[i,1:nrep,year])
    }
  }
  fit <- sum(fsum[1:nsite,1:nyear])
  fit.new <- sum(nsum[1:nsite,1:nyear])
  
}) ## end of nimble code chunk



######################################################################################################
########## CREATE INPUT DATA FOR NIMBLE - GENERIC PART     -----------------------
#######################################################################################################
#### data and inits that will be the same for all species
#### DISTINGUISH CONSTANTS AND DATA
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~

trend.constants <- list(nsite=nsites,
                        nrep=3,
                        #primocc=seq(2011:YEAR),
                        primocc=as.vector(scale(c(2011:YEAR))),
                        nyear=nyears,
                        elev=siteCov$elev,
                        treeheight=siteCov$tree,
                        canopy=siteCov$canopy,
                        rain=rain,
                        wind=wind,
                        day=day,
                        ridge=siteCov$ridge,
                        time=time,
                        ACT=ACT)





####   DEFINE INITIAL VALUES----     ################################
## MUST BE FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

inits.trend <- list(#N = Nst,
  trend=rnorm(1,0,1),
  trend2=rnorm(1,0,1),
  loglam.int = rnorm(1,0,1),
  logitp.int = rnorm(1,0,1),
  lam.site = rnorm(nsites,0,1),
  lam.year = rnorm(nyears,0,1),
  sigma.site = 1,
  sigma.year=1,
  logit.p0 = rnorm(nyears,0,1),
  sigma.year.p0=1,
  beta.canopy=rnorm(1,0,1),
  beta.treeheight=rnorm(1,0,1),
  beta.elev=rnorm(1,0,1),
  bwind=-1,
  brain=-1,
  bridge=-1,
  btime=-1,
  b2time=-1,
  bday=1,
  bact=2
)
#inits.trend$lam.site<-rnorm(nsites,0,inits.trend$sigma.site)
#inits.trend$lam.year<-rnorm(nyears,0,inits.trend$sigma.year)





####   DEFINE RUN SETTINGS AND OUTPUT DATA----     ################################

# Define parameters to be monitored
parameters.trend <- c("fit", "fit.new","trend","trend2","totalN","anndet")  #


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 15   #150000
n.burnin <- 10  #100000
n.chains <- 3


# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------

### fill in array for bird data and initial values
bird_s<-SURVEYDATA[,c(1,2,3,4,14)] %>%
  arrange(Point,year,Count) %>%
  rename(N=5) %>%
  mutate(N=if_else(is.na(VisitID),NA,N)) %>%  ### RE-INTRODUCE THE NAs for COUNTS THAT DID NOT TAKE PLACE #####
dplyr::select(Point,year,Count,N)
BIRD.y<-array(NA, dim=c(nsites,3,nyears))
for (y in 2011:YEAR){
  x<-bird_s %>%
    dplyr::filter(year==y) %>%
    dplyr::select(Point, Count, N) %>%
    tidyr::spread(key=Count, value=N) %>%
    dplyr::arrange(Point)
  yc<-match(y,c(2011:YEAR))						## translates the year (2011, 2012, etc.) into consecutive number (1,2,...) for array dimensions
  BIRD.y[,,yc]<-as.matrix(x[,2:4])
}
trend.data <- list(M = BIRD.y)
test <- nimbleModel(code = trend.model,
                    constants=trend.constants,
                    data = trend.data,
                    inits = inits.trend,
                    calculate=TRUE)

# USE TEST VALUES TO SUPPLEMENT INITS
# inits.trend$lp = array(rnorm(trend.constants$nsite*trend.constants$nrep*trend.constants$nyear, c(test$mu.lp), inits.trend$sigma.p),
#                        dim= c(trend.constants$nsite, trend.constants$nrep,trend.constants$nyear)) 
inits.trend$p = array(runif(trend.constants$nsite*trend.constants$nrep*trend.constants$nyear, 0.3,0.7),
                      dim= c(trend.constants$nsite, trend.constants$nrep,trend.constants$nyear)) 
test$calculate()
# inits.trend$p0

# test$initializeInfo()
# 
# # Missing values (NAs) or non-finite values were found in model variables:
# # M, p, lp, anndet,
# # E, E2, M.new, E.new, E2.new, fit, fit.new.
# 
# ### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# test$calculate()
# is.na(trend.constants) ## check whether there are NA in the data
# test$calculate(nodes="lam.site") # this seems to be ok
# test$logProb_lam.site
# test$logProb_lam.year
# test$logProb_p0
# log(test$logProb_p0/(1-test$logProb_p0))
# test$logProb_btime
# test$logProb_bday
# test$logProb_bridge
# test$logProb_bact
# test$logProb_M  # this is NA for all values
# test$logProb_lp  # this is NA for all values
# test$logProb_p  # this is NA for all values
# test$initializeInfo()
# #help(modelInitialization)
# 
# ### make sure that none of the logProbs result in NA or -Inf as the model will not converge
# configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits










######################################################################################
#############  START THE PARALLEL LOOP OVER EVERY SPECIES          ############################
######################################################################################
# trendout<-data.frame(species=SPECIES, timeframe=sprintf("2011-%i",YEAR), trend=0, lower95CI=-1, upper95CI=1, pval=0, slope=1)
# annestimates<-data.frame(species=rep(SPECIES, each=nyears), Year=seq(2011,YEAR), trend=0, lower95CI=-1, upper95CI=1, detprob=0, detproblower95CI=-1, detprobupper95CI=1) 

#setup parallel backend to use 8 processors - collapsed my machine
cl<-makeCluster(4)
registerDoParallel(cl)
Result <- foreach(s=SPECIES, .packages=c('nimble',"tidyverse","MCMCvis","tidyverse","dplyr","data.table")) %dopar% {		#.combine = rbind,
  
  # for (s in SPECIES){
  # s="MTOR"
  
  ######################################################################################
  #############  TAKE SUBSET OF DATA FOR FOCAL SPECIES AND SORT THE TABLES    ###################
  ######################################################################################
  
  bird_s<-SURVEYDATA[,c(1,2,3,4,match(s,colnames(SURVEYDATA)))] %>%
    arrange(Point,year,Count) %>%
    rename(N=5) %>%
    #mutate(N=if_else(is.na(VisitID),NA,N)) %>%  ### RE-INTRODUCE THE NAs for COUNTS THAT DID NOT TAKE PLACE #####
  #mutate(N=if_else(N>12,12,N)) %>%  ### truncate data to no more than 12 individuals (3 per quadrant) to avoid big outliers
  dplyr::select(Point,year,Count,N)
  
  
  
  ###############################################################################
  ############## CREATE BIRD DATA INPUT MATRIX   ################################
  ###############################################################################
  #### FILL THE MISSING DATA WITH MEAN VALUES FOR INITS
  ## https://groups.google.com/g/nimble-users/c/wCwacQPLR2w?pli=1
  
  ### create array to be filled with data
  BIRD.y<-array(NA, dim=c(nsites,3,nyears))
  inits.y<-array(NA, dim=c(nsites,3,nyears))
  inits.new<-array(NA, dim=c(nsites,3,nyears))
  
  ### fill in array for bird data and initial values
  for (y in 2011:YEAR){
    x<-bird_s %>%
      dplyr::filter(year==y) %>%
      dplyr::select(Point, Count, N) %>%
      tidyr::spread(key=Count, value=N) %>%
      dplyr::arrange(Point)
    yc<-match(y,c(2011:YEAR))						## translates the year (2011, 2012, etc.) into consecutive number (1,2,...) for array dimensions
    BIRD.y[,,yc]<-as.matrix(x[,2:4])
    
    x<-bird_s %>%
      mutate(N=ifelse(is.na(N),median(bird_s$N, na.rm=T),N)) %>%   ### fill in missing values
      dplyr::filter(year==y) %>%
      dplyr::select(Point, Count, N) %>%
      tidyr::spread(key=Count, value=N) %>%
      dplyr::arrange(Point)
    inits.y[,,yc]<-as.matrix(x[,2:4])
    inits.new[,,yc]<-as.matrix(x[,2:4]) 
  }
  
  #### GET THE MAXIMUM COUNT PER POINT PER YEAR FOR INITIAL VALUES
  Nst<-as.matrix(bird_s %>%
                   mutate(N=ifelse(is.na(N),median(bird_s$N, na.rm=T),N)) %>%   ### fill in missing values - switch to max if there is invalid parent error
                   group_by(Point, year) %>%
                   summarise(K=max(N, na.rm=T)) %>%
                   spread(key=year,value=K, fill=max(bird_s$N,na.rm=T)) %>%
                   ungroup() %>%
                   arrange(Point) %>%
                   dplyr::select(-Point))
  
  
  
  ######################################################################################################
  ########## CREATE INPUT DATA FOR NIMBLE ------------------------
  #######################################################################################################
  
  #### DISTINGUISH CONSTANTS AND DATA
  # Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
  # Data are values that you might want to change, basically anything that only appears on the left of a ~
  
  trend.data <- list(M = BIRD.y)
  
  
  ####   ADD INITIAL VALUES----     ################################
  ## MUST ADD Nst TO INITIAL VALUESBE FOR ALL PARAMETERS
  ## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8
  # inits.trend$lp = array(rnorm(trend.constants$nsite*trend.constants$nrep*trend.constants$nyear, c(test$mu.lp), inits.trend$sigma.p),
  #                        dim= c(trend.constants$nsite, trend.constants$nrep,trend.constants$nyear))
  inits.trend$p = array(runif(trend.constants$nsite*trend.constants$nrep*trend.constants$nyear, 0.3,0.7),
                        dim= c(trend.constants$nsite, trend.constants$nrep,trend.constants$nyear)) 
  inits.trend$N = Nst
  inits.trend$M = inits.y
  inits.trend$M.new = inits.new
  
  #allchaininits.trend <- list(inits.trend, inits.trend, inits.trend)
  
  
  
  
  ###############################################################################
  ####   RUN THE MODEL IN NIMBLE  --------------------###########################
  ###############################################################################
  
  
  ### this takes 3-5 hrs for 250000 iterations and converges for most species
  TRENDMOD <- nimbleMCMC(code = trend.model,
                         constants=trend.constants,
                         data = trend.data,
                         inits = inits.trend,
                         monitors = parameters.trend,
                         thin=4,
                         niter = n.iter,
                         nburnin = n.burnin,
                         nchains = n.chains,
                         progressBar = getNimbleOption("MCMCprogressBar"),
                         summary=T)
  
  
  saveRDS(TRENDMOD,sprintf("output/%s_trend_model_nimble.rds",s))
  
  
  
  ###############################################################################
  ####   EVALUATE MODEL FIT WITH BAYESIAN P VALUE   #############################
  ###############################################################################
  ### COMBINE SAMPLES ACROSS CHAINS
  MCMCout<-as_tibble(rbind(TRENDMOD$samples[[1]],TRENDMOD$samples[[2]],TRENDMOD$samples[[3]]))
  
  ### PLOT AND SAVE GoF PLOT
  ylow<-round((min(MCMCout$fit,MCMCout$fit.new, na.rm=T)-50)/1000,1)*1000
  yup<-round((max(MCMCout$fit,MCMCout$fit.new, na.rm=T)+50)/1000,1)*1000
  pdf(sprintf("output/%s_trendmodel_fit2024.pdf",s), width=10, height=10, title="")
  plot(MCMCout$fit, MCMCout$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(ylow, yup), ylim = c(ylow, yup))
  abline(0, 1, lwd = 2, col = "black")
  dev.off()
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # EXAMINE OUTPUT AND DIAGNOSTICS WITH MCMCvis
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  out<- as.data.frame(MCMCsummary(TRENDMOD$samples, params=c("trend","trend2","totalN","anndet")))
  out$parameter<-row.names(out)
  out$species<-s
  out$BayesP<-mean(MCMCout$fit > MCMCout$fit.new)
  out$GoFSlope<-mean(MCMCout$fit) / mean(MCMCout$fit.new)
  names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
  fwrite(out,sprintf("output/%s_trend_estimates2024.csv",s))
  
  
  # MCMCplot(TRENDMOD$samples, params=c("trend","totalN","anndet"))
  # ggsave(sprintf("output/%s_trend_estimates.jpg"), height=9, width=9)
  # 
  # ## look at the chains and whether they mixed well
  # chainsPlot(TRENDMOD$samples,var=c("trend"))
  # chainsPlot(TRENDMOD$samples,var=c("totalN"))
  # chainsPlot(TRENDMOD$samples,var=c("anndet"))
  
  
  
  
  
  
  # 
  # ###############################################################################
  # ####   PRINT AND SAVE MODEL OUTPUT                #############################
  # ###############################################################################
  # 
  # # Summarize posteriors
  # #print(model, dig = 3)
  # 
  # write.table(model$summary,sprintf("%s_abund_estimates2023_p%f.csv",s,pval), sep=",")
  # 
  # 
  # trendout[trendout$species==s,3]<-round(out[1,4],3)
  # trendout[trendout$species==s,4]<-round(out[1,3],3)
  # trendout[trendout$species==s,5]<-round(out[1,5],3)
  # # trendout[trendout$species==s,6]<-pval
  # # trendout[trendout$species==s,7]<-slope
  # 
  # annestimates[annestimates$species==s,3]<-round(model$summary[2:(nyears+1),5],3)
  # annestimates[annestimates$species==s,4]<-round(model$summary[2:(nyears+1),3],3)
  # annestimates[annestimates$species==s,5]<-round(model$summary[2:(nyears+1),7],3)
  # annestimates[annestimates$species==s,6]<-round(model$summary[(nyears+4):(dim(model$summary)[1]-1),5],3)
  # annestimates[annestimates$species==s,7]<-round(model$summary[(nyears+4):(dim(model$summary)[1]-1),3],3)
  # annestimates[annestimates$species==s,8]<-round(model$summary[(nyears+4):(dim(model$summary)[1]-1),7],3)
  # 
  # ###############################################################################
  # ####   CREATE TREND PLOT AND SAVE AS PDF          #############################
  # ###############################################################################
  # 
  # trendlabel<- paste("Trend: ",trendout[trendout$species==s,3]," (",trendout[trendout$species==s,4]," - ",trendout[trendout$species==s,5],")", sep="")
  # ggplot(annestimates[annestimates$species==s,], aes(x=Year,y=trend)) +
  #   geom_line(colour="indianred", linewidth=1.5) +
  #   geom_ribbon(aes(ymin = lower95CI, ymax = upper95CI), fill="indianred", alpha = 0.2) +
  # 
  #   ## format axis ticks
  #   scale_x_continuous(name="Year", breaks=seq(2011,2023,2), labels=as.character(seq(2011,2023,2)))+
  #   #scale_y_continuous(name="Number of Birds at 67 Sampling Points", breaks=seq(0,4000,500), labels=as.character(seq(0,4000,500)))+
  #   ylab(sprintf("Number of %s at %i sampling points",s,nsites)) +
  # 
  #   annotate("text", x = -Inf, y = Inf, label = trendlabel, vjust = 2, hjust = -0.1,size=6, color="black") +
  # 
  #   ## beautification of the axes
  #   theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         axis.text=element_text(size=18, color="black"),
  #         axis.title=element_text(size=18),
  #         strip.text.x=element_text(size=18, color="black"),
  #         axis.title.y=element_text(margin=margin(0,20,0,0)),
  #         strip.background=element_rect(fill="white", colour="black"))
  # 
  # ggsave(sprintf("MONTSERRAT_%s_abund_plot2023.pdf",s), width=12, height=9)
  # 
  # 
  # 
  
  ##### END THE LOOP ACROSS SPECIES ####
}



############################################################################################################
####   COLLATE DATA FROM INDIVIDUAL MODEL OUTPUT FILES          #############################
############################################################################################################
#setwd("C:\\STEFFEN\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat\\output")
# setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat\\output")
# setwd('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/output') # for Filibert
#setwd("C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat\\Montserrat\\output")
nsites<-length(unique(siteCov$Point))
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")
allout<-list.files(pattern="trend_estimates2024.csv")
annestimates<-tibble()
trendout<-tibble()
for (f in allout){
  x<-fread(f)
  trendout<-trendout %>% bind_rows(x %>% filter(parameter %in% c("trend","trend2")))
  annestimates<-annestimates %>% bind_rows(x %>% filter(substr(parameter,1,6)=="totalN"))
}


write.table(annestimates,"Annual_estimates2024.csv", row.names=F, sep=",")
write.table(trendout,"Trend_estimates2024.csv", row.names=F, sep=",")




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




