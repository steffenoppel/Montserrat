#### Script for data preparation for yearly analysis - brings together data from database and new data from SURVEY123 #####
#### this script is part of the automated workflow for an annual report for the Centre Hills Forest Bird Monitoring in Montserrat ####
#### Script written by Filibert Heim, filibert.heim@posteo.de, in Nov 2024 with important parts of the script token and adapted from Steffen Oppel 

# at point 6. Prepare siteCovs you can include additional environmental/habitat variables as siteCovs!
# yearly data preparation workflow for all sort of analysis is implemented until line 250 where this data is exported (the further code is preparation for yearly N-mix models in NIMBLE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
# load packages
library(tidyverse) 
library(data.table)
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

# create year, choose species for analysis and select points that might cause isses for removing them 
# YEAR <- year(Sys.Date()) # set the most recent year 
SPECIES <- c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA') # set species for which we actually want to conduct analysis 
removal<-c(99,76)	# these are the two Points that cause error in JAGS, but here 

# set working directory
# getwd()
# setwd('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/') # this works when applying the script local but not on git actions - there data/ file has to be used  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Load data from database (until 2024) and Survey123 upload (after 2024)  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data from Montserrat Access data base which has been saved in a RData file by using the script 'Montserrat_forest_bird_monitoring_data_extraction_from_db.R'
load(file = 'data/Montserrat_forest_bird_monitoring_data_2011_2024.RData')

# load data from Survey123 .csv files and remove all data which already exists in data base
unzip(zipfile = 'data/Montserrat_Forest_Bird_Survey.zip', exdir = 'data/') # unzip the zip file in which ESRI's Survey123 App provides the data 
import <- fread("data/Montserrat_Forest_Bird_Survey_0.csv") # read in data from Survey123 App 
# import <- import %>% filter(!(GlobalID %in% tblVisit$GlobalID)) %>% # I'm not totally sure if this step is needed, revise this later 
  # filter(Editor!="edwards.alice_Ext")
head(import)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Manipulate data from the Survey123 app --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.1. TROUBLESHOOT AND WEED OUT ERRORS WITH POINT LABEL--------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### FIND POINTS WITH ><3 counts

impVisit %>% group_by(year,Point,Count) %>%
  summarise(n=length(unique(Time))) %>%
  filter(n>1)


### FIND POINTS WITH ><3 counts

CheckErrors<-impVisit %>% group_by(year,Route,Point) %>%
  summarise(n=length(unique(Time))) %>%
  filter(n!=3) %>%
  mutate(Problem=ifelse(n>3,"TooMany","TooFew")) %>%
  arrange(desc(Problem))

for(l in 1:dim(CheckErrors[CheckErrors$Problem=="TooMany",])[1]){
  
  ## go through case by case and find points with >3 visits
  case<-impVisit %>% filter(year==CheckErrors$year[l],Point==CheckErrors$Point[l])
  
  ## identify the duplicate
  duplicates<-case %>% group_by(year,Route,Point,Date) %>%
    summarise(n=length(VisitID)) %>%
    filter(n>1)
  
  ## duplicate VisitIDs
  dupIDs<-case %>% filter(Date==duplicates$Date, Point==duplicates$Point)
  
  ## check whether there is a missing count on the same route and day
  pot.missing<-CheckErrors %>% filter(Problem=="TooFew") %>%
    filter(Route==duplicates$Route)
  
  ## check whether the missing count is usually before or after the duplicate count
  check_sequence<-impVisit %>% filter(year==CheckErrors$year[l],Point %in% c(duplicates$Point,pot.missing$Point)) %>%
    arrange(Count,Time) %>%
    group_by(Count) %>%
    mutate(SEQ=seq_along(Time)) %>%
    ungroup() %>%
    filter(Point %in% c(pot.missing$Point)) %>%
    summarise(seq=min(SEQ))
      
  
  ## CHANGE THE POINT NUMBER IN THE VISIT DATA
  
  impVisit$Point[impVisit$VisitID %in% dupIDs$VisitID][check_sequence$seq]<-pot.missing$Point
  
}



### TEST WHETHER IT HAS WORKED

impVisit %>% group_by(year,Route,Point) %>%
  summarise(n=length(unique(Time))) %>%
  filter(n!=3) %>%
  mutate(Problem=ifelse(n>3,"TooMany","TooFew")) %>%
  arrange(desc(Problem))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.2. TROUBLESHOOT AND WEED OUT ERRORS WITH COUNT NUMBER--------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### FIND POINTS WITH incomplete set of 'counts'

CheckErrors<-impVisit %>% group_by(year,Route,Point) %>%
  summarise(n=length(unique(Count))) %>%
  filter(n<3)


for(l in 1:dim(CheckErrors)[1]){
  
  ## go through case by case and create sequential count
  case<-impVisit %>% filter(year==CheckErrors$year[l],Point==CheckErrors$Point[l]) %>%
    arrange(Time) %>%
    mutate(NEWCount=seq_along(Time))
  
  ## CHANGE THE COUNT NUMBER IN THE VISIT DATA
    impVisit$Count[match(case$VisitID,impVisit$VisitID)]<-case$NEWCount
  
}



### TEST WHETHER IT HAS WORKED
impVisit %>% group_by(year,Route,Point) %>%
  summarise(n=length(unique(Count))) %>%
  filter(n<3)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Combine new data from Survey123 and data until 2024 --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# combine bird data 
birds <- rbind(tblBirdData, impBirdData)

# combine visit data 
tblVisit <- rbind(impVisit, tblVisit)

# get rid of object that arent needed anymore 
needed_objects <- c('birds', 'tblVisit', 'Point_hab', 'tblSpecies', 'Point_tree', 'select', 'filter', 'rename', 'YEAR', 'removal', 'SPECIES')
rm(list = setdiff(ls(), needed_objects)) # remove all unneeded objects 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Overview of the data --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get an overview of what is set 
# YEAR # the most recent year data is available for 
SPECIES # all species for which we want to prepare the data 
removal # there are two points which might cause errors in JAGS, however, those are still included 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Prepare siteCovs --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  select(Point, Elevation, Location, DBH, treeheight, Canopy_cover) # select the habitat variables that are relevant

siteCov <- siteCov %>% 
  mutate(Point = as.numeric(Point), 
         Location = factor(Location)) # location as factor with 4 levels flat-open, midslope, ridge and valley

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Prepare obsCovs --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# tidy everything a bit up, remove unneeded columns and give the remaining a new order
head(tblVisit)
obsCov <- tblVisit %>% 
  select(-Route, -season, -GlobalID, -day, -month, -jday) %>% # remove unneeded columns
  select(VisitID, year, Point, Count, everything()) %>% # reorder columns 
  arrange(year, Point, Count)
head(obsCov)

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

# create variable YEAR as most recent year where data exists!
YEAR <- max(obsCov$year) # there should not be any NAs

# there are cases (23) where not all surveys were done in each year and point - include missing surveys by joining a df that includes all theoretically performed surveys to join with actual obsCov
# create df with all surveys that should have been done each year for each point
#obsCov_theor <- data.frame(year = rep(2011:YEAR, each = length(siteCov$Point)), 
#                           Point = rep(siteCov$Point, each = 3, times = length(2011:YEAR)), 
#                           Count = rep(1:3, times = length(2011:YEAR)*length(siteCov$Point)))

obsCov_theor <- expand.grid(year = c(2011:YEAR), 
                           Point = siteCov$Point, 
                           Count = c(1:3))

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Prepare bird countdata --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# digression to calculate activity covariate for obsCov BEFORE ALL OTHER SPECIES ARE ELIMINATED #####
ACT <- birds %>% 
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

countdata <- obsCov %>% 
  left_join(ACT, by="VisitID") %>%
  left_join(birdmat, by="VisitID")
dim(countdata)
dim(obsCov)


# calculate the theoretically number of observations that should have been done since 2011 at the selected points
length(siteCov$Point)*3*length(2011:YEAR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8. Check the data and export prepared data as RData  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
needed_objects <- c('countdata', 'obsCov', 'siteCov', 'tblSpecies', 'YEAR', 'SPECIES', 'removal')
rm(list = setdiff(ls(), needed_objects)) # remove all unneeded objects 

# save the prepared data for other possible analysis 
save.image('data/MONTSERRAT_ANNUAL_DATA_INPUT.RData')
# save.image(paste0('data/MONTSERRAT_ANNUAL_DATA_INPUT', YEAR, '.RData')) # for storing the file with a year flag



