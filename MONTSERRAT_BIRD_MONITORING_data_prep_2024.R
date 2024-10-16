######################################################################################
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  DATA PREPARATION FOR MULTI-YEAR ANALYSIS   ############################
#############  filibert.heim@posteo.de and steffen.oppel@gmail.com      ##############
######################################################################################

##### Montserrat forest bird surveys data preparation ####
##### written in September 2024 by Filibert Heim  #### 

##### 1: load packages and make other preparations ---------- ####

# load packages
library(RODBC)
library(tidyverse)
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

# fill in the most recent year and set the species you want to analyse
YEAR <- 2024
SPECIES<-c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA')
removal<-c(99,76)																			## removes only two points that cause error in JAGS

# set working directory
getwd()
setwd('C:/Users/filib/Documents/Praktika/Sempach/Montserrat') # for Filibert
setwd('C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/UKOT/Montserrat/Analysis/Population_status_assessment/AnnualMonitoring/Montserrat') # for Steffen



#### 2: load data from ms access database ####

# connect with database to load data from query and afterwards close connection again 
db <- odbcConnectAccess2007('data/Montserrat_Birds_2024.accdb') # change name of the db to the actual one 
tblVisit <- sqlFetch(db, 'tblVisit') # observation level covariates
birds <- sqlFetch(db, 'tblBirdData') # bird count data
Point_hab <- sqlFetch(db, 'point_habitat_data')
Point_tree <- sqlFetch(db, 'point_tree_data')
species <- sqlFetch(db, 'lkSpecies')
odbcClose(db)

save.image("data/Montserrat_Birds_2024_accdb.RData")
load("data/Montserrat_Birds_2024_accdb.RData")

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
  select(-Redundant, -Cloud, -Observers, -NumObservers, -Route, -season, -GlobalID, -CensusID, -day, -month, -jday) %>% # remove unneeded columns
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
ACT %>% filter(is.na(activity))

# first, check for duplicate entries in the data base bird count data for each species - here we will only care for all data between 2011 and the most recent year and SPECIES
birds %>% 
  filter(VisitID %in% obsCov$VisitID, Species %in% SPECIES) %>% # remove all observations from before 2011 which were done with another sampling scheme
  group_by(VisitID, Species) %>%
  summarise(n=length(Number)) %>%
  filter(n>1) %>%
  filter(!is.na(VisitID)) %>%
  ungroup()
# okay, there are several duplicates to deal with - NEED TO BE SUMMED UP DUE TO HISTORIC DATA ENTRY (distance sampling)

# sum the numbers up, because a distance sampling framework was applied until 2011 and the duplicate entries refer to different distance bands, for all other cases: we cant solve it so we treat them the same 
birds <- birds %>%
  filter(VisitID %in% obsCov$VisitID, Species %in% SPECIES) %>% # remove all observations from before 2011 which were done with another sampling scheme
  group_by(VisitID, Species) %>%
  summarise(Number = sum(Number)) %>%
  filter(!is.na(VisitID)) %>% ungroup()

## create a matrix with 0 counts
birdmat<-birds %>% spread(key=Species, value=Number, fill=0)

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

#### 6: save workspace for analysis #### 

# fill in the objects you want to save from the environment for further analysis
needed_objects <- c('countdata', 'obsCov', 'siteCov', 'species', 'YEAR', 'SPECIES')
rm(list = setdiff(ls(), needed_objects)) # remove all unneeded objects 

# save workspace to wd
#setwd() # set the file you want to store the data if it differs from the current wd
save.image(paste0('data/MONTSERRAT_ANNUAL_DATA_INPUT', YEAR, '.RData'))
save.image(paste0('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/Montserrat_git/MONTSERRAT_ANNUAL_DATA_INPUT', YEAR, '.RData'))


