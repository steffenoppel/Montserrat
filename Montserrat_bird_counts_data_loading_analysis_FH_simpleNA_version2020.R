
##### Montserrat Forest Bird Counts - data import and analysis with 'unmarked' #####
##### written in June 2023 by  Filibert Heim, filibert.heim@posteo.de          #####

##### 1: load required packages ####

#library(RODBC)
library(reshape)
library(tidyverse)
library(lubridate)
library(unmarked)
library(AICcmodavg) # package for model selection and goodness-of-fit tests
library(MuMIn)
filter<-dplyr::filter
select<-dplyr::select
rename<-dplyr::rename

##### 2: set YEAR and SPECIES that should be analysed ####
YEAR <- 2023 # fill in the most recent year that data is available for
SPECIES <- c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH',
             'SNPI','CAEL','BANA') # fill in species the analysis should be made for

###### 2.1: load the prepared data and make last preparation #####

# load prepared data from Steffen 
setwd('C:/Users/sop/Documents/Steffen/RSPB/Montserrat')
# load(file = 'MONTSERRAT_ANNUAL_DATA_INPUT2023.RData')
# ls()
# 
# # load rain and wind data from the database for obsCovs and yearlySiteCovs
# setwd('C:/Users/sop/Documents/Steffen/RSPB/Montserrat')
# mt_db <- odbcConnectAccess2007('Montserrat_Birds_2023.accdb')
# rain_y <- sqlQuery(mt_db, 'SELECT * FROM yearlyRainCov')
# siteCov_db <- sqlQuery(mt_db, 'SELECT * FROM SiteCov')
# obsCov_db <- sqlQuery(mt_db, 'SELECT * FROM obsCovariates_simple')
# obsCov_db$Point <- as.character(obsCov_db$Point)
# obsCov_db$Wind <- as.factor(obsCov_db$Wind)
# odbcClose(mt_db)
# save.image('MONTSERRAT_ANNUAL_DATA_INPUT2023.RData')
load(file = 'MONTSERRAT_ANNUAL_DATA_INPUT2023.RData')
filter<-dplyr::filter
select<-dplyr::select
rename<-dplyr::rename
ls()

###### 2.2: check the data and remove unneeded stuff #####

head(activity) # (tidyverse-)data.frame with calculated bird-activity data
head(add) # I don't know what kind of data there is inside
head(birds) # data like its provided by the database
head(birds2023) # bird data from 2023
head(COUNTDATA) # completely prepared bird observation data for analysis!
head(MTOR_count) # RODBC connection to data base
head(nyears) # number of years with bird data from the monitoring 
head(obsCov) # prepared obsCov data, rename columns and remove 'Rain', 'Year' (2nd) - but also without wind data!
head(pointmax) # (tidyverse-)data.frame with extracted maximal activity per point
head(siteCov) # prepared siteCov data for analysis, but there are mistakes in column names!
head(siteCov_db) # loaded siteCovs from database to check colnames of siteCov
head(species) # species list and abundance over all surveys
bird_names <- species
rm(species)
print(SPECIES)
head(startdate) # first date in every year to calculate numeric dates
head(summary2023) # hohle abundance over all visits in 2023
head(SURVEYDATA) # nearly the same data.frame like obsCOV
head(surveys2023) # obsCov from all visits in 2023
head(table1) # I don't really know what kind of data this data.frame contains
head(totals2023) # I don't really know what kind of data this data.frame contains
print(y) # last recent year 
head(rain_y)

###### 2.3: relable variate names/columns #####

obsCov <- select(obsCov, -Rain, -Year)
names(obsCov)[1:10] <- c('year','point','count','obs','skill','rain','wind','day',
                         'time','activity')

names(siteCov)[c(1:9, 21)] <- c('point', 'eastings', 'northings', 'habitat', 'dbh', 'distance', 
                         'treeheight', 'elevation', 'slope', 'alt') # go on with renaming

names(COUNTDATA) <- c('species', 'year', 'point', 'count', 'N')

###### 2.4: digression to calculate missing wind values as obsCov ####

obsCov_db <- obsCov_db %>% 
  select(-Time, -Rain, -MinOfskill, -Date, -CountOfObservername) %>%
  rename(point = Point, count = Count)
obsCov_db$Wind <- ifelse(obsCov_db$Wind %in% c('light','calm'),'low','high')
obsCov <- obsCov %>% 
  left_join(x = obsCov, y = obsCov_db, by = c('year', 'point', 'count')) %>%
  select(-wind) %>% 
  rename(wind = Wind) %>%
  mutate(wind = as.character(wind)) %>% 
  replace(is.na(.), 'low') %>%
  mutate(wind = as.factor(wind))

###### 2.5: remove unneeded object from global environment

rm(list = c('surveys2023', 'pointmax', 'table1', 'totals2023', 'birds', 
            'birds2023', 'activity', 'summary2023', 'nyears', 'MTOR_count', 
            'SURVEYDATA', 'rain_y', 'add'))

##### 3: provide an unmarkedMultFrame for multi-season analysis with 'colext()' ####

numPrimary <- length(unique(COUNTDATA$year)) # calculates number of years

occdataColext <- COUNTDATA %>% filter(species == 'MTOR') %>% # choose species
  left_join(obsCov %>% mutate(point = as.numeric(point)), by = c('point', 'count', 'year')) %>%
  mutate(N = ifelse(obs == 'NA', 'NA', N)) %>%
  select(-obs, -skill, -rain, -day, -time, -activity, -wind) %>%
  mutate(occupancy = ifelse(N > 0, 1, 0)) %>% # convert in detection/non-detection
  mutate(season = paste(year, count, sep = '_')) %>% # connect year, count to string
  select(-species, -N, -year, -count) %>% # remove unneeded columns 
  spread(key = season, value = occupancy) %>% # spread a key-value pair across multiple columns
  arrange(point) # sort rows in order of point number

siteCovsColext <- siteCov %>% mutate(point = as.numeric(point)) %>% 
  arrange(point) # sort rows in order of point number

obsCovsColext <- COUNTDATA %>% filter(species == 'MTOR') %>%
  left_join((obsCov %>% mutate(point = as.numeric(point))), by = c('point', 'count', 'year')) %>% 
  select(-N, -species, -skill) %>% # remove columns
  arrange(point, year, count) # sort rows in order of year, point, count
  
# mutate(year = ifelse(obs == 'NA', 'NA', year), # fill NA in fields from rows where no surveys were made
#      point = ifelse(obs == 'NA', 'NA', point), 
#       count = ifelse(obs == 'NA', 'NA', count))

###### 3.2: check dimensions ####

dim(occdataColext)
dim(siteCovsColext)
dim(obsCovsColext)/(3*numPrimary)

###### 3.3: create unmarkedMultFrame ####

umf <- unmarkedMultFrame(y = occdataColext[,2:40], siteCovs = siteCovsColext, 
                         obsCovs = obsCovsColext, numPrimary = numPrimary)
summary(umf) # looks good
str(umf) # looks good, but the values have to be scaled

###### 3.4: scale numeric variables to fitting problems ####

siteCovs(umf)[c(2,3,5:17,19,21:38)] <- scale(siteCovs(umf)[c(2,3,5:17,19,21:38)])
obsCovs(umf)[c(5:7)] <- scale(obsCovs(umf)[c(5:7)])
str(umf)

##### 4: multi-season, dynamic one species occupancy model ####

###### 4.1: fit models for detection probabilty p() first for modSel ####

# fit models for detection probabilty p() manually, go on with fitList(), modSel()
fm1 <- colext(~1, ~1, ~1, ~1, data = umf, se =T)
fm2 <- colext(~1, ~1, ~1, ~day, data = umf, se =T)
fm3 <- colext(~1, ~1, ~1, ~day+time, data = umf, se =T)
fm4 <- colext(~1, ~1, ~1, ~day+time+rain, data = umf, se =T)
fm5 <- colext(~1, ~1, ~1, ~day+time+rain+wind, data = umf, se =T)
fm6 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se =T)
fm7 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se =T)

# organise models in fitList before printing the modSel-table
MTOR_p_fitList <- fitList('psi(.)g(.)e(.)p(.)' = fm1, 
                         'psi(.)g(.)e(.)p(~day)' = fm2, 
                         'psi(.)g(.)e(.)p(~day+time)' = fm3, 
                         'psi(.)g(.)e(.)p(~day+time+rain)' = fm4,
                         'psi(.)g(.)e(.)p(~day+time+rain+wind)' = fm5,
                         'psi(.)g(.)e(.)p(~day+time+rain+wind+activity)' = fm6,
                         'psi(.)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm7)
(MTOR_p_modSel <- modSel(MTOR_p_fitList))
# best model: psi(.)g(.)e(.)p(~day+time+rain+wind+activity+ridge), go on with this 

# fit global model for detection probability p(), modSel with dredge(), get.models()

# small global model, covariates only additively
MTOR_p_global_model <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se =T)
# MTOR_p_global_model_dredge <- dredge(MTOR_p_global_model, trace = 2)
# print(MTOR_p_global_model_dredge)
# get.models(MTOR_p_global_model_dredge, subset = delta <= 2)
# # best model: psi(.)g(.)e(.)p(~activity+day+rain+ridge+time+1) - same like with modSel
# 
# # big global model, compute all interactions that are possible
# MTOR_p_global_model_full <- colext(~1, ~1, ~1, ~day*time+rain*wind*activity*ridge, data = umf, se =T)
# MTOR_p_global_model_full_dredge <- dredge(MTOR_p_global_model, trace = 2)
# print(MTOR_p_global_model_full_dredge)
# get.models(MTOR_p_global_model_full_dredge, subset = delta <= 2)
# # best model: psi(.)g(.)e(.)p(~activity+day+rain+ridge+time+1) - same like with modSel
# 
# ###### 4.2: fit models for initial occupancy psi() first for modSel ####
# 
# fm8 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T) # intercept only model
# # add alt, I(alt^2), I(alt^3)
# fm9 <- colext(~alt, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm10 <- colext(~alt+I(alt^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm11 <- colext(~alt+I(alt^2)+I(alt^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # add treeheight
# fm12 <- colext(~treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm13 <- colext(~alt+treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm14 <- colext(~alt+I(alt^2)+treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm15 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # add I(treeheight^2)
# fm16 <- colext(~treeheight+I(treeheight^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm17 <- colext(~alt+treeheight+I(treeheight^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm18 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm19 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # add I(alt^3)
# fm20 <- colext(~treeheight+I(treeheight^2)+I(treeheight^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm21 <- colext(~alt+treeheight+I(treeheight^2)+I(treeheight^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm22 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+I(treeheight^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm23 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2)+I(treeheight^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # add dbh
# fm24 <- colext(~dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm25 <- colext(~alt+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm26 <- colext(~alt+I(alt^2)+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm27 <- colext(~alt+I(alt^2)+I(alt^3)+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm28 <- colext(~treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm29 <- colext(~alt+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm30 <- colext(~alt+I(alt^2)+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm31 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm32 <- colext(~alt+treeheight+I(treeheight^2)+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm33 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm34 <- colext(~treeheight+I(treeheight^2)+I(treeheight^3)+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm35 <- colext(~alt+treeheight+I(treeheight^2)+I(treeheight^3)+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm36 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm37 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # add I(dbh^2)
# fm38 <- colext(~dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm39 <- colext(~alt+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm40 <- colext(~alt+I(alt^2)+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm41 <- colext(~alt+I(alt^2)+I(alt^3)+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm42 <- colext(~treeheight+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm43 <- colext(~alt+treeheight+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm44 <- colext(~alt+I(alt^2)+treeheight+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm45 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm46 <- colext(~alt+treeheight+I(treeheight^2)+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm47 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm48 <- colext(~treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm49 <- colext(~alt+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm50 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm51 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # add I(dbh^3)
# fm52 <- colext(~dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm53 <- colext(~alt+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm54 <- colext(~alt+I(alt^2)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm55 <- colext(~alt+I(alt^2)+I(alt^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm56 <- colext(~treeheight+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm57 <- colext(~alt+treeheight+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm58 <- colext(~alt+I(alt^2)+treeheight+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm59 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm60 <- colext(~alt+treeheight+I(treeheight^2)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm61 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm62 <- colext(~treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm63 <- colext(~alt+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm64 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm65 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# 
# # fit some more models with p(day+time+rain+wind+activity+ridge) and interactions
# fm66 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # hypothesis: initial occupancy depends on combination/interaction of alt and treeheight
# fm67 <- colext(~alt:treeheight:dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # hypothesis: initial occupancy depends on combination/interaction of alt and treeheight and dbh
# fm68 <- colext(~alt:treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # hypothesis: initial occupancy depends on combination/interaction of alt and treeheight plus dbh
# fm69 <- colext(~alt+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # hypothesis: initial occupancy depends on alt and treeheight and dbh
# fm70 <- colext(~alt+treeheight+alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# # hypothesis: initial occupancy depends on alt and treeheight additively and the interaction between alt and treeheight
# 
# # fit some random models with p(~day+time+rain+activity+ridge)
# fm71 <- colext(~1, ~1, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm72 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm73 <- colext(~alt:treeheight:dbh, ~1, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm74 <- colext(~alt:treeheight+dbh, ~1, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm75 <- colext(~alt:treeheight:dbh, ~1, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# 
# # put the fitted models in a fitList() and rank them by AIC in modSel()
# MTOR_psi_fitList <- fitList(fm8, fm9, fm10, fm11, fm12, fm13, fm14, fm15, fm16, 
#                             fm17, fm18, fm19, fm20, fm21, fm22, fm23, fm24, fm25, 
#                             fm26, fm27, fm28, fm29, fm30, fm31, fm32, fm33, fm34, 
#                             fm35, fm36, fm37, fm38, fm39, fm40, fm41, fm42, fm43, 
#                             fm44, fm45, fm46, fm47, fm48, fm49, fm50, fm51, fm52, 
#                             fm53, fm54, fm55, fm56, fm57, fm58, fm59, fm60, fm61, 
#                             fm62, fm63, fm64, fm65, fm66, fm67, fm68, fm69, fm70, 
#                             fm71, fm72, fm73, fm74, fm75) # do there is an alternative way to do this automatically?
# print(MTOR_psi_modSel <- modSel(MTOR_psi_fitList)) # extract the rows with delta < 2
# # best models: fm72 (and fm37, fm66, fm74)
# 
# # alternative workflow for initial occupancy with dredge()
# MTOR_psi_global_model <- colext(~alt*treeheight*dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# MTOR_psi_global_model_dredge <- dredge(MTOR_psi_global_model, trace = 2)
# print(MTOR_psi_global_model_dredge)
# get.models(MTOR_psi_global_model_dredge, subset = delta<2)
# # best models: psi(~1)g(~1)e(~1)p(~activity+day+rain+ridge+time+1), 
# #              psi(~1)g(~1)e(~1)p(~activity+day+rain+ridge+1),
# #              psi(~alt+treeheight+alt:treeheight+1)g(~1)e(~1)p(~activity+day+rain+ridge+time+1),
# # why there are different results than with the manually fitted models?

###### 4.3: fit models for extinction and colonization probability (e) for modSel ####

# with p(day+time+rain+wind+activity+ridge)
fm76 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm77 <- colext(~alt:treeheight, ~alt, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm78 <- colext(~alt:treeheight, ~1, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm79 <- colext(~alt:treeheight, ~alt, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)

# # with p(~day+time+rain+activity+ridge)
# fm80 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm81 <- colext(~alt:treeheight, ~alt, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm82 <- colext(~alt:treeheight, ~1, ~alt, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm83 <- colext(~alt:treeheight, ~alt, ~alt, ~day+time+rain+activity+ridge, data = umf, se = T)
# 
# # with psi(~1) and p(~day+time+rain+activity+ridge) - in dredge() this models were high ranked
# fm84 <- colext(~1, ~1, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm85 <- colext(~1, ~alt, ~1, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm86 <- colext(~1, ~1, ~alt, ~day+time+rain+activity+ridge, data = umf, se = T)
# fm87 <- colext(~1, ~alt, ~alt, ~day+time+rain+activity+ridge, data = umf, se = T)
# 
# # with psi(~1) and p(day+time+rain+wind+activity+ridge)
# fm88 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm89 <- colext(~1, ~alt, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm90 <- colext(~1, ~1, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# fm91 <- colext(~1, ~alt, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)

MTOR_g_e_fitList <- fitList(constant=fm76, expansion=fm77, contraction=fm78, shift=fm79,null=fm1)
                            # fm80, fm81, fm82, fm83, fm84, 
                            # fm85, fm86, fm87, fm88, fm89, fm90, fm91, fm1)
print(MTOR_g_e_modSel <- modSel(MTOR_g_e_fitList))

(MTOR_fm79_gof <- mb.gof.test(fm79, nsim = 100)) 


# best model over all (lowest AIC) is fm83 (psi(~alt:treeheight)g(~alt)e(~alt)p(~day+time+rain+activity+ridge)
# but fm81 (psi(~alt:treeheight)g(~alt)e(~1)p(~day+time+rain+activity+ridge) is 
# only deltaAIC 0.3 between them  

###### 4.5: alternative workflow with global model and MUMIn::dredge()

# # fit most complex model as global model, dregde to get best (AICc - change to AIC by rank = AIC) model
# MTOR_global_model <- colext(psiformula = ~alt*treeheight*dbh, 
#                             gammaformula = ~alt,
#                             epsilonformula = ~alt, 
#                             pformula = ~ day+time+wind+rain+activity, data = umf, se = T)
# MTOR_global_model_dredge <- dredge(MTOR_global_model, trace = 2)
# print(MTOR_global_model_dredge)
# get.models(MTOR_global_model_dredge, subset = delta<2)
# 
# 


# best models: psi(~1)g(~alt+1)e(~1)p(~activity+day+rain+time+1) AIC:3533.28, 
#              psi(~1)g(~alt+1)e(~1)p(~activity+day+rain+1) AIC:3534.49, 
#              psi(~1)g(~alt+1)e(~alt+1)p(~activity+day+rain+time+1) AIC:3533.882,
#              psi(~1)g(~alt+1)e(~alt+1)p(~activity+day+rain+1) AIC: 3534.71,
#              psi(~alt+1)g(~alt+1)e(~alt+1)p(~activity+day+rain+1) AIC: 3534.654
# can't understand why they are different to my manually fitted ones

##### 5: Goodness-of-fit-test (GOF) to check for any assumption violations ####

# Mackenzie-Bailey Goodness-of-Fit tests to simulate detection data on this specific model 
# to compare the binomial distribution with the real data/model with a simple 
# chi-square test

# mb-GOF - dont really know how to interpret (c-hat near 1)
(MTOR_fm83_gof <- mb.gof.test(fm83, nsim = 10)) # if there is everything alright - nsim = 1000
# c-hat (variance of  should be around 1 - why?
(MTOR_fm81_gof <- mb.gof.test(fm81, nsim = 10))

# GOF from cornell lab of ornitology and the best practise for ebird data
# (code from: https://doi90.github.io/lodestar/fitting-occupancy-models-with-unmarked.html)
fitstats <- function(fm83, method = "nonparboot") {
  observed <- getY(fm83@data)
  expected <- fitted(fm83)
  resids <- residuals(fm83, method = "nonparboot")
  sse <- sum(resids^2, na.rm = TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm = TRUE)
  out <- c(SSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(out)}

# calculation with unmarked::parboot
(pb <- parboot(fm83, fitstats, nsim = 10, report = TRUE, 
               method = "nonparboot"))

# plotting
par(mfrow = c(3,1))
plot(pb, main = "", xlab = c("SSE", "Chisq", "FT"))
# dev.off()

##### 6: view results and extract data ####

###### 6.1: from the best single model ####

# view models
summary(fm83) # best model based on AIC
summary(fm81) # second best model based on AIC

# return predicted values from model, specify by 'psi', 'col', 'ext', or 'det'
predict(fm83, type = 'psi')
predict(fm83, type = 'col')
# length(predict(fm83, type = 'col')$Predicted)/86*13 # what happens here?
predict(fm83, type = 'ext')
predict(fm83, type = 'det')
# I think thats not really what I like to do...

# extract coefficient and SE estimates from model
(data.frame(coef(fm83))) # estimates only 
SE(fm83) # standard errors only 

###### 6.2: model averaging of the best performing model-set ####

# model average on fm83 and fm81 (this were the best performing ones by AIC)
(MTOR_fm83_fm81_modavr <- model.avg(c(fm83, fm81), fit = T)) # fill in the list the models that should be averaged

# model average on dredged models with delta < 2
MTOR_modavr <- model.avg(get.models(MTOR_psi_global_model_dredge, subset = delta<2), fit = T) 
# here fill in the get.models-output
print(MTOR_modavr)
coef(MTOR_modavr) # function to print coefficients


##### 7: extraction of values and plotting ####

###### 7.3: model predictions from the best model fm83 ####
?predict

###### 7.2: model averages from fm83 and fm81 ####

# calculate averaged model predictions
MTOR_fm83_fm81_modavr_predict_psi <- modavgPred(c(fm83, fm81), modnames = c('fm83', 'fm81'), parm.type = 'psi', newdata = umf@siteCovs)
# some data to predict have to be added, but what kind of data?
MTOR_fm83_fm81_modavr_predict_psi

# store data
MTOR_fm83_fm81_modavr_predict_psi_data <- data.frame(psi_predicted = MTOR_fm83_fm81_modavr_predict_psi$mod.avg.pred,
                                                     lower = MTOR_fm83_fm81_modavr_predict_psi$lower.CL, 
                                                     upper = MTOR_fm83_fm81_modavr_predict_psi$upper.CL, 
                                                     point = siteCovsColext[,1], 
                                                     alt = siteCovsColext$alt, 
                                                     dbh = siteCovsColext$dbh, 
                                                     treeheight = siteCovsColext$treeheight)
MTOR_fm83_fm81_modavr_predict_psi_data # data.frame for initial occupancy predictions

# plotting, improve points (work with northings and eastings from siteCovs)
MTOR_fm83_fm81_modavr_predict_psi_data <- MTOR_fm83_fm81_modavr_predict_psi_data %>%
  mutate(x = rep(1:43, 2), # Adding pseudo-coordinates as if 100 sites are in 10 x 10 grid 
         y = rep(1:2, each = 43))

ggplot(data = MTOR_fm83_fm81_modavr_predict_psi_data, aes(x = x, y = y, fill = psi_predicted)) +
  geom_point(size = 10, pch = 22) +
  theme_classic() +
  scale_fill_gradient(low = "blue", high = "yellow", name = "Predicted\noccupancy") +
  ggtitle("Predicted occupancy")

# relationship with covariates 

###### 8: notes

