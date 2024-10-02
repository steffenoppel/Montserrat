
##### Montserrat Forest Bird Counts - data import and analysis with 'unmarked' #####
##### written in June 2023 by  Filibert Heim, filibert.heim@posteo.de          #####

##### 1: load required packages ####

library(RODBC)
library(reshape)
library(tidyverse)
library(lubridate)
filter<-dplyr::filter
select<-dplyr::select
rename<-dplyr::rename
library(unmarked)
# install.packages('AICcmodavg')
library(AICcmodavg) # package for model selection and GOF 
# install.packages('MuMIn')
library(MuMIn)

##### 2: set YEAR and SPECIES that should be analysed ####
YEAR <- 2023 # fill in the most recent year that data is available for
SPECIES <- c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH',
             'SNPI','CAEL','BANA') # fill in species the analysis should be made for

###### 2.1: load the prepared data and make last preparation #####

# load prepared data from Steffen 
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat/GitHub_Steffen_Montserrat/')
load(file = 'MONTSERRAT_ANNUAL_DATA_INPUT2023.RData')
ls()

# load rain and wind data from the database for obsCovs and yearlySiteCovs
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/Datengrundlage/')
mt_db <- odbcConnectAccess2007('Montserrat_Birds_2023.accdb')
rain_y <- sqlQuery(mt_db, 'SELECT * FROM yearlyRainCov')
siteCov_db <- sqlQuery(mt_db, 'SELECT * FROM SiteCov')
obsCov_db <- sqlQuery(mt_db, 'SELECT * FROM obsCovariates_simple')
obsCov_db$Point <- as.character(obsCov_db$Point)
obsCov_db$Wind <- as.factor(obsCov_db$Wind)
odbcClose(mt_db)

###### 2.2: check the data and remove unneeded objects #####

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

###### 2.4: digression to calculate missing wind / PETH values and yearly siteCovs ####

# fill in NA's in COUNTDATA for 2020 that didn't take place 
# COUNTDATA$N[COUNTDATA$year == '2020'] <- NA

# obsCov: wind 
obsCov_db <- obsCov_db %>% 
  select(-Time, -Rain, -MinOfskill, -Date, -CountOfObservername) %>%
  rename(point = Point, count = Count)
obsCov_db$Wind <- ifelse(obsCov_db$Wind %in% c('light','calm'),'low','high')

obsCov <- obsCov %>% 
  left_join(x = obsCov, y = obsCov_db, by = c('year', 'point', 'count')) %>%
  select(-wind) %>% 
  rename(wind = Wind) %>%
  mutate(wind = as.character(wind)) %>% # attention! thats not good - but I see 
  # no other possibility than replace all NA's by 'low' wind conditions, I also 
  # could search for data from nearby points and the same day // now used 'low'
  # because its the most common value # plot(obsCovsColext$wind)
  replace(is.na(.), 'low') %>%
  mutate(wind = as.factor(wind))

# obsCov: PETH
COUNTDATA_PETH <- COUNTDATA %>%
  subset(species == 'PETH') %>%
  rename(PETH_N = N)
COUNTDATA_PETH$point <- as.character(COUNTDATA_PETH$point)
obsCov <- obsCov %>%
  left_join(x = obsCov, y = COUNTDATA_PETH, by = c('year', 'point', 'count')) %>%
  select(-species)
# obsCov[is.na(obsCov)] <- 0 # not sure, if this is okay like this or simply not 
# needed because this data is removed automatically later

# yearlySiteCov: PETH (mean and max)
COUNTDATA_PETH_mean <- COUNTDATA_PETH %>% 
  aggregate(PETH_N ~ year + point, FUN = mean) %>% 
  rename(PETH_mean = PETH_N)
COUNTDATA_PETH_max <- COUNTDATA_PETH %>% 
  aggregate(PETH_N ~ year + point, FUN = max) %>% 
  rename(PETH_max = PETH_N)
yearlySiteCovPETH_max <- COUNTDATA_PETH_max %>% 
  spread(key = year, value = PETH_max) %>% # should year 2020 be populated with NA's for missing surveys?
  mutate(point = as.numeric(point)) %>%
  # mutate('2020' = NA, .after = '2019') %>% # something went wrong, thats why this 
  # line had to be added
  arrange(point)
yearlySiteCovPETH_mean <- COUNTDATA_PETH_mean %>% 
  spread(key = year, value = PETH_mean) %>% # should year 2020 be populated with NA's for missing surveys?
  mutate(point = as.numeric(point)) %>%
  # mutate('2020' = NA, .after = '2019') %>% # something went wrong, thats why this 
  # line had to be added
  arrange(point)

# yearlySiteCov: rain_y
yearlySiteCovRain_yearly <- rain_y %>% 
  subset(subset = Point %in% siteCov$point) %>%
  spread(key = Year, value = SumOfrainfall) %>%
  mutate('2020' = 'NA', '2021' = 'NA', '2022' = 'NA', '2023' = 'NA') %>% 
  rename(point = Point) %>%
  mutate_all(~replace(.,is.na(.),0)) %>% # replace all 'NA', remove, when data available
  select(-'2020', -'2021', -'2022', -'2023') %>%
  mutate('2020' = 0, '2021' = 0, '2022' = 0, '2023' = 0)
yearlySiteCovRain_yearly[,c(1,11:14)] <- as.numeric(unlist(yearlySiteCovRain_yearly[,c(1,11:14)]))
str(yearlySiteCovRain_yearly)

###### 2.5: remove unneeded objects from global environment
rm(MTOR_count)
rm(nyears)
rm(summary2023)
rm(list = c('surveys2023', 'pointmax', 'table1', 'totals2023', 'birds', 
            'birds2023', 'activity'))

##### 3: provide an unmarkedFrame for multi season analysis with 'colext()' ####

numPrimary <- YEAR - 2010 

occdataColext <- COUNTDATA %>% filter(species == 'MTOR') %>%
  mutate(occupancy = ifelse(N > 0, 1, 0)) %>% # convert in detection/non-detection
  mutate(season = paste(year, count, sep = '_')) %>% # connect year, count to string
  select(-species, -N, -year, -count) %>% # remove unneeded columns 
  spread(key = season, value = occupancy) %>% # spread a key-value pair across multiple columns
  arrange(point)

siteCovsColext <- siteCov %>% mutate(point = as.numeric(point)) %>% 
  arrange(point) # sort rows in order of point number

yearlySiteCovColext <- list(PETH_max = yearlySiteCovPETH_max[,2:14], 
                            PETH_mean = yearlySiteCovPETH_mean[,2:14], 
                            rain_yearly = yearlySiteCovRain_yearly[,2:14])

obsCovsColext <- COUNTDATA %>% filter(species == 'MTOR') %>%
  left_join((obsCov %>% mutate(point = as.numeric(point))), by = c('point', 'count', 'year')) %>% 
  select(-N, -species, -skill) %>% # remove columns
  arrange(point, year, count) # sort rows in order of year, point, count

###### 3.2: check dimensions ####

dim(occdataColext)
dim(siteCovsColext)
dim(obsCovsColext)/(3*numPrimary)
dim(yearlySiteCovColext$PETH_max)
dim(yearlySiteCovColext$PETH_mean)
dim(yearlySiteCovColext$rain_yearly)

###### 3.3: create unmarkedMultFrame ####

umf <- unmarkedMultFrame(y = occdataColext[,2:40], siteCovs = siteCovsColext, 
                         obsCovs = obsCovsColext, numPrimary = numPrimary, 
                         yearlySiteCovs = yearlySiteCovColext)
summary(umf) # looks good
str(umf) # looks good, but the values have to be scaled

###### 3.4: scale numeric variables to fitting problems ####

siteCovs(umf)[c(2,3,5:17,19,21:38)] <- scale(siteCovs(umf)[c(2,3,5:17,19,21:38)])
obsCovs(umf)[c(5:7)] <- scale(obsCovs(umf)[c(5:7)])
yearlySiteCovs(umf)[c(1:3)] <- scale(yearlySiteCovs(umf)[c(1:3)])
str(umf)

##### 4: multi-season, dynamic one species occupancy model ####

###### 4.1: fit models for detection probabilty p() first for modSel ####

# colext(~1, ~1, ~1, ~day + time + wind + ridge + rain + activity, data = umf, se = T) most complex model

# fit models for detection probabilty p() manually, go on with fitList(), modSel()
fm1 <- colext(~1, ~1, ~1, ~1, data = umf, se =T)
fm2 <- colext(~1, ~1, ~1, ~day, data = umf, se =T)
fm3 <- colext(~1, ~1, ~1, ~day+time, data = umf, se =T)
fm4 <- colext(~1, ~1, ~1, ~day+time+rain, data = umf, se =T)
fm5 <- colext(~1, ~1, ~1, ~day+time+rain+wind, data = umf, se =T)
fm6 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se =T)
fm7 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se =T)

MTORp_fitList <- fitList('psi(.)g(.)e(.)p(~day)' = fm2, 
                         'psi(.)g(.)e(.)p(~day+time)' = fm3, 
                         'psi(.)g(.)e(.)p(~day+time+rain)' = fm4,
                         'psi(.)g(.)e(.)p(~day+time+rain+wind)' = fm5,
                         'psi(.)g(.)e(.)p(~day+time+rain+wind+activity)' = fm6,
                         'psi(.)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm7)
(MTORp_modSel <- modSel(MTORp_fitList))

# fit global model for detection probability p(), modSel with dredge(), get.models()
MTOR_p_global_model <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se =T)
MTOR_p_global_model_dredge <- dredge(MTOR_p_global_model, trace = 2)
print(MTOR_p_global_model_dredge)
get.models(MTOR_p_global_model_dredge, subset = delta <2)

###### 4.2: fit models for initial occupancy psi() first for modSel ####

fm8 <- colext(~1, ~1, ~1, ~1, data = umf, se = T) # intercept only model
# add alt, I(alt^2), I(alt^3)
fm9 <- colext(~alt, ~1, ~1, ~1, data = umf, se = T)
fm10 <- colext(~alt+I(alt^2), ~1, ~1, ~1, data = umf, se = T)
fm11 <- colext(~alt+I(alt^2)+I(alt^3), ~1, ~1, ~1, data = umf, se = T)
# add treeheight
fm12 <- colext(~treeheight, ~1, ~1, ~1, data = umf, se = T)
fm13 <- colext(~alt+treeheight, ~1, ~1, ~1, data = umf, se = T)
fm14 <- colext(~alt+I(alt^2)+treeheight, ~1, ~1, ~1, data = umf, se = T)
fm15 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight, ~1, ~1, ~1, data = umf, se = T)
# add I(treeheight^2)
fm16 <- colext(~treeheight+I(treeheight^2), ~1, ~1, ~1, data = umf, se = T)
fm17 <- colext(~alt+treeheight+I(treeheight^2), ~1, ~1, ~1, data = umf, se = T)
fm18 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2), ~1, ~1, ~1, data = umf, se = T)
fm19 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2), ~1, ~1, ~1, data = umf, se = T)
# add I(alt^3)
fm20 <- colext(~treeheight+I(treeheight^2)+I(treeheight^3), ~1, ~1, ~1, data = umf, se = T)
fm21 <- colext(~alt+treeheight+I(treeheight^2)+I(treeheight^3), ~1, ~1, ~1, data = umf, se = T)
fm22 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+I(treeheight^3), ~1, ~1, ~1, data = umf, se = T)
fm23 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2)+I(treeheight^3), ~1, ~1, ~1, data = umf, se = T)
# add dbh
fm24 <- colext(~dbh, ~1, ~1, ~1, data = umf, se = T)
fm25 <- colext(~alt+dbh, ~1, ~1, ~1, data = umf, se = T)
fm26 <- colext(~alt+I(alt^2)+dbh, ~1, ~1, ~1, data = umf, se = T)
fm27 <- colext(~alt+I(alt^2)+I(alt^3)+dbh, ~1, ~1, ~1, data = umf, se = T)
fm28 <- colext(~treeheight+dbh, ~1, ~1, ~1, data = umf, se = T)
fm29 <- colext(~alt+treeheight+dbh, ~1, ~1, ~1, data = umf, se = T)
fm30 <- colext(~alt+I(alt^2)+treeheight+dbh, ~1, ~1, ~1, data = umf, se = T)
fm31 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+dbh, ~1, ~1, ~1, data = umf, se = T)
fm32 <- colext(~alt+treeheight+I(treeheight^2)+dbh, ~1, ~1, ~1, data = umf, se = T)
fm33 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+dbh, ~1, ~1, ~1, data = umf, se = T)
fm34 <- colext(~treeheight+I(treeheight^2)+I(treeheight^3)+dbh, ~1, ~1, ~1, data = umf, se = T)
fm35 <- colext(~alt+treeheight+I(treeheight^2)+I(treeheight^3)+dbh, ~1, ~1, ~1, data = umf, se = T)
fm36 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh, ~1, ~1, ~1, data = umf, se = T)
fm37 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh, ~1, ~1, ~1, data = umf, se = T)
# add I(dbh^2)
fm38 <- colext(~dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm39 <- colext(~alt+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm40 <- colext(~alt+I(alt^2)+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm41 <- colext(~alt+I(alt^2)+I(alt^3)+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm42 <- colext(~treeheight+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm43 <- colext(~alt+treeheight+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm44 <- colext(~alt+I(alt^2)+treeheight+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm45 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm46 <- colext(~alt+treeheight+I(treeheight^2)+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm47 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm48 <- colext(~treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm49 <- colext(~alt+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm50 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
fm51 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2), ~1, ~1, ~1, data = umf, se = T)
# add I(dbh^3)
fm52 <- colext(~dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm53 <- colext(~alt+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm54 <- colext(~alt+I(alt^2)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm55 <- colext(~alt+I(alt^2)+I(alt^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm56 <- colext(~treeheight+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm57 <- colext(~alt+treeheight+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm58 <- colext(~alt+I(alt^2)+treeheight+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm59 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm60 <- colext(~alt+treeheight+I(treeheight^2)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm61 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm62 <- colext(~treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm63 <- colext(~alt+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm64 <- colext(~alt+I(alt^2)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)
fm65 <- colext(~alt+I(alt^2)+I(alt^3)+treeheight+I(treeheight^2)+I(treeheight^3)+dbh+I(dbh^2)+I(dbh^3), ~1, ~1, ~1, data = umf, se = T)

# put the fitted models in a fitList() and rank them by AIC
MTOR_psi_fitList <- fitList(fm8, fm9, fm10, fm11, fm12, fm13, fm14, fm15, fm16, 
                            fm17, fm18, fm19, fm20, fm21, fm22, fm23, fm24, fm25, 
                            fm26, fm27, fm28, fm29, fm30, fm31, fm32, fm33, fm34, 
                            fm35, fm36, fm37, fm38, fm39, fm40, fm41, fm42, fm43, 
                            fm44, fm45, fm46, fm47, fm48, fm49, fm50, fm51, fm52, 
                            fm53, fm54, fm55, fm56, fm57, fm58, fm59, fm60, fm61, 
                            fm62, fm63, fm64, fm65) # do there is an alternative way to do this automatically?
print(MTOR_psi_modSel <- modSel(MTOR_psi_fitList)) # extract the rows with delta < 2

# alternative workflow for initial occupancy with dredge()
MTOR_psi_global_model <- colext(~alt*treeheight*dbh, ~1, ~1, ~1, data = umf, se = T)
MTOR_psi_global_model_dredge <- dredge(MTOR_psi_global_model, trace = 2)
print(MTOR_psi_global_model_dredge)
get.models(MTOR_psi_global_model_dredge, subset = delta<2)

###### 4.3: fit models for extinction probability (e) for modSel ####

###### 4.4: fit models for colonization probability (p) for modSel ####

###### 4.5: alternative workflow with global model for model
#           fitting and model comparison with MUMIn::dredge()

# fit most complex model as global model, dregde to get best (AIC) model
global.model_MTOR <- colext(psiformula = ~ alt + treeheight * dbh, 
                            gammaformula = ~ alt + rain_yearly + PETH_max * PETH_mean + treeheight * dbh,
                            epsilonformula = ~ alt + rain_yearly + PETH_max * PETH_mean + treeheight * dbh, 
                            pformula = ~ day * time + wind * rain + activity + PETH_N, data = umf, se = T)

# fit small global model as global model, dredge() and the use get.model() to get best (AIC) model(s) 
global.model_MTOR_small <- colext(psiformula = ~ alt, 
                            gammaformula = ~ alt * rain_yearly * PETH_mean,
                            epsilonformula = ~ alt * rain_yearly * PETH_mean, 
                            pformula = ~ day * time * wind * rain * activity, 
                            data = umf, se = T)
global.model_MTOR_small_dredge <- dredge(global.model_MTOR_small, trace = 2)


# description: 
## psiformula - initial occupancy, without yearlySiteCovs!
## pformula - detection probability
## gammaformula - colonization probability
## epsilonformula - extinction probability

# Aikaike Information Criterion to measure the quality of fitted 
# values (but also accounts for the number of variables)
# -> the highter the number, the better the model fits the data 

##### 5: Goodness-of-fit-test (GOF) to check for any assumption violations ####

# Mackenzie-Bailey Goodness-of-Fit tests to simulate detection data on this specific model 
# to compare the binomial distribution with the real data/model with a simple 
# (pearsons?-)chi-square test

global.model_MTOR_small_GOF <- mb.gof.test(global.model_MTOR_small, nsim = 10, 
                                           plot.hist = FALSE)

test.model_MTOR_GOF <- mb.gof.test(test.model_MTOR, nsim = 50, plot.hist = TRUE)
print(test.model_MTOR_GOF)

### ERROR! Error in dimnames(x) <- dn : 
# length of 'dimnames' [2] not equal to array extent
# In addition: There were 50 or more warnings (use warnings() to see the first 50)


###### 4.6: model averaging of the best performing model-set 

# model average on models with delta < 2
MTOR_modavr <- model.avg(get.models(MTOR_psi_global_model_dredge, subset = delta<2), fit = T) 
# here fill in the get.models-output
print(MTOR_modavr)
coef(MTOR_modavr) # function to print coefficients

###### 6.2: view results

summary(global.model_MTOR)
summary(test.model_MTOR)


###### 7: extraction of values and plotting ####

predict()

###### 8: notes

