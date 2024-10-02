
##### Montserrat Forest Bird Counts - data import and analysis with 'unmarked' #####
##### written in June 2023 by  Filibert Heim, filibert.heim@posteo.de          #####

##### 1: load required packages ####

# install.packages('scales')
library(scales)
library(reshape)
library(lubridate)
library(unmarked)
library(AICcmodavg) # package for model selection and goodness-of-fit tests
library(MuMIn)
library(tidyverse)
# rename <- dplyr::rename
# select <- dplyr::select
# filter <- dplyr::filter

##### 2: load the prepared data and make last preparation #####

# set working directory and load prepared data from Steffen
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat')
setwd('C:/Users/sop/Documents/Steffen/RSPB/Montserrat')
load(file = 'MONTSERRAT_ANNUAL_DATA_INPUT2023.RData')
filter<-dplyr::filter
select<-dplyr::select
rename<-dplyr::rename
ls()

###### 2.1: set YEAR and SPECIES that should be analysed ####
YEAR <- 2023 # fill in the most recent year that data is available for
SPECIES <- c('MTOR') # fill in species the analysis should be made for
# SPECIES: 'MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA'

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

###### 2.3: remove unneeded object from global environment

rm(list = c('surveys2023', 'pointmax', 'table1', 'totals2023', 'birds', 
            'birds2023', 'activity', 'summary2023', 'nyears', 'MTOR_count', 
            'SURVEYDATA', 'rain_y', 'add', 'startdate', 'starttime', 'y', 
            'removal'))

###### 2.4: relable variate names/columns #####

obsCov <- select(obsCov, -Rain, -Year)
names(obsCov)[1:10] <- c('year','point','count','obs','skill','rain','wind','day',
                         'time','activity')

names(siteCov)[c(1:9, 21)] <- c('point', 'eastings', 'northings', 'habitat', 'dbh', 'distance', 
                         'treeheight', 'elevation', 'slope', 'alt') # go on with renaming

names(COUNTDATA) <- c('species', 'year', 'point', 'count', 'N')

###### 2.5: digression to calculate missing wind values as obsCov ####

obsCov_db <- obsCov_db %>% 
  select(-Time, -Rain, -MinOfskill, -Date, -CountOfObservername) %>%
  dplyr::rename(point = Point, count = Count)
obsCov_db$Wind <- ifelse(obsCov_db$Wind %in% c('light','calm'),'low','high')
obsCov <- obsCov %>% 
  left_join(x = obsCov, y = obsCov_db, by = c('year', 'point', 'count')) %>%
  select(-wind) %>% 
  rename(wind = Wind) %>%
  mutate(wind = as.character(wind)) %>% 
  replace(is.na(.), 'low') %>%
  mutate(wind = as.factor(wind))

###### 2.6: digression to postpone surveys if some weren't made ######

# filter for SPECIES and align with obsCov to get information from surveys that weren't made
COUNTDATA <- COUNTDATA %>% filter(species == SPECIES) %>% # choose species
  left_join(obsCov %>% mutate(point = as.numeric(point)), by = c('point', 'count', 'year')) %>%
  mutate(N = ifelse(obs == 'NA', 'NA', N)) %>% # insert 'NA' to COUNTDATA were no surveys were made
  select(-obs, -skill, -rain, -time, -activity, -wind) # remove unneeded columns

# transfer information of surveys that weren't made also to obsCov data
obsCov <- COUNTDATA %>%
  left_join(obsCov %>% mutate(point = as.numeric(point)), by = c("point", 'count', 'year')) %>%
  select(-species, -N, -day.x) %>%
  rename(day = day.y)

# filter COUNTDATA for all the surveys in one year at one point, sort by day and rename counts
COUNTDATA <- COUNTDATA %>% 
  group_by(year, point) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(count = as.numeric(rep(1:3))) %>% 
  select(-day) %>% 
  ungroup()

# filter obsCov for all the surveys in one year at one point, sort by day and rename counts
obsCov <- obsCov %>% 
  group_by(year, point) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(count = as.numeric(rep(1:3))) %>% 
  ungroup()

##### 3: provide and prepare data for unmarkedMultFrame ####

numPrimary <- length(unique(COUNTDATA$year)) # calculates number of years

occdataColext <- COUNTDATA %>%
  mutate(occupancy = ifelse(N > 0, 1, 0)) %>% # convert in detection/non-detection
  mutate(season = paste(year, count, sep = '_')) %>% # connect year, count to string
  select(-species, -N, -year, -count) %>% # remove unneeded columns 
  spread(key = season, value = occupancy) %>% # spread a key-value pair across multiple columns
  arrange(point) # sort rows in order of point number

siteCovsColext <- siteCov %>% mutate(point = as.numeric(point)) %>% 
  arrange(point) # sort rows in order of point number

obsCovsColext <- COUNTDATA %>% filter(species == SPECIES) %>%
  left_join((obsCov %>% mutate(point = as.numeric(point))), by = c('point', 'count', 'year')) %>% 
  select(-N, -species, -skill) %>% # remove columns
  arrange(point, year, count) # sort rows in order of year, point, count

yearlySiteCovsYear <- COUNTDATA %>% filter(species == SPECIES) %>% 
  group_by(point,year) %>%
  summarise(occu = as.factor(mean(year) - 2010)) %>%
  spread(key = year, value = occu) %>% # spread a key-value pair across multiple columns
  arrange(point)
yearlySiteCovsColext <- list(year=yearlySiteCovsYear[,2:14])
  
###### 3.2: check dimensions ####

dim(occdataColext)
dim(siteCovsColext)
dim(obsCovsColext)/(3*numPrimary)
dim(yearlySiteCovsYear)

###### 3.3: create unmarkedMultFrame ####

umf <- unmarkedMultFrame(y = occdataColext[,2:40], siteCovs = siteCovsColext,yearlySiteCovs=yearlySiteCovsColext, 
                         obsCovs = obsCovsColext, numPrimary = numPrimary)
summary(umf) # looks good
str(umf) # looks good, but the values have to be scaled

###### 3.4: scale numeric variables to fitting problems ####

siteCovs(umf)[c(2,3,5:17,19,21:38)] <- scale(siteCovs(umf)[c(2,3,5:17,19,21:38)])
obsCovs(umf)[c(5:7)] <- scale(obsCovs(umf)[c(5:7)])
str(umf)
summary(umf)

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
# best submodel for p(): ~day+time+rain+wind+activity+ridge, go on with this 

###### 4.2: fit models for initial occupancy psi() first for modSel ####

fm8 <- colext(~alt, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm9 <- colext(~treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm10 <- colext(~alt+treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm11 <- colext(~dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm12 <- colext(~alt+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm13 <- colext(~treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm14 <- colext(~alt+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm15 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm16 <- colext(~alt:treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)

# put the fitted models in a fitList() and rank them by AIC in modSel()
MTOR_psi_fitList <- fitList(fm7, fm8, fm9, fm10, fm11, fm12,fm13, fm14, fm15, fm16) # do there is an alternative way to do this automatically?
print(MTOR_psi_modSel <- modSel(MTOR_psi_fitList)) 
# best sub-model for psi(): ~alt:treeheight, go on with this

###### 4.3: fit models for extinction and colonization probability for modSel ####

fm17 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm18 <- colext(~alt:treeheight, ~alt, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm19 <- colext(~alt:treeheight, ~1, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm20 <- colext(~alt:treeheight, ~alt, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm21 <- colext(~alt:treeheight, ~year, ~year, ~day+time+rain+wind+activity+ridge, data = umf, se = T)  ## this model will exclude the possibility that observed changes are just annual changes

# put the fitted models in a fitList() and rank them by AIC in modSel()
MTOR_g_e_fitList <- fitList(constant = fm17, expansion = fm18, contraction = fm19,
                            shift = fm20, year = fm21, null = fm1)
print(MTOR_g_e_modSel <- modSel(MTOR_g_e_fitList))
# shift model (fm20) has lowest AIC (but not far away from expansion model (fm18))

best_model <- fm20 # store best model as object for ongoing code

##### 6: Goodness-of-Fit Test ####

# mb-GOF
(MTOR_best_model_gof <- mb.gof.test(best_model, nsim = 15)) ## p-value of 0.1 is ok

# GOF from cornell lab of ornitology and the best practise for ebird data
# (code from: https://doi90.github.io/lodestar/fitting-occupancy-models-with-unmarked.html)
fitstats <- function(best_model, method = "nonparboot") {
  observed <- getY(best_model@data)
  expected <- fitted(best_model)
  resids <- residuals(best_model, method = "nonparboot")
  sse <- sum(resids^2, na.rm = TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm = TRUE)
  out <- c(SSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(out)}

# calculation with unmarked::parboot
(pb <- parboot(best_model, fitstats, nsim = 100, report = TRUE, 
               method = "nonparboot")) ## prob of 0.9 is ok

# plot  gof-test
par(mfrow = c(3,1))
plot(pb, main = "", xlab = c("SSE", "Chisq", "FT"))
# dev.off()

##### 7: view results and extract data ####

###### 7.1: extract data from the best single model - no model averaging needed in this thesis ####

# view models
summary(best_model) # best model based on AIC

# extract coefficient and SE estimates from model
data.frame(coef(best_model)) # estimates only 
SE(best_model) # standard errors only 

# extract everything
toExportpsi <- summary(best_model@estimates@estimates$psi) %>%
  mutate(Parameter=names(best_model@estimates@estimates$psi@estimates)) %>%
  mutate(Component="Initial Occupancy")

toExportdet <- summary(best_model@estimates@estimates$det) %>%
  mutate(Parameter=names(best_model@estimates@estimates$det@estimates)) %>%
  mutate(Component="Detection probability")

toExportcol <- summary(best_model@estimates@estimates$col) %>%
  mutate(Parameter=names(best_model@estimates@estimates$col@estimates)) %>%
  mutate(Component="Colonization")

toExportext <- summary(best_model@estimates@estimates$ext) %>%
  mutate(Parameter=names(best_model@estimates@estimates$ext@estimates)) %>%
  mutate(Component="Extinction")

toExport<-bind_rows(toExportpsi,toExportdet,toExportcol,toExportext)
fwrite(toExport,sprintf("%s_ColExt_parameter_estimates.csv",SPECIES))


##### 8: extraction of values and plotting ####

###### 8.1 plot relationship between colonization and extinction with altitude #####
# code mainly taken from: https://cran.r-project.org/web/packages/unmarked/vignettes/colext.html and steffen

# create input data.frame with values for the prediction, 
nd <- data.frame(day = 0, time = 0, rain = 0, wind = 'low', activity = 1, 
                   ridge = 0, treeheight = mean(siteCovsColext$treeheight), # mean treeheight used (steffen used 15, why?)
                   alt = seq(from = min(umf@siteCovs$alt), 
                                  to = max(umf@siteCovs$alt), by = 0.05),
                   altitude = rescale(seq(from = min(umf@siteCovs$alt), to = max(umf@siteCovs$alt), by = 0.05), 
                                 to = c(min(siteCovsColext$alt), max(siteCovsColext$alt)), # to = output range as c()
                                 from = c(min(umf@siteCovs$alt), max(umf@siteCovs$alt)))) # from = input range as c()
# alt: scaled altitude for prediction, altitude: rescaled altitude for plotting

# predict values by input data from nd, add column with rescaled altitude
pred_ext <- predict(best_model, type = 'ext', newdata = nd) %>%
  mutate(Type = 'Extinction', Altitude = nd$altitude)
pred_col <- predict(best_model, type = 'col', newdata = nd) %>%
  mutate(Type = 'Colonization', Altitude = nd$altitude)

# plot colonization~altitude and extinction~altidude together
bind_rows(pred_ext,pred_col) %>%
  ggplot(aes(x = Altitude, y = Predicted, colour = Type, fill = Type)) +
  geom_line(linewidth = 3) +
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.2) +
  
  labs(x="Elevation (m above sea level)", y="Probability") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=20, color="black"), 
        axis.title=element_text(size=24),
        legend.text=element_text(size=20),
        legend.title = element_text(size=22),
        legend.position=c(0.2,0.88),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))



###### 8.1 plot estimated occupancy for time series ####

# using empirical bayes estimates of occupancy averaged across all points
occupancy <- ranef(best_model)

# manipulate data structure, summarize and store data in data.frame
occupancy_data <- as.data.frame(bup(occupancy, stat = 'mean')) %>% # posterior mean by stat = ''
  gather(key = 'Year',value = 'occu') %>% # transfer from wide to long format 
  mutate(Year = as.numeric(str_replace(Year,'V',''))) %>% # delete V before number of season
  mutate(Year = Year + 2010) %>% # calculate year 
  group_by(Year) %>%
  summarise(Occupancy = mean(occu)) # calculate the mean occupancy about all points in one year 

# calculate confidence intervals and store them in occupancy_data data.frame for plotting
occupancy_confint <- confint(occupancy, level = 0.95) # 95% CI
occupancy_data$lower_cl <- apply(occupancy_confint[,1,],2, mean)
occupancy_data$upper_cl <- apply(occupancy_confint[,2,],2, mean)

# romove data from 2020 because in this year surveys didn't take place 
occupancy_data <- occupancy_data %>% filter(!Year == '2020')

# plot occupancy over time from occupancy_data data.frame
ggplot(data = occupancy_data, aes(x = Year, y = Occupancy)) +
  geom_point(size = 3, pch = 16, colour="firebrick") +
  geom_errorbar(aes(ymin = lower_cl,ymax = upper_cl), width = 0.2) +
  scale_y_continuous(limits = c(0, 1.1), breaks=seq(0,1,0.2),labels=seq(0,1,0.2))  +
  scale_x_continuous(limits = c(2010, 2024), breaks=seq(2011,2023,2),labels=seq(2011,2023,2)) +
  labs(x="Year", y="Mean occupancy probability") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=20, color="black"), 
        axis.title=element_text(size=24),
        legend.text=element_text(size=20),
        legend.title = element_text(size=22),
        legend.position=c(0.2,0.88),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))

##### 9: notes ####

# how to export the stuff? 
# how to make the plots more fancy?
# should I add intercepts and stuff like that to my results in the thesis?
# how to export model comparison tables?
