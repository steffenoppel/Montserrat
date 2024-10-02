
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

yearlySiteCovsYear <- COUNTDATA %>% filter(species == 'MTOR') %>% # choose species
  group_by(point,year) %>%
  summarise(occu=as.factor(mean(year)-2010)) %>%
  spread(key = year, value = occu) %>% # spread a key-value pair across multiple columns
  arrange(point)

yearlySiteCovsColext <- list(year=yearlySiteCovsYear[,2:14])
  
###### 3.2: check dimensions ####

dim(occdataColext)
dim(siteCovsColext)
dim(obsCovsColext)/(3*numPrimary)

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
# best model: psi(.)g(.)e(.)p(~day+time+rain+wind+activity+ridge), go on with this 


### EXPLORE MODELS FOR INITIAL OCCUPANCY
fm9 <- colext(~alt, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm12 <- colext(~treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm13 <- colext(~alt+treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm24 <- colext(~dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm25 <- colext(~alt+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm28 <- colext(~treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm29 <- colext(~alt+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm66 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm68 <- colext(~alt:treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)

# # put the fitted models in a fitList() and rank them by AIC in modSel()
MTOR_psi_fitList <- fitList(fm7, fm9, fm12, fm13, fm24, fm25,fm28, fm29, fm66, fm68) # do there is an alternative way to do this automatically?
print(MTOR_psi_modSel <- modSel(MTOR_psi_fitList)) # extract the rows with delta < 2
# # best models: fm66


###### 4.3: fit models for extinction and colonization probability (e) for modSel ####

# with p(day+time+rain+wind+activity+ridge)
fm76 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm77 <- colext(~alt:treeheight, ~alt, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm78 <- colext(~alt:treeheight, ~1, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm79 <- colext(~alt:treeheight, ~alt, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
year <- colext(~alt:treeheight, ~year, ~year, ~day+time+rain+wind+activity+ridge, data = umf, se = T)  ## this model will exclude the possibility that observed changes are just annual changes

MTOR_g_e_fitList <- fitList(constant=fm76, expansion=fm77, contraction=fm78, shift=fm79,null=fm1,year=year)
print(MTOR_g_e_modSel <- modSel(MTOR_g_e_fitList))


# Mackenzie-Bailey Goodness-of-Fit tests to simulate detection data on this specific model 
# to compare the binomial distribution with the real data/model with a simple 
# chi-square test

# mb-GOF - dont really know how to interpret (c-hat near 1)
(MTOR_fm79_gof <- mb.gof.test(fm79, nsim = 100)) ## p-value of 0.1 is ok

# GOF from cornell lab of ornitology and the best practise for ebird data
# (code from: https://doi90.github.io/lodestar/fitting-occupancy-models-with-unmarked.html)
fitstats <- function(fm79, method = "nonparboot") {
  observed <- getY(fm79@data)
  expected <- fitted(fm79)
  resids <- residuals(fm79, method = "nonparboot")
  sse <- sum(resids^2, na.rm = TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm = TRUE)
  out <- c(SSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(out)}

# calculation with unmarked::parboot
(pb <- parboot(fm79, fitstats, nsim = 100, report = TRUE, 
               method = "nonparboot")) ## prob of 0.9 is ok

# plotting
par(mfrow = c(3,1))
plot(pb, main = "", xlab = c("SSE", "Chisq", "FT"))
# dev.off()

##### 6: view results and extract data ####

###### 6.1: from the best single model ####

# view models
summary(fm79) # best model based on AIC

# extract coefficient and SE estimates from model
(data.frame(coef(fm79))) # estimates only 
SE(fm79) # standard errors only 


##### 7: extraction of values and plotting ####



##### PLOT COLONISATION AND EXTINCTION WITH ELEVATION

# relationship with covariates 
## https://cran.r-project.org/web/packages/unmarked/vignettes/colext.html
## predict returns the predictions along with standard errors and confidence intervals. These can be used to create plots. The with function is used to simplify the process of requesting the columns of data.frame returned by predict.

nd <- data.frame(day=0,time=0,rain=0,wind="low",activity=1,ridge=0,treeheight=15,
                 alt=seq(-1.5266,2.3832,0.065),  # the scaled version used in the model
                 Elevation=seq(min(siteCov$alt),max(siteCov$alt),10))
E.ext <- predict(fm79, type='ext', newdata=nd) %>%
  mutate(Type="Extinction", Elevation=nd$Elevation)
E.col <- predict(fm79, type='col', newdata=nd) %>%
  mutate(Type="Colonization", Elevation=nd$Elevation)

bind_rows(E.ext,E.col) %>%
  

  ggplot(aes(x = Elevation, y = Predicted, colour=Type, fill=Type)) +
  geom_line(linewidth = 3) +
  geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.2) +
  theme_classic() +
  ggtitle("Relationship of colonisation and extinction probability with elevation")



#### EXTRACT REALISED ESTIMATES OF OCCUPANCY AND PLOT TIME SERIES ####
## using empirical Bayes estimates of occupancy averaged across all points
out<-ranef(fm79)
occuplot<-as.data.frame(bup(out, stat="mean")) %>%          # Posterior mean
  gather(key="Year",value="occu") %>%
  mutate(Year=as.numeric(str_replace(Year,"V",""))) %>%
  mutate(Year=Year+2010) %>%
  group_by(Year) %>%
  summarise(Occupancy=mean(occu))

occuconfint<-confint(out, level=0.95) # 95% CI
occuplot$lcl<-apply(occuconfint[,1,],2,mean)
occuplot$ucl<-apply(occuconfint[,2,],2,mean)



# Plot Occupancy over time
ggplot(data = occuplot, aes(x = Year, y = Occupancy)) +
  geom_point(size = 3, pch = 16) +
  geom_errorbar(aes(ymin=lcl,ymax=ucl), width=0.2) +
  theme_classic() +
  ggtitle("Predicted occupancy")
