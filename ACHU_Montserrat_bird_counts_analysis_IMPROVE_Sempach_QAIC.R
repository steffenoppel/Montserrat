
##### Montserrat Forest Bird Counts - data import and analysis with 'unmarked' #####
##### written in Sep 2024 by  Filibert Heim, filibert.heim@posteo.de          #####

##### 1: load required packages ####

# install.packages('scales')
library(scales)
# install.packages('data.table')
library(data.table)
library(reshape)
library(lubridate)
library(unmarked)
library(AICcmodavg) # package for model selection and goodness-of-fit tests
library(MuMIn)
library(tidyverse)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename

##### 2: load the prepared data and make last preparation #####

# set working directory and load prepared data 
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat')
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat')
setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")

load(file = 'data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData') # change to the current year (most recent year with prepared data)

###### 2.1: set YEAR and SPECIES and rm_point that should be analysed ####
# YEAR <- 2024 # set the most recent year
SPECIES # all species prepared data is available for 
SPECIES <- c('ACHU') # fill in SPECIES the analysis should be made for
rm_point<-c(99,76) # remove two points which are not independent

###### 2.2: check the data and remove unneeded stuff #####

head(countdata) # table with all birdcounts
head(obsCov) # all observation level Covs
rm(obsCov) # they are not needed since all information is available in countdata
head(siteCov) # site level Covs
(species_names <- species) # species codes and their English names
rm(species)

##### 3: prepare data for unmarkedMultFrame ####

###### 3.1: prepare siteCovs ####

head(siteCov)
siteCov <- siteCov %>% 
  rename_with(tolower) %>% 
  rename(canopy = canopy_cover, alt = elevation) %>% 
  filter(!point %in% rm_point) %>%  # remove all 2 points that are not independent
  arrange(point)

###### 3.2: prepare obsCov #####

obsCov <- countdata %>% 
  select(year, Point, Count, Rain, Wind, day, time, activity) %>% 
  rename_with(tolower) %>% 
  filter(!point %in% rm_point) %>% 
  arrange(point, year, count)

###### 3.3: prepare yearlySiteCovs ####

yearlySiteCovYear <- obsCov %>% 
  group_by(point,year) %>%
  summarise(season = as.numeric(mean(year) - 2010)) %>% # code year as.numeric
  spread(key = year, value = season) %>% # spread a key-value pair across multiple columns
  arrange(point)
yearlySiteCov <- list(year=yearlySiteCovYear[,2:ncol(yearlySiteCovYear)])

###### 3.3: prepare countdata ####

head(countdata)
countdata <- countdata %>% select(year, Point, Count, all_of(SPECIES)) %>% # select all the columns, including the correct species for analysis
  rename(n = all_of(SPECIES)) %>% # this filters for the one species that
  rename_with(tolower) %>% 
  filter(!point %in% rm_point) # remove the points that should not be analysed

occdata <- countdata %>% 
  mutate(occupancy = ifelse(n > 0, 1, 0)) %>% # convert in detection/non-detection
  mutate(season = paste(year, count, sep = '_')) %>% # connect year, count to string
  select(-n, -year, -count) %>% # remove unneeded columns 
  spread(key = season, value = occupancy) %>% # spread a key-value pair across multiple columns
  arrange(point) # sort rows in order of point number

###### 3.4: prepare numPrimary (number of survey years) #### 

numPrimary <- length(unique(countdata$year)) # calculates number of years

###### 3.6: check dimensions and input data.frames ####

head(siteCov)
head(obsCov)
head(yearlySiteCov)
head(occdata)

dim(occdata)
dim(siteCov)
nrow(obsCov)/(3*numPrimary)
dim(yearlySiteCovYear)

###### 3.3: create unmarkedMultFrame ####

umf <- unmarkedMultFrame(y = occdata[,2:ncol(occdata)], siteCovs = siteCov, yearlySiteCovs = yearlySiteCov, 
                         obsCovs = obsCov, numPrimary = numPrimary)
summary(umf) # looks good
str(umf) # looks good, but pay attention: not all variables from this data set are coded in an appropriate way (as.factor())

###### 3.4: scale numeric variables to fitting problems ####

siteCovs(umf)[c(2,4:6)] <- scale(siteCovs(umf)[c(2,4:6)]) # scale elevation, dbh and teeheight
obsCovs(umf)[c(1,4,6:8)] <- scale(obsCovs(umf)[c(1,4,6:8)]) # scale activity, rain, day and time (all numeric variables)
yearlySiteCovs(umf)[1] <- scale(yearlySiteCovs(umf)[1])
str(umf)
summary(umf)

##### 4: multi-season, dynamic single species occupancy model ####

###### 4.1: build global model and conduct MacKenzie-Bailey-Goodness-of-Fit Test ####

# built global model and perform gof
global_model <- colext(~alt:treeheight+dbh, ~year+alt, ~year+alt, ~day+time+I(time^2)+rain+wind+activity+location, data = umf, se = T) # global model which is a year corrected shift model
gof <- mb.gof.test(global_model, nsim = 100, parallel = T) # perform gof with 1000 sim
# saveRDS(best_model_gof, file = sprintf('output/data/GOF/%s_gof_mb.rds', SPECIES)) # save gof as RDS

fitstats <- function(global_model) {
  observed <- getY(global_model@data)
  expected <- fitted(global_model)
  resids <- residuals(global_model)
  sse <- sum(resids^2,na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

pb <- parboot(global_model, fitstats, nsim=20, report=1)  ### increase nsim
pb


# check gof
print(gof)
c_hat <- gof$c.hat.est # save c-hat estimates for usage of the QAIC
p_value = gof$p.value # save p_value for modsel table

###### 4.2: fit models for detection probability p() first for modSel ####

# fit models for detection probability p() manually, go on with fitList(), modSel()
fm1 <- colext(~1, ~1, ~1, ~1, data = umf, se = T)
fm2 <- colext(~1, ~1, ~1, ~day, data = umf, se = T)
fm3 <- colext(~1, ~1, ~1, ~day+time+I(time^2), data = umf, se = T)
fm4 <- colext(~1, ~1, ~1, ~day+time+I(time^2)+rain, data = umf, se = T)
fm5 <- colext(~1, ~1, ~1, ~day+time+rain+wind, data = umf, se = T)
fm6 <- colext(~1, ~1, ~1, ~day+time+I(time^2)+rain+wind+activity, data = umf, se = T)
fm6b <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm7 <- colext(~1, ~1, ~1, ~day+time+I(time^2)+rain+wind+activity+location, data = umf, se = T)

p_fitList <- list(fm1, fm2, fm3, fm4, fm5, fm6, fm7, fm6b)
names(p_fitList) <- lapply(p_fitList, function(x) formula(x)) # set formulas as model names 
(p_modSel_df <- aictab(cand.set = p_fitList, c.hat = c_hat) %>% # create a model comparison table with QAICc
  mutate(step = 'p'))
# best submodel for p(): day+time+rain+wind+activity, QAICc difference to second best 3.88 - go on with this fm6 best one

###### 4.2: fit models for initial occupancy psi() first for modSel ####

fm8 <- colext(~alt, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
# add treeheight
fm9 <- colext(~treeheight, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm10 <- colext(~alt+treeheight, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
# add dbh
fm11 <- colext(~dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm12 <- colext(~alt+dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm13 <- colext(~treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm14 <- colext(~alt+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
# add interaction between alt and treeheight
fm15 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm16 <- colext(~alt:treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)

# put the fitted models in a fitList() and rank them by QAICc in modSel()
psi_fitList <- list(fm6, fm8, fm9, fm10, fm11, fm12, fm13, fm14, fm15, fm16) # don't forget to include the best model from the last modeling step!
names(psi_fitList) <- lapply(psi_fitList, function(x) formula(x)) # set formulas as model names 
(psi_modSel_df <- aictab(cand.set = psi_fitList, c.hat = c_hat) %>% 
    mutate(step = 'psi'))
# best sub-model for psi(): ~1, QAICc difference to second best is 1.47 (~treeheight), then QAIC difference is 1.82 (~alt + treeheight) - go on with this best one (fm6)

###### 4.3: fit models for extinction and colonisation probability for modSel ####

fm17 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T) # constant model
fm18 <- colext(~1, ~alt, ~1, ~day+time+rain+wind+activity, data = umf, se = T) # expansion model
fm19 <- colext(~1, ~1, ~alt, ~day+time+rain+wind+activity, data = umf, se = T) # contraction model 
fm20 <- colext(~1, ~alt, ~alt, ~day+time+rain+wind+activity, data = umf, se = T) # shift model
fm21 <- colext(~1, ~year, ~year, ~day+time+rain+wind+activity, data = umf, se = T)  ## year model, this model will exclude the possibility that observed changes are just annual changes
fm22 <- colext(~alt+treeheight, ~year+alt, ~year, ~day+time+rain+wind+activity+location, data = umf, se = T)  # corrected year - expansion
fm23 <- colext(~alt+treeheight, ~year, ~year+alt, ~day+time+rain+wind+activity+location, data = umf, se = T)  # corrected year - contraction
fm24 <- colext(~alt+treeheight, ~year+alt, ~year+alt, ~day+time+rain+wind+activity+location, data = umf, se = T)  # corrected year - shift, also  global model

# put the fitted models in a list and rank them by QAIC in aictab
g_e_fitList <- list(constant = fm17, expansion = fm18, contraction = fm19,
                       shift = fm20, year = fm21, nullmodel = fm1, 
                       year_expansion = fm22, year_contraction = fm23, year_shift = fm24)
names(g_e_fitList) <- lapply(g_e_fitList, function(x) formula(x)) # set formulas as model names
model_names <- aictab(cand.set = list(constant = fm17, expansion = fm18, contraction = fm19, shift = fm20, year = fm21, nullmodel = fm1, year_expansion = fm22, year_contraction = fm23, year_shift = fm24), c.hat = c_hat)[, 1] # get model names in the correct order 
(g_e_modSel_df <- aictab(cand.set = g_e_fitList, c.hat = c_hat) %>%
  mutate(step = 'g_e',  # Add the step indicator
         model = model_names))  # include short model names
# best sub-model for g_e(): ~1, QAICc difference to second best expansion is 2.3 (~alt ~1), then QAIcC difference is 2.5 (~year ~year) year - go on with this best one (fm17)

##### 5: Explore best model and export the first things ####

###### 5.1: Best model ####
best_model <- fm17 # save best model as best_model
# saveRDS(best_model, file = sprintf('output/data/best_models/%s_best_model.rds', SPECIES)) # save model on local storage
summaryOD(best_model, c.hat = c_hat) # adjusted summary statistics with c-hat, there are NaNs for col!
names(best_model) # get names from the submodels

###### 5.2: Model selection tables ####
modSel_export <- rbind(as.data.frame(g_e_modSel_df) %>% rename(formula = Modnames) %>% select(model, formula, step, everything()), 
                       as.data.frame(psi_modSel_df) %>% rename(formula = Modnames) %>% mutate(model = NA) %>% select(model, formula, step, everything()), 
                       as.data.frame(p_modSel_df) %>% rename(formula = Modnames) %>% mutate(model = NA) %>% select(model, formula, step, everything()))
modSel_export <- modSel_export %>% mutate(species = SPECIES, p_value = p_value) %>%  # add species name and p_value from gof saved earlier
  select(species, everything()) %>% # change order of columns
  arrange(match(step, c('g_e','psi','p')), QAICc) # check for correct order
# write.csv(modSel_export, file = sprintf('output/data/modSel/%s_modSel_full.csv', SPECIES)) # export as .csv file

##### 6: draw inference from the best model ####

###### 6.1: Extract Occupancy values for population trajectory ####

# using empirical bayes estimates of occupancy averaged across all points
occupancy <- ranef(best_model)

# manipulate data structure of occupancy estimates, summarize and store data in data.frame
occupancy_data <- as.data.frame(bup(occupancy, stat = 'mean')) %>% # posterior mean by stat = ''
  gather(key = 'Year',value = 'occu') %>% # transfer from wide to long format 
  mutate(Year = as.numeric(str_replace(Year,'V',''))) %>% # delete V before number of season
  mutate(Year = Year + 2010) %>% # calculate year 
  group_by(Year) %>%
  summarise(Occupancy = mean(occu))

# calculate confidence intervals and store them in occupancy_data data.frame for plotting
occupancy_confint <- confint(occupancy, level = 0.95) # 95% CI
occupancy_data$lower_cl <- apply(occupancy_confint[,1,],2, mean)
occupancy_data$upper_cl <- apply(occupancy_confint[,2,],2, mean)

# remove data from 2020 because in this year surveys didn't take place 
occupancy_data <- occupancy_data %>% filter(!Year == '2020')

# make a quick plot out of curiosity
occupancy_data %>% ggplot() + 
  geom_point(mapping = aes(x = Year, y = Occupancy), size = 3, col = 'firebrick') + 
  geom_errorbar(aes(ymin = lower_cl,ymax = upper_cl, x = Year), width = 0.2) + 
  scale_y_continuous(limits = c(0, 1.05), breaks=seq(0,1,0.2),labels=seq(0,1,0.2))  +
  scale_x_continuous(limits = c(2011, 2024), breaks=seq(2011,2023,2),labels=seq(2011,2023,2)) + 
  labs(y="Mean occupancy", title= paste0(SPECIES, ' Occupancy Trajectory')) + 
  theme_bw()

###################################################
# try other ways to calculate yearly occupancy 
##################################################

# nonparametric bootstrap to obtain SEs and then calculate inflated CIs 
best_model <- nonparboot(best_model, B = 300) # perform nonparametric bootstrap, higher B would be better (<1000)

# smoothed 
occu_sm <- smoothed(best_model)[2,] # finite sample quantity, proportion of sampled sites that are occupied (not the entire population of sites)
occu_sm_l <- occu_sm - (1.96*best_model@smoothed.mean.bsse[2,]) # calculate CI from SE by taking estimates and calculate +/- 2 times SE (and Merc wrote something about backTransform/plogis, but I didn't quite understand that)
occu_sm_u <- occu_sm + (1.96*best_model@smoothed.mean.bsse[2,])

# projected - population-level estimates of occupancy probability
occu_proj <- projected(best_model)[2,]
occu_proj_l <- occu_proj - (1.96*best_model@projected.mean.bsse[2,]) # calculate CI from SE by taking estimates and calculate +/- 2 times SE (and Marc wrote something about backTransform/plogis, but I didn't quite understand that)
occu_proj_u <- occu_proj + (1.96*best_model@projected.mean.bsse[2,])


# make a plot to cisualise the difference between the methods
occupancy_data %>% mutate(Type = 'ranef()') %>% 
  rbind(data.frame(Year = 2011:2024, Occupancy = occu_sm, upper_cl = occu_sm_l, lower_cl = occu_sm_u, Type = 'smoothed() +/- 1.96SE')) %>% 
  rbind(data.frame(Year = 2011:2024, Occupancy = occu_proj, upper_cl = occu_proj_l, lower_cl = occu_proj_u, Type = 'projected() +/- 1.96SE')) %>% 
  filter(!Year == 2020) %>% 
  ggplot() + 
  geom_line(mapping = aes(x = Year, y = Occupancy, colour = Type)) +
  geom_point(mapping = aes(x = Year, y = Occupancy, colour = Type)) +
  geom_ribbon(mapping = aes(x = Year, ymin = lower_cl, ymax = upper_cl, fill = Type), alpha = 0.2) + 
  theme_bw() + 
  scale_fill_viridis_d(alpha=0.1,begin=0,end=1,direction=-1) + # set begin/end to 0.98 for yellow/ext // 0 for purple/col
  scale_color_viridis_d(alpha=1,begin=0,end=1,direction=-1) + 
  labs(title = paste0(SPECIES, ' Comparison Mean Occupancy Extraction Methods'))

# save occupancy_data for later plotting
# fwrite(occupancy_data, file = sprintf('output/data/ranef_occupancy_data/%s_occupancy_data.csv', SPECIES))

###### 6.2: Make predictions for colonisation and extinction if elevation is included as predictor in the best model ####

# create input df with for prediction, 
nd <- data.frame(day = 0, time = 0, rain = 0, wind = 1, activity = max(umf@obsCovs$activity, na.rm = T), # use maximum bird activity, lowest wind speed, mean of time and day = 0 
                 location = 'midslope', treeheight = 0, # location with valley or midslope used, mean scaled treeheight used, should be 0 
                 dbh = 0, canopy = 0, # mean scaled dbh and canopy used, should be 0 
                 year = 0, # mean of scaled year, should be 0
                 alt = seq(from = min(umf@siteCovs$alt), # alt: scaled altitude for prediction,
                           to = max(umf@siteCovs$alt), by = 0.02),
                 elevation = rescale(seq(from = min(umf@siteCovs$alt), to = max(umf@siteCovs$alt), by = 0.02), # elevation: rescaled altitude for plotting and an exact ecological meaning
                                    to = c(min(siteCov$alt), max(siteCov$alt)), # to = output range as c()
                                    from = c(min(umf@siteCovs$alt), max(umf@siteCovs$alt)))) # from = input range as c()

# predict values for col and ext using modavgPred() on a list that just contains the best_modl, input data from nd
pred_col <- as.data.frame(modavgPred(cand.set = list(best_model), newdata = nd, c.hat = c_hat, parm.type = 'gamma')) %>%  # for parm.type choose one of 'psi', 'gamma', 'epsilon', 'detect'
  mutate(Type = 'Colonisation', Elevation = nd$elevation, Species = SPECIES) %>% # add Elevation and Species for easier plotting
  rename(Predicted = mod.avg.pred, SE = uncond.se, lower = lower.CL, upper = upper.CL) %>% 
  select(Predicted, SE, lower, upper, Type, Elevation, Species) # select only the needed columns 

pred_ext <- as.data.frame(modavgPred(cand.set = list(best_model), newdata = nd, c.hat = c_hat, parm.type = 'epsilon')) %>%  # for parm.type choose one of 'psi', 'gamma', 'epsilon', 'detect'
  mutate(Type = 'Extinction', Elevation = nd$elevation, Species = SPECIES) %>% # add Elevation and Species for easier plotting 
  rename(Predicted = mod.avg.pred, SE = uncond.se, lower = lower.CL, upper = upper.CL) %>% 
  select(Predicted, SE, lower, upper, Type, Elevation, Species) # select only the needed columns 

pred_ext_col <- bind_rows(pred_ext, pred_col) # connect tables 
# fwrite(pred_ext_col, file = sprintf('output/data/pred_col_ext/%s_pred_ext_col.csv', SPECIES)) # only export pred_col

# quickly plot col-ext dynamics against elevation 
pred_ext_col %>%
  ggplot(aes(x = Elevation, y = Predicted, colour = Type, fill = Type)) +
  geom_line(linewidth = 2) +
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.2) +
  labs(x="Elevation", y="Predicted probability", title=paste0(SPECIES, ' Elevational Range Dynamics')) +
  theme_bw()

###### 6.3: Extract all coefficients ####

# this calculates all estimates with inflated SE and CI using c-hat 
summary(best_model) # uninflated effect sizes for comparison
best_model_summary <- summaryOD(best_model, c.hat = c_hat) # summary statistics corrected for overdispersion by c-hat on logit-scale
(estimates_best_model <- as.data.frame(best_model_summary$outMat) %>% 
    mutate(species = SPECIES) %>% 
    rename(SE = se, lower = lowlim, upper = upplim))
summary(best_model) # uninflated effect sizes for comparison

# export coefficient table 
# fwrite(estimates_best_model, sprintf('output/data/best_model_output_estimates/%s_best_model_output_estimates_all.csv', SPECIES))

#########################################################################################################################
########################################################################################################################



# save occupancy_data for later plotting
fwrite(occupancy_data, file = sprintf('output/data/ranef_occupancy_data/%s_occupancy_data.csv', SPECIES))

##### 8: export some prepared data #### 

saveRDS(umf, file = sprintf('output/data/prepared_data/%s_unmarked_mult_frame.rds', SPECIES))

##### 9: notes ####

