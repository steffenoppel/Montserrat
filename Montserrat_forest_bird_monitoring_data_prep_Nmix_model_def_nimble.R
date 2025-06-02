#### Scrip for data preparation for yearly Nmix models in NIMBLE #####
#### this script is part of the automated workflow for an annual report for the Centre Hills Forest Bird Monitoring in Montserrat ####
#### Script assembled by Filibert Heim, filibert.heim@posteo.de, in Nov 2024 but mainly taken and from a script from Steffen Oppel 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages 
library(tidyverse)
library(data.table)
library(nimble)

# load data from general data preparation script 'Montserrat_forest_bird_monitoring_data_prep.R'
load(file = 'data/MONTSERRAT_ANNUAL_DATA_INPUT.RData')

# save SPECIES as full names  
SPECIES # print species codes 
fullnames <- c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler", # save full species names in the same order as SPECIES 
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")

# define dimensions of arrays 
nsites <- length(unique(siteCov$Point))
nyears<-length(unique(countdata$year))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Create site covariate data input matrix --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

siteCov<-siteCov %>% arrange(Point) %>%
  mutate(ridge=ifelse(Location=="Ridge",1,0)) %>%
  dplyr::select(Point,treeheight,Elevation,Canopy_cover,ridge) %>%
  mutate(tree=scale(treeheight)[,1], elev=scale(Elevation)[,1],canopy=scale(Canopy_cover)[,1])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Create observation covariate input matrix --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SURVEYDATA<-countdata %>%
  arrange(Point,year,Count) %>%
  mutate(activity=ifelse(is.na(activity),mean(activity, na.rm=T),activity)) %>%
  mutate(time=scale(time),day=scale(day),activity=scale(activity))

### create array for each covariate

wind<-array(NA, dim=c(nsites,3,nyears))
rain<-array(NA, dim=c(nsites,3,nyears))
time<-array(NA, dim=c(nsites,3,nyears))
day<-array(NA, dim=c(nsites,3,nyears))
ACT<-array(NA, dim=c(nsites,3,nyears))				## REPLACED ON 2 MAY WITH RAINFALL AMOUNT

# fill in array for each covariate
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. replace all NA's in covariates otherwise 'undefinded node' error --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (d in 1:nyears){	# replace missing dates with mean for each survey round in each year
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Write the NIMBLE model and set inits  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Create input data fro NIMBLE - generic part  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Define run settings and output data  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define parameters to be monitored
parameters.trend <- c("fit", "fit.new","trend","trend2","totalN","anndet")  #


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 150000   #150000
n.burnin <- 100000  #100000
n.chains <- 3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8. Preliminary test of NIMBLE model to identify problems --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9. Export workspace for running the models in extra scripts  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# fill in the objects you want to save from the environment for further analysis
# needed_objects <- c() #### FILL IN ALL OBJECTS WHICH ARE NEEDED FOR THE NEXT STEP OF RUNNING THE MODELS FOR EACH SPECIES 
# rm(list = setdiff(ls(), needed_objects)) # remove all unneeded objects 

# save the prepared data for other possible analysis 
save.image('data/Montserrat_forest_bird_monitoring_yearly_NIMBLE_model_data.RData')








