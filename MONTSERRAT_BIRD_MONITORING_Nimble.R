######################################################################################
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  ANALYSIS OF ANNUAL SURVEYS   ##########################################
#############  steffen.oppel@gmail.com      ##########################################
######################################################################################

### based on multi-year trend model of Kery et al. 2009
### REQUIRES DATA PREPARATION FROM "MONTSERRAT_BIRD_MONITORING_data_prep.R"
### updated in 2024 to use Nimble and parallelize loop over species

######################################################################################
#############  Load required packages       ##########################################
######################################################################################

library(tidyverse)
library(data.table)
library(dplyr)
library(dtplyr)
library(lubridate)
library(ggplot2)
# library(knitr)
# library(rmarkdown)
library(MCMCvis)
library(nimble)
library(basicMCMCplots) # for trace plots called chainsPlot
library(parallel)
library(foreach)
library(doParallel)

######################################################################################
#############  Set your working directory (path where the database is)       #########
######################################################################################

#setwd("C:\\STEFFEN\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")
setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")
#setwd("C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat\\Montserrat")



######################################################################################
#############  load the pre-prepared dataset					     #########
######################################################################################

load("data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData")
#load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_ANNUAL_DATA_INPUT2023.RData")
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
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
  loglam~dunif(-5,5)          ##  mean abundance prior
  trend~dunif(-5,5)         ##  trend prior
  beta.elev~dunif(-2,2)
  #beta.rain~dunif(-2,2)
  beta.canopy~dunif(-2,2)
  beta.treeheight~dunif(-2,2)
  bwind~dunif(-2,0)   ## wind can only have negative effect on detection
  brain~dunif(-2,0)   ## rain can only have negative effect on detection
  btime~dunif(-2,2)
  bday~dunif(-2,2)
  bridge~dunif(-2,2)
  bact~dunif(-2,2)
  
  ## SITE RANDOM EFFECT ##
  for(i in 1:nsite){
    lam.site[i]~dnorm(loglam,tau=tau.site)    ## site-specific random effect with hierarchical centering from Kery email 5 June 2018
  }
  tau.site<-1/(sigma.site*sigma.site)
  sigma.site~dunif(0,2)
  
  ## YEAR RANDOM EFFECT FOR ABUNDANCE AND ANNUALLY VARYING DETECTION PROBABILITY ##
  for(year in 1:nyear){
    p0[year]~dunif(0.01,0.99)## detection probability
    logitp0[year]<-log(p0[year]/(1-p0[year]))
    lam.year[year]~dnorm(trend*primocc[year],tau=tau.year)    ## year-specific random effect with hierarchical centering from Kery email 5 June 2018
  }
  tau.lp<-1/(sigma.p*sigma.p)
  sigma.p~dunif(0,2)
  tau.year<-1/(sigma.year*sigma.year)
  sigma.year~dunif(0,2)
  
  
  ######### State and observation models ##############
  for(year in 1:nyear){
    for(i in 1:nsite){
      log(lambda[i,year])<- lam.year[year]+
        #beta.rain*rain[year]+
        beta.elev*elev[i]+
        beta.treeheight*treeheight[i]+
        beta.canopy*canopy[i]+
        lam.site[i]
      N[i,year]~dpois(lambda[i,year])
      
      for(t in 1:nrep){
        M[i,t,year]~dbin(p[i,t,year],N[i,year])
        p[i,t,year] <- exp(lp[i,t,year])/(1+exp(lp[i,t,year]))
        lp[i,t,year] ~ dnorm(mu.lp[i,t,year], tau=tau.lp)
        mu.lp[i,t,year]<-logitp0[year] +
          btime*time[i,t,year]+
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
                        primocc=seq(2011:YEAR),
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
                    trend=runif(1,-2,2),
                    loglam = runif(1,-2,2),
                    sigma.site = runif(1,0,2),
                    sigma.year=runif(1,0,2),
                    sigma.p=runif(1,0,2),
                    beta.canopy=runif(1,-2,2),
                    #beta.rain=runif(1,-2,2),
                    beta.treeheight=runif(1,-2,2),
                    beta.elev=runif(1,-2,2),
                    bwind=-1,
                    brain=-1,
                    bridge=-1,
                    btime=-1,
                    bday=1,
                    bact=2,
                    p0 = runif(nyears,0.1,0.9))
inits.trend$lam.site<-rnorm(nsites,inits.trend$loglam,inits.trend$sigma.site)
inits.trend$lam.year<-rnorm(nyears,(inits.trend$trend*seq(1:(nyears))),inits.trend$sigma.year)





####   DEFINE RUN SETTINGS AND OUTPUT DATA----     ################################

# Define parameters to be monitored
parameters.trend <- c("fit", "fit.new","trend","totalN","anndet")  #


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 250000
n.burnin <- 150000
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
inits.trend$lp = array(rnorm(trend.constants$nsite*trend.constants$nrep*trend.constants$nyear, c(test$mu.lp), inits.trend$sigma.p),
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
# #s="MTOR"

######################################################################################
#############  TAKE SUBSET OF DATA FOR FOCAL SPECIES AND SORT THE TABLES    ###################
######################################################################################

bird_s<-SURVEYDATA[,c(1,2,3,4,match(s,colnames(SURVEYDATA)))] %>%
  arrange(Point,year,Count) %>%
  rename(N=5) %>%
  mutate(N=if_else(is.na(VisitID),NA,N)) %>%  ### RE-INTRODUCE THE NAs for COUNTS THAT DID NOT TAKE PLACE #####
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
    mutate(N=ifelse(is.na(N),median(bird_s$N, na.rm=T),NA)) %>%   ### fill in missing values
    dplyr::filter(year==y) %>%
    dplyr::select(Point, Count, N) %>%
    tidyr::spread(key=Count, value=N) %>%
    dplyr::arrange(Point)
  inits.y[,,yc]<-as.matrix(x[,2:4])
  
  x<-bird_s %>%
    mutate(N=ifelse(is.na(N),median(bird_s$N, na.rm=T),N)) %>%   ### fill in missing values
    dplyr::filter(year==y) %>%
    dplyr::select(Point, Count, N) %>%
    tidyr::spread(key=Count, value=N) %>%
    dplyr::arrange(Point)
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
inits.trend$lp = array(rnorm(trend.constants$nsite*trend.constants$nrep*trend.constants$nyear, c(test$mu.lp), inits.trend$sigma.p),
                       dim= c(trend.constants$nsite, trend.constants$nrep,trend.constants$nyear))
inits.trend$N = Nst
inits.trend$M = inits.y
inits.trend$M.new = inits.new

allchaininits.trend <- list(inits.trend, inits.trend, inits.trend)




###############################################################################
####   RUN THE MODEL IN NIMBLE  --------------------###########################
###############################################################################


### this takes 14-15 mins for 50000 iterations and converges in that time
TRENDMOD <- nimbleMCMC(code = trend.model,
                            constants=trend.constants,
                            data = trend.data,
                            inits = allchaininits.trend,
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
ylow<-round((min(MCMCout$fit)-50)/1000,1)*1000
yup<-round((max(MCMCout$fit)+50)/1000,1)*1000
pdf(sprintf("output/%s_trendmodel_fit2024.pdf",s), width=10, height=10, title="")
plot(MCMCout$fit, MCMCout$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(ylow, yup), ylim = c(ylow, yup))
abline(0, 1, lwd = 2, col = "black")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMINE OUTPUT AND DIAGNOSTICS WITH MCMCvis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out<- as.data.frame(MCMCsummary(TRENDMOD$samples, params=c("trend","totalN","anndet")))
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



write.table(annestimates,"output/Annual_estimates2024.csv", row.names=F, sep=",")
write.table(trendout,"output/Trend_estimates2024.csv", row.names=F, sep=",")




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
  mutate(col=ifelse(lower95CI<0,ifelse(upper95CI<0,"darkred","black"),ifelse(upper95CI>0,"forestgreen","black"))) %>%
  mutate(col=ifelse(species=="CAEL","darkred",col))
annestimates %>% filter(Year!=2020) %>%
  mutate(col = as.factor(trendout$col[match(species,trendout$species)])) %>%


ggplot()+
geom_line(aes(x=Year, y=trend,col=col), size=1)+facet_wrap(~fullspec, ncol=2, scales="free_y")+
geom_point(aes(x=Year, y=trend,col=col), size=2)+
#geom_ribbon(data=annestimates,aes(x=Year, ymin=lower95CI,ymax=upper95CI),alpha=0.2)+
geom_errorbar(aes(x=Year, ymin=lower95CI,ymax=upper95CI,col=col), width=.1) +

## remove the legend
theme(legend.position="none")+
guides(fill=FALSE)+
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


surveys2023<-obsCov %>% filter(year==2023)
birds2023<-birds %>% filter(Year==2023) %>% rename(N=SumOfNumber1) %>% mutate(Point=as.integer(as.character(Point))) %>% filter(!Species %in% c('UNK','NA'))
summary2023<-surveys2023 %>% select(Year, Point, Count, Date) %>%
  left_join(birds2023, by=c('Year','Point','Count')) %>%
  group_by(Count, Species) %>%
  summarise(N=sum(N, na.rm=T))

totals2023<-summary2023 %>% group_by(Count) %>%
  summarise(N=sum(N), n_spec=length(unique(Species)))

table1<-summary2023 %>% spread(key=Count, value=N, fill = 0) %>%
  filter(!is.na(Species)) %>%
  mutate(Species=species$Species[match(Species,species$SpeciesCode)]) %>%
  arrange(desc(`1`))

fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")
# summary<-COUNTDATA %>% group_by(Species,Year) %>%
#   summarise(mean=mean(N, na.rm=T), sd=sd(N, na.rm=T)) %>%
#   filter(Year!=2020) %>% ### no counts were done in 2020
#   mutate(Species=fullnames[match(Species,SPECIES)])
# summary


table2<-trendout %>%
  mutate(dir=ifelse(lower95CI<0,ifelse(upper95CI<0,"decrease","stable"),ifelse(upper95CI>0,"increase","stable"))) %>%
  mutate(dir=ifelse(species=="CAEL","(decrease)",dir)) %>%
  mutate(conf=paste(trend," (",lower95CI," - ",upper95CI,")", sep="")) %>%
    select(fullspec,dir,conf,pval)

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









