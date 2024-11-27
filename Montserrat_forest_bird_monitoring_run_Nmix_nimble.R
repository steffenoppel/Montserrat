#### Scrip to run Nmix models for the forest bird monitoring Montserrat in NIMBLE across all species #####
#### this script is part of the automated workflow for an annual report for the Centre Hills Forest Bird Monitoring in Montserrat ####
#### Script assembled by Filibert Heim, filibert.heim@posteo.de, in Nov 2024 but mainly takan and from a script from Steffen Oppel 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages 
library(tidyverse)
library(data.table) 
library(lubridate)
library(MCMCvis)
library(nimble)
library(parallel)
library(doParallel)

# load data which has been prepared in the script 'Montserrat_forest_bird_monitoring_data_prep_for_Nmix_nimble.R'
load(file = 'data/Montserrat_forest_bird_monitoring_yearly_NIMBLE_model_data.RData')

# set species s for which the model should run 
s <- commandArgs(trailingOnly = TRUE)[1] # this runs only in the github actions workflow  

# set clusters that should be used in a parallel manner 
cl<-makeCluster(4)
registerDoParallel(cl)

# print a message for the github actions workflow 
cat("Running Nmixture models in NIMBLE for the species", s, 'with', n.iter, 'iterations after', n.burnin, 'iterations in the burnin for overall', n.chains, 'chains.', "\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Start to run model for one SPECIES  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Take subset of data for focal species and sort the tables  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
bird_s<-SURVEYDATA[,c(1,2,3,4,match(s,colnames(SURVEYDATA)))] %>%
  arrange(Point,year,Count) %>%
  rename(N=5) %>%
  dplyr::select(Point,year,Count,N)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Create bird data input matrix  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Create input data for nimble  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

#### DISTINGUISH CONSTANTS AND DATA
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~
  
trend.data <- list(M = BIRD.y)
  
  
####   ADD INITIAL VALUES----     ################################
## MUST ADD Nst TO INITIAL VALUESBE FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

inits.trend$p = array(runif(trend.constants$nsite*trend.constants$nrep*trend.constants$nyear, 0.3,0.7),
                        dim= c(trend.constants$nsite, trend.constants$nrep,trend.constants$nyear)) 
inits.trend$N = Nst
inits.trend$M = inits.y
inits.trend$M.new = inits.new

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Finally run the model in NIMBLE  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

### this takes 3-5 hrs on Steffens laptop for 250000 iterations and converges for most species
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
  
# save Trendmodel to output
saveRDS(TRENDMOD,sprintf("output/%s_trend_model_nimble.rds",s))
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Evaluate model fit with bayesian p-value  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

### COMBINE SAMPLES ACROSS CHAINS
MCMCout<-as_tibble(rbind(TRENDMOD$samples[[1]],TRENDMOD$samples[[2]],TRENDMOD$samples[[3]]))
  
### PLOT AND SAVE GoF PLOT
# ylow<-round((min(MCMCout$fit,MCMCout$fit.new, na.rm=T)-50)/1000,1)*1000
# yup<-round((max(MCMCout$fit,MCMCout$fit.new, na.rm=T)+50)/1000,1)*1000
# pdf(sprintf("output/%s_trendmodel_fit2024.pdf",s), width=10, height=10, title="")
# plot(MCMCout$fit, MCMCout$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(ylow, yup), ylim = c(ylow, yup))
# abline(0, 1, lwd = 2, col = "black")
# dev.off()
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Examine output and diagnostics with MCMCvis  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out<- as.data.frame(MCMCsummary(TRENDMOD$samples, params=c("trend","trend2","totalN","anndet")))
out$parameter<-row.names(out)
out$species<-s
out$BayesP<-mean(MCMCout$fit > MCMCout$fit.new)
out$GoFSlope<-mean(MCMCout$fit) / mean(MCMCout$fit.new)
names(out)[c(3,4,5)]<-c('lcl','median', 'ucl')
fwrite(out,sprintf("output/%s_trend_estimates2024.csv",s))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8. Export model object and trend estimates  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# save model to output
saveRDS(TRENDMOD,sprintf("output/%s_trend_model_nimble.rds",s))

# save trendestimates to output
fwrite(out,sprintf("output/%s_trend_estimates2024.csv",s))

