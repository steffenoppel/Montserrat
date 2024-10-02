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
library(dtplyr)
library(lubridate)
library(ggplot2)
library(knitr)
library(rmarkdown)
library(MCMCvis)
library(nimble)
library(basicMCMCplots) # for trace plots called chainsPlot
library(parallel)


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

### summarise rainfall from Jan to March, productivity from PREVIOUS year will affect count in current year
rain<-fread("data/MontserratRain2005_2023.csv",fill=TRUE) %>%
  dplyr::filter(Variable=="RainMM") %>%
  dplyr::filter(YEAR %in% seq(2010,2024,1)) %>%
  dplyr::select(-Variable,-Total) %>%
  gather(key="Month", value="mm",-YEAR) %>%
  dplyr::filter(Month %in% c('JAN','FEB','MAR')) %>%
  group_by(YEAR) %>%
  summarise(rain=sum(mm)) %>%
  mutate(rain=scale(rain)[,1])

###############################################################################
############## CREATE SITE COVARIATE DATA INPUT MATRIX   ######################
###############################################################################
nsites<-length(unique(siteCov$Point))

siteCov<-siteCov %>% arrange(Point) %>%
  dplyr::select(Point,treeheight,Elevation,Canopy_cover,ridge) %>%
  mutate(tree=scale(treeheight)[,1], elev=scale(Elevation)[,1],canopy=scale(Canopy_cover)[,1])





###############################################################################
############## CREATE OBSERVATION COVARIATE DATA INPUT MATRIX   ###############
###############################################################################
SURVEYDATA<-SURVEYDATA %>%
  arrange(Point,Year,Count) %>%
  mutate(time=scale(time),Day=scale(Day))

## SORT THE TABLE SO IT HAS THE SAME ORDER AS THE BIRD DATA
obsCov<-obsCov %>%
  arrange(Point,Year,Count)


### only needs standardisation if measured in mm, not as 0/1 variable
#meant<-mean(SURVEYDATA$rain, na.rm = TRUE)
#sdt<-sd(SURVEYDATA$rain, na.rm = TRUE)
#SURVEYDATA$rain<-(SURVEYDATA$rain-meant)/sdt


### create array for each covariate

wind<-array(NA, dim=c(nsites,3,nyears))
time<-array(NA, dim=c(nsites,3,nyears))
ACT<-array(NA, dim=c(nsites,3,nyears))				## REPLACED ON 2 MAY WITH RAINFALL AMOUNT


### fill in array for each covariate
for (y in 2011:YEAR){
  obsC<-subset(SURVEYDATA, Year==y)
  y<-match(y,c(2011:YEAR))						## translates the year (2011, 2012, etc.) into consecutive number (1,2,...) for array dimensions
  x<-obsC %>% dplyr::select(Point, Count, time) %>% tidyr::spread(key=Count, value=time) %>% dplyr::arrange(Point)
  time[,,y]<-as.matrix(x[,2:4])
  
  x<-obsC %>% dplyr::select(Point, Count, wind) %>% tidyr::spread(key=Count, value=wind) %>% dplyr::arrange(Point)
  wind[,,y]<-as.matrix(x[,2:4])
  
  x<-obsC %>% dplyr::select(Point, Count, ACT) %>% tidyr::spread(key=Count, value=ACT) %>% dplyr::arrange(Point)
  ACT[,,y]<-as.matrix(x[,2:4])

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
  wind[is.na(wind[,1,d]),1,d]<-0
  wind[is.na(wind[,2,d]),2,d]<-0
  wind[is.na(wind[,3,d]),3,d]<-0
}






######################################################################################
#############  WRITE THE NIMBLE MODEL AND SET INITS  ############################
######################################################################################

# Specify model in NIMBLE format
trend.model<-nimbleCode({
  
  
  ####  Priors ########
  loglam~dunif(-5,5)          ##  mean abundance prior
  trend~dunif(-10,10)         ##  trend prior
  beta.elev~dunif(-5,5)
  beta.canopy~dunif(-5,5)
  beta.treeheight~dunif(-5,5)
  bwind~dunif(-5,5)
  btime~dunif(-5,5)
  bridge~dunif(-5,5)
  bact~dunif(-5,5)
  
  ## SITE RANDOM EFFECT ##
  for(i in 1:nsite){
    lam.site[i]~dnorm(loglam,tau=tau.site)    ## site-specific random effect with hierarchical centering from Kery email 5 June 2018
  }
  tau.site<-1/(sigma.site*sigma.site)
  sigma.site~dunif(0,10)
  
  ## YEAR RANDOM EFFECT FOR ABUNDANCE AND ANNUALLY VARYING DETECTION PROBABILITY ##
  for(year in 1:nyear){
    p0[year]~dunif(0,1)## detection probability
    logitp0[year]<-log(p0[year]/(1-p0[year]))
    lam.year[year]~dnorm(trend*primocc[year],tau=tau.year)    ## year-specific random effect with hierarchical centering from Kery email 5 June 2018
  }
  tau.lp<-1/(sigma.p*sigma.p)
  sigma.p~dunif(0,10)
  tau.year<-1/(sigma.year*sigma.year)
  sigma.year~dunif(0,10)
  
  
  ######### State and observation models ##############
  for(year in 1:nyear){
    for(i in 1:nsite){
      log(lambda[i,year])<- lam.year[year]+beta.rain*rain[year]+beta.elev*elev[i]+beta.treeheight*treeheight[i]+beta.canopy*canopy[i]+lam.site[i]
      N[i,year]~dpois(lambda[i,year])
      
      for(t in 1:nrep){
        M[i,t,year]~dbin(p[i,t,year],N[i,year])
        p[i,t,year] <- exp(lp[i,t,year])/(1+exp(lp[i,t,year]))
        lp[i,t,year] ~ dnorm(mu.lp[i,t,year], tau=tau.lp)
        mu.lp[i,t,year]<-logitp0[year] + btime*time[i,t,year]+ bridge*ridge[i]+ bwind*wind[i,t,year]+ bact*ACT[i,t,year]
        
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
  fit <- sum(E2[1:nsite,1:nrep,1:nyear])# Sum up squared residuals for actual data set
  fit.new <- sum(E2.new[1:nsite,1:nrep,1:nyear]) # Sum up for replicate data sets
}) ## end of nimble code chunk



# 5.1. Define initial values

# Initial values for some parameters
## MUST BE FOR ALL PARAMETERS
## NIMBLE CAN HAVE CONVERGENCE PROBLEMS IF DIFFERENT INITS ARE SPECIFIED: https://groups.google.com/g/nimble-users/c/dgx9ajOniG8

inits.trend <- list(N = Nst,
                    trend=runif(1,-2,2),
                    loglam = runif(1,-2,2),
                    sigma.site = runif(1,0,2),
                    sigma.year=runif(1,-2,2),
                    sigma.p=runif(1,-2,2),
                    beta.canopy=runif(1,-2,2),
                    beta.rain=runif(1,-2,2),
                    beta.treeheight=runif(1,-2,2),
                    beta.elev=runif(1,-2,2),
                    bwind=runif(1,-2,2),
                    bridge=runif(1,-2,2),
                    btime=runif(1,-2,2),
                    bact=runif(1,-2,2),
                    p0 = runif(nyears,0,1))



# 5.3. Define parameters to be monitored
parameters.trend <- c("trend","totalN","fit", "fit.new","anndet")


# MCMC settings
# number of posterior samples per chain is n.iter - n.burnin
n.iter <- 50000
n.burnin <- 25000
n.chains <- 4
inits.trend <- list(inits.trend, inits.trend, inits.trend,inits.trend)


# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = trend.model,
                    constants=trend.constants,
                    data = trend.data,
                    inits = inits.trend,
                    calculate=TRUE)

test$initializeInfo()

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
test$calculate()
test$calculate(nodes="lam.site") # this causes a NaN when the nest probability (z[i] = 0)
test$logProb_lam.site
test$logProb_lam.year
test$logProb_tau.lp
test$logProb_tau.year
test$logProb_M
test$logProb_lambda



test$logProbs$anndet
# Missing values (NAs) or non-finite values were found in model variables:
# lam.site, lam.year, tau.lp, tau.year,
# lambda, M, p, lp, mu.lp, anndet, eval, sd.resi,
# E, E2, M.new, E.new, E2.new, fit, fit.new.

test$initializeInfo()
#help(modelInitialization)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
configureMCMC(test) # check that the samplers used are ok - all RW samplers need proper inits


  










######################################################################################
#############  START THE LOOP OVER EVERY SPECIES          ############################
######################################################################################
trendout<-data.frame(species=SPECIES, timeframe=sprintf("2011-%i",YEAR), trend=0, lower95CI=-1, upper95CI=1, pval=0, slope=1)
annestimates<-data.frame(species=rep(SPECIES, each=nyears), Year=seq(2011,YEAR), trend=0, lower95CI=-1, upper95CI=1, detprob=0, detproblower95CI=-1, detprobupper95CI=1) 


for (s in SPECIES){


######################################################################################
#############  TAKE SUBSET OF DATA FOR FOCAL SPECIES AND SORT THE TABLES    ###################
######################################################################################

bird_s<-subset(COUNTDATA, Species==s)
bird_s<-bird_s[order(bird_s$Point,bird_s$Year, bird_s$Count, decreasing=F),] 


### RE-INTRODUCE THE NAs for COUNTS THAT DID NOT TAKE PLACE #####
bird_s$N[is.na(SURVEYDATA$time)]<-NA


###############################################################################
############## CREATE BIRD DATA INPUT MATRIX   ################################
###############################################################################

### create array to be filled with data
BIRD.y<-array(NA, dim=c(nsites,3,nyears))

#### GET THE MAXIMUM COUNT PER POINT PER YEAR FOR INITIAL VALUES
Nst<-as.matrix(bird_s %>%
                 mutate(N=ifelse(is.na(N),median(bird_s$N, na.rm=T),N)) %>%   ### fill in missing values - switch to max if there is invalid parent error
                 group_by(Point, Year) %>%
                 summarise(K=max(N, na.rm=T)) %>%
                 spread(key=Year,value=K, fill=max(bird_s$N,na.rm=T)) %>%
                 ungroup() %>%
                 arrange(Point) %>%
                 dplyr::select(-Point))




######################################################################################################
########## CREATE INPUT DATA FOR NIMBLE ------------------------
#######################################################################################################

#### DISTINGUISH CONSTANTS AND DATA
# Constants are values that do not change, e.g. vectors of known index values or the indices used to define for loops
# Data are values that you might want to change, basically anything that only appears on the left of a ~
R = nrow(BIRD.y)
T = ncol(BIRD.y)
nyears = dim(BIRD.y)[3]
trend.constants <- list(nsite=dim(BIRD.y)[1],
                            nrep=ncol(BIRD.y),
                            primocc=seq(1:(dim(BIRD.y)[3])),
                            nyear=dim(BIRD.y)[3],
                            elev=siteCov$elev,
                            treeheight=siteCov$tree,
                            canopy=siteCov$canopy,
                            rain=rain$rain,
                            wind=wind,
                            ridge=siteCov$ridge,
                            time=time,
                            ACT=ACT)


trend.data <- list(M = BIRD.y)






###############################################################################
####   DEFINE RUN SETTINGS AND OUTPUT DATA     ################################
###############################################################################

params <- c("trend","totalN","fit", "fit.new","anndet")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 25000
nc <- 4



###############################################################################
####   RUN THE MODEL IN PARALLEL JAGS                   ################################
###############################################################################

model <- jagsUI(bugs.data, inits, params,
		"C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat\\TRENDMODEL_ACT_RandomYear.jags",
		n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=4)



###############################################################################
####   EVALUATE MODEL FIT WITH BAYESIAN P VALUE   #############################
###############################################################################

ylow<-round((min(model$sims.list$fit)-50)/1000,1)*1000
yup<-round((max(model$sims.list$fit)+50)/1000,1)*1000

# pdf(sprintf("MONTSERRAT_%s_model_fit2023.pdf",s), width=10, height=10, title="")
# plot(model$sims.list$fit, model$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(ylow, yup), ylim = c(ylow, yup))
# abline(0, 1, lwd = 2, col = "black")
# dev.off()
pval<-mean(model$sims.list$fit.new > model$sims.list$fit)
slope<-mean(model$mean$fit) / mean(model$mean$fit.new)





###############################################################################
####   PRINT AND SAVE MODEL OUTPUT                #############################
###############################################################################

# Summarize posteriors
#print(model, dig = 3)

write.table(model$summary,sprintf("%s_abund_estimates2023_p%f.csv",s,pval), sep=",")


trendout[trendout$species==s,3]<-round(model$summary[1,5],3)
trendout[trendout$species==s,4]<-round(model$summary[1,3],3)
trendout[trendout$species==s,5]<-round(model$summary[1,7],3)
trendout[trendout$species==s,6]<-pval
trendout[trendout$species==s,7]<-slope

annestimates[annestimates$species==s,3]<-round(model$summary[2:(nyears+1),5],3)
annestimates[annestimates$species==s,4]<-round(model$summary[2:(nyears+1),3],3)
annestimates[annestimates$species==s,5]<-round(model$summary[2:(nyears+1),7],3)
annestimates[annestimates$species==s,6]<-round(model$summary[(nyears+4):(dim(model$summary)[1]-1),5],3)
annestimates[annestimates$species==s,7]<-round(model$summary[(nyears+4):(dim(model$summary)[1]-1),3],3)
annestimates[annestimates$species==s,8]<-round(model$summary[(nyears+4):(dim(model$summary)[1]-1),7],3)

###############################################################################
####   CREATE TREND PLOT AND SAVE AS PDF          #############################
###############################################################################

trendlabel<- paste("Trend: ",trendout[trendout$species==s,3]," (",trendout[trendout$species==s,4]," - ",trendout[trendout$species==s,5],")", sep="")
ggplot(annestimates[annestimates$species==s,], aes(x=Year,y=trend)) +
  geom_line(colour="indianred", linewidth=1.5) +
  geom_ribbon(aes(ymin = lower95CI, ymax = upper95CI), fill="indianred", alpha = 0.2) +

  ## format axis ticks
  scale_x_continuous(name="Year", breaks=seq(2011,2023,2), labels=as.character(seq(2011,2023,2)))+
  #scale_y_continuous(name="Number of Birds at 67 Sampling Points", breaks=seq(0,4000,500), labels=as.character(seq(0,4000,500)))+
  ylab(sprintf("Number of %s at %i sampling points",s,nsites)) +
  
  annotate("text", x = -Inf, y = Inf, label = trendlabel, vjust = 2, hjust = -0.1,size=6, color="black") +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"),
        axis.title.y=element_text(margin=margin(0,20,0,0)), 
        strip.background=element_rect(fill="white", colour="black"))

ggsave(sprintf("MONTSERRAT_%s_abund_plot2023.pdf",s), width=12, height=9)




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









