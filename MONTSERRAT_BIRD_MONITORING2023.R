######################################################################################
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  ANALYSIS OF ANNUAL SURVEYS   ##########################################
#############  steffen.oppel@gmail.com      ##########################################
######################################################################################

### created 13 April 2014
### based on multi-year trend model of Kery et al. 2009
### code adapted from AQWA monitoring of Oppel et al. 2014
### REQUIRES DATA PREPARATION FROM "MONTSERRAT_BIRD_MONITORING_data_prep.R"
### last update 2 May 2014: include 'rain' (rainfall in previous 24 hrs at nearest rain station) as obsCovariate instead of 'Day'
### updated to 2016 data on 24 May 2016 (without rainfall data but with simple 0/1 rain variable)

### modified 13 April 2017 to fix figures for paper

### modified 24 Dec 2017 to incorporate suggestion by David Borchers to include activity as observation-level covariate

### modified 10 sEPT 2021 to GENERALISE output for multiple years



######################################################################################
#############  Load required packages       ##########################################
######################################################################################

library(tidyverse)
library(reshape)
library(R2jags)
library(jagsUI)
#install.packages("rjags")
library(rjags)
library(ggplot2)
library(knitr)
library(rmarkdown)


######################################################################################
#############  Set your working directory (path where the database is)       #########
######################################################################################

#setwd("C:\\STEFFEN\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring")
#setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring")
setwd("C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat")





######################################################################################
#############  load the pre-prepared dataset					     #########
######################################################################################

load("MONTSERRAT_ANNUAL_DATA_INPUT2023.RData")
#load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_ANNUAL_DATA_INPUT2023.RData")
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")

###############################################################################
############## CREATE SITE COVARIATE DATA INPUT MATRIX   ######################
###############################################################################
nsites<-length(unique(siteCov$Point))

#### CAST THE SITE COV DATA FRAME INTO MATRIX WITH 1 COLUMN PER YEAR
siteCovList<-array(NA, dim=c(nsites,nyears, 3))
selectedCovs<-c(5,6,15)

for (col in 1:length(selectedCovs)){
  dis<-siteCov[,c(1,selectedCovs[col])]
  
  ### STANDARDIZE COVARIATES FOR WINBUGS
  
  meant<-mean(dis[,2], na.rm = TRUE)
  sdt<-sd(dis[,2], na.rm = TRUE)
  dis[,2]<-(dis[,2]-meant)/sdt
  
  dis<-dis[order(dis$Point,decreasing=F),]
  siteCovList[,,col]<-as.matrix(rep(dis[,2],nyears), dimnames=NULL)
}


names(siteCov)[c(5,6,15)]
treeheight<-siteCovList[,,1]
elev<-siteCovList[,,2]
canopy<-siteCovList[,,3]






###############################################################################
############## CREATE OBSERVATION COVARIATE DATA INPUT MATRIX   ###############
###############################################################################
SURVEYDATA<-SURVEYDATA[order(SURVEYDATA$Point,SURVEYDATA$Year, SURVEYDATA$Count, decreasing=F),] 

## SORT THE TABLE SO IT HAS THE SAME ORDER AS THE BIRD DATA

obsCov<-obsCov[order(obsCov$Point,obsCov$Year, obsCov$Count, decreasing=F),] 
head(obsCov)


### STANDARDIZE COVARIATES FOR WINBUGS

meant<-mean(SURVEYDATA$time, na.rm = TRUE)
sdt<-sd(SURVEYDATA$time, na.rm = TRUE)
SURVEYDATA$time<-(SURVEYDATA$time-meant)/sdt

meant<-mean(SURVEYDATA$Day, na.rm = TRUE)
sdt<-sd(SURVEYDATA$Day, na.rm = TRUE)
SURVEYDATA$Day<-(SURVEYDATA$Day-meant)/sdt


### only needs standardisation if measured in mm, not as 0/1 variable
#meant<-mean(SURVEYDATA$rain, na.rm = TRUE)
#sdt<-sd(SURVEYDATA$rain, na.rm = TRUE)
#SURVEYDATA$rain<-(SURVEYDATA$rain-meant)/sdt



### create array for each covariate

ridge<-array(NA, dim=c(nsites,3,nyears))
wind<-array(NA, dim=c(nsites,3,nyears))
time<-array(NA, dim=c(nsites,3,nyears))
ACT<-array(NA, dim=c(nsites,3,nyears))				## REPLACED ON 2 MAY WITH RAINFALL AMOUNT


### fill in array for each covariate
for (y in 2011:YEAR){
  obsC<-subset(SURVEYDATA, Year==y)
  y<-match(y,c(2011:YEAR))						## translates the year (2011, 2012, etc.) into consecutive number (1,2,...) for array dimensions
  x<-cast(obsC, Point ~ Count, value='time')
  x2<-as.matrix(x[,2:4])
  time[,,y]<-x2
  
  x<-cast(obsC, Point ~ Count, value='wind')
  wind[,,y]<-as.matrix(x[,2:4])
  
  x<-cast(obsC, Point ~ Count, value='ACT')			## REPLACED ON 2 MAY WITH RAINFALL AMOUNT - changed to ACTIVITY ON 24 Dec 2017
  ACT[,,y]<-as.matrix(x[,2:4])
  
  ridge[,,y]<-matrix(rep(siteCov$ridge,3), ncol=3)		### site-level obs covariates are constant across years and counts
}





###############################################################################
####   REPLACE ALL NA IN COVARIATES otherwise "undefined node" error    #######
###############################################################################

siteCovList[is.na(siteCovList)]<-0

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


#### TROUBLESHOOT MISMATCH ERROR ###
# errorfind<- SURVEYDATA %>% left_join(bird_s, by=c('Year','Point','Count')) %>%
#   dplyr::filter(is.na(Species))


###############################################################################
############## CREATE BIRD DATA INPUT MATRIX   ################################
###############################################################################

### create array to be filled with data
BIRD.y<-array(NA, dim=c(nsites,3,nyears))
Nst<-array(NA, dim=c(nsites,nyears))

#### CAST THE MOLTEN DATA FRAME INTO MATRIX WITH 1 COLUMN PER COUNT and fill in array
for (y in 1:nyears){
b<-subset(bird_s, Year==c(2011:YEAR)[y])
b$YearPoint<-paste(b$Point,b$Year,sep = "_")		# creates a matching expression for each transect and count
dis<-cast(b, YearPoint~Count, value="N")							# pivot table to create data frame with one line per transect per year, and each column reflecting the observations per distance band on each count survey
dis<-dis[order(dis$YearPoint,decreasing=F),]
BIRD.y[,,y]<-as.matrix(dis[,2:4], dimnames=NULL)
dis[is.na(dis)]<-1
Nst[,y]<-apply(dis[,2:4], MARGIN=1, FUN=max, na.rm=T)+1
}



######################################################################################################
########## CREATE INPUT DATA FOR JAGS
#######################################################################################################

# check data dimensions
#dim(BIRD.y)
#dim(treeheight)
#dim(elev)
#dim(canopy)
#dim(day)
#dim(wind)
#dim(ridge)
#dim(time)



### Bundle data into a single list passed on to JAGS

R = nrow(BIRD.y)
T = ncol(BIRD.y)
nyears = dim(BIRD.y)[3]
bugs.data<-list(M = BIRD.y,
                nsite=nsites,
                nrep=T,
                primocc=seq(1:nyears),
                nyear=nyears,
                elev=elev,
                treeheight=treeheight,
                canopy=canopy,
                wind=wind,
                ridge=ridge,
                time=time,
                ACT=ACT)


###############################################################################
####   SET INITIAL VALUES FOR THE MODEL RUN    ################################
###############################################################################

#Nst <- apply(BIRD.y, c(1, 3), max) + 1
#Nst[is.na(Nst)] <- 1

inits <- function(){list(N = Nst,
loglam = runif(1, -3, 3),
sigma.site = runif(1, 0, 1),
beta.canopy=runif(1,-5,5),
beta.treeheight=runif(1,-5,5),
beta.elev=runif(1,-5,5),
bwind=runif(1,-5,5),
bridge=runif(1,-5,5),
btime=runif(1,-5,5),
bact=runif(1,-5,5),
p0 = runif(nyears, 0, 1))}

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
		"C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat\\TRENDMODEL_ACT_RandomYear.jags",
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



write.table(annestimates,"Annual_estimates2023.csv", row.names=F, sep=",")
write.table(trendout,"Trend_estimates2023.csv", row.names=F, sep=",")




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

ggsave("Montserrat_ForestBird_Trends_2023.pdf", width=13, height=16)






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
                  output_file = "Montserrat_ForestBird_AnnualSummary2023.html",
                  output_dir = 'C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring')

rmarkdown::render('C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat\\Annual_abundance_report_modelled.Rmd',
                  output_file = "Montserrat_ForestBird_AnnualSummary2023.html",
                  output_dir = 'C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat')









