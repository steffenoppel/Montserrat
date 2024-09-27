######################################################################################
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  DATA PREPARATION FOR MULTI-YEAR ANALYSIS   ############################
#############  steffen.oppel@gmail.com      ##########################################
######################################################################################


YEAR<-2024									## ENTER THE YEAR OF THE MOST RECENT MONITORING
SPECIES<-c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA')				## ENTER THE SPECIES FOR WHICH YOU WANT ANALYSES TO BE DONE


######################################################################################
#############  Load required packages       ##########################################
######################################################################################

library(RODBC)
library(reshape)
library(tidyverse)
library(lubridate)
filter<-dplyr::filter
select<-dplyr::select
rename<-dplyr::rename
library(rmarkdown)
library(knitr)


######################################################################################
#############  Import data from the queries in the Access Database        ############
######################################################################################
#setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Montserrat\\Raw_Data")
#setwd("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data")
#setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data")
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data\\RODBC_annual_count_import.r")), wait = TRUE, invisible = FALSE)
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data\\RODBC_annual_count_import.r")), wait = TRUE, invisible = FALSE)


#setwd("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring")
#load("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_RODBC_DATA_INPUT.RData")
# setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Raw_Data")
setwd('C:/Users/filib/Documents/Praktika/Sempach/Montserrat') # for Filibert
setwd()
MTOR_count <- odbcConnectAccess2007('Montserrat_Birds_2024.accdb') # change name of the db to the actual one 
birds <- sqlQuery(MTOR_count, "SELECT * FROM Total_count_allpoints_by_species_by_year_vertical")
obsCov <- sqlQuery(MTOR_count, "SELECT * FROM obsCovariates_simple")
siteCov <- sqlQuery(MTOR_count, "SELECT * FROM SiteCov")
species <- sqlQuery(MTOR_count, "SELECT * FROM lkSpecies")
odbcClose(MTOR_count)

nyears<-YEAR-2010




######################################################################################
#############  PROVIDE BASIC SUMMARY FOR ANNUAL NEWS ARTICLE        ##################
######################################################################################



surveys2024<-obsCov %>% filter(year==2024) %>% rename(Year=year)
birds2024<-birds %>% filter(Year==2024) %>% rename(N=SumOfNumber1) %>% mutate(Point=as.integer(as.character(Point))) %>% filter(!Species %in% c('UNK','NA'))
summary2024<-surveys2024 %>% select(Year, Point, Count, Date) %>%
  left_join(birds2024, by=c('Year','Point','Count')) %>%
  group_by(Count, Species) %>%
  summarise(N=sum(N, na.rm=T))
  
totals2024<-summary2024 %>% group_by(Count) %>%
  summarise(N=sum(N), n_spec=length(unique(Species)))

table1<-summary2024 %>% spread(key=Count, value=N, fill = 0) %>%
  filter(!is.na(Species)) %>%
  mutate(Species=species$Species[match(Species,species$SpeciesCode)]) %>%
  arrange(desc(`1`))




######################################################################################
#############  REMOVE POINTS THAT ARE NOT INDEPENDENT         ########################
######################################################################################
## these points are too close to other points in the forest and may lead to double-counts
## this is arguably not critical for trend monitoring, but was for total population size extrapolation
 
#removal<-c(6,8,109,46,103,92,122,125,72,69,70,90,121,115,124,126,1,2,16,34,10,19,45,50,73,75,81,113, 99,76)		## exhaustive list
#removal<-c(6,8,109,46,103,92,122,125,72,69,70,90,121,115,124,126,1,2,16,34,10,19,45,50,73,75,81,113, 176:300)		## exhaustive list for ANDY (includes 76,99, 174, and 175)
removal<-c(99,76)																			## removes only two points that cause error in JAGS

siteCov<-siteCov[!(siteCov$Point %in% removal),]
obsCov<-obsCov[!(obsCov$Point %in% removal),]
birds<-birds[!(birds$Point %in% removal),]





######################################################################################
#############  RELABEL HABITAT VARIABLE NAMES AND FILL IN ZEROs          #############
######################################################################################


names(siteCov)[2:7]<-c('hab','dbh','dist','treeheight','alt','slope')
siteCov[is.na(siteCov)]<-0
siteCov$ridge<-ifelse(siteCov$Location=="Ridge",1,0)				## convert location to a two-level factor
head(siteCov)



######################################################################################
#############  CREATE NEW NUMERIC VARIABLES FOR RAIN AND WIND AND DAY    #############
######################################################################################

obsCov$rain<- 1
obsCov$rain[is.na(obsCov$Rain)] <- 0

obsCov$wind<-as.ordered(match(obsCov$Wind,levels(obsCov$Wind)))
obsCov$wind<-ifelse(obsCov$wind<3,0,1)					### convert wind into two-level factor yes and no
obsCov$Wind<-NULL
names(obsCov)[c(7,8)]<- c("obs", "skill")
head(obsCov)


### CREATE A NUMERIC VARIABLE FOR THE DAY
obsCov$Year<-as.numeric(format(strptime(obsCov$Date, format="%Y"),format="%Y"))
obsCov$Day<-1
for (y in 2011:YEAR){
startdate<-as.Date(sprintf("%i-03-01 00:00",y))  				#set a startdate, here 1 March (must be smaller than first survey day)
obsCov$Day[obsCov$Year==y]<-as.numeric(as.Date(obsCov$Date[obsCov$Year==y])-startdate) 	#calculate the difference between the start date and the survey date
}

obsCov$Date<-NULL

starttime<-min(as.numeric(obsCov$Time-60)) 				#sets the first survey time-1 min as the min time (here 5:49 am)
obsCov$time<-(as.numeric(obsCov$Time)- starttime)/60  		#calculate the difference between 5:49 am and the time of the survey in minutes (default is in seconds)
obsCov$Time<-NULL

head(obsCov)


################ CREATE A BIRD ACTIVITY COVARIATE (ALL SPECIES) FOR EACH POINT COUNT ##########
#### NEEDS TO BE SCALED FOR EACH POINT TO OVERCOME SPATIAL DIFFERENCES #########

head(birds)
pointmax<- birds %>% group_by(Year,Count,Point) %>%
	summarise(activity=sum(SumOfNumber1, na.rm=T)) %>%
	group_by(Point) %>%
	summarise(max_act=max(activity))

activity<- birds %>% group_by(Year,Count,Point) %>%
	summarise(activity=sum(SumOfNumber1, na.rm=T)) %>%
	inner_join(pointmax) %>%
	mutate(ACT=activity/max_act) %>%
	select(Year,Count,Point,ACT)
names(activity)[1]<-'year'

obsCov<-merge(obsCov,activity, by=c('year','Point','Count'), all.x=T)
head(obsCov)


######################################################################################
#############  ADD ZERO COUNTS FOR TARGET SPECIES FROM CUT DOWN POINT 15 #############
######################################################################################
add<-data.frame(Year=2015, Species=rep(SPECIES,each=3), Point='15', Count=rep(1:3,length(SPECIES)), SumOfNumber1=0)
birds<-rbind(birds,add)

head(obsCov)
add<-data.frame(year=2015, Point='15', Count=rep(1:3), Rain=NA, obs=2, skill=1, rain=0, wind=0, Year=2015, Day=37, time=150, ACT=0.2)
obsCov<-rbind(obsCov,add)



######################################################################################
#############  CREATE BLANK ARRAYS FOR INPUT DATA          ###########################
######################################################################################

COUNTDATA<-data.frame(Year=rep(seq(2011,YEAR,1),length(SPECIES)*3*length(unique(siteCov$Point))),
			Species=rep(SPECIES, each=length(seq(2011:YEAR))*3*length(unique(siteCov$Point))),
			Point=rep(rep(unique(siteCov$Point),each=length(seq(2011:YEAR))), length(SPECIES)*3),
			Count=rep(rep(1:3, each=length(seq(2011:YEAR))*length(unique(siteCov$Point))),length(seq(SPECIES))))

SURVEYDATA<-data.frame(Year=rep(seq(2011,YEAR,1),3*length(unique(siteCov$Point))),
			Point=rep(rep(unique(siteCov$Point),each=length(seq(2011:YEAR))), 3),
			Count=rep(1:3, each=length(seq(2011:YEAR))*length(unique(siteCov$Point))))




######################################################################################
#############  FILL IN THE ACTUAL DATA FOR THE FULL ARRAY     ########################
######################################################################################

COUNTDATA<-merge(COUNTDATA, birds, by=c("Species", "Year", "Point", "Count"), all.x=T)
names(COUNTDATA)[5]<-"N"
COUNTDATA$N[is.na(COUNTDATA$N)]<-0										### sets all NA to 0, which includes counts that did not take place
SURVEYDATA<-merge(SURVEYDATA, obsCov, by=c("Year", "Point", "Count"), all.x=T)


### troubleshoot errors of dimension mismatches ##
# s<-"MTOR"
# bird_s<-COUNTDATA %>% filter(Species==s) %>%
#   arrange(Point,Year,Count)
# dim(bird_s)
# dim(SURVEYDATA)
# # SURVEYDATA <- SURVEYDATA %>% mutate(ID=paste(Year,Point,Count, sep="_"))
# # bird_s <- bird_s %>% mutate(ID=paste(Year,Point,Count, sep="_"))
# # 
# # setdiff(SURVEYDATA$ID,bird_s$ID)
# # duplicates<-SURVEYDATA %>% group_by(ID) %>% summarise(N=length(Day)) %>% filter(N>1)
# # SURVEYDATA %>% filter(ID %in% duplicates$ID)
# # obsCov %>% mutate(count=1) %>%
# #   group_by(year,Point,Count) %>%
# #   summarise(n=sum(count))%>%
# #   dplyr::filter(n>1)


######################################################################################
#############  SAVE WORKSPACE TO RUN ANALYSIS IN CLOUD OR ON 64-bit SYSTEM    ########
######################################################################################
#save.image("C:\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_ANNUAL_DATA_INPUT2023.RData")
#save.image("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\MONTSERRAT_ANNUAL_DATA_INPUT2023.RData")
save.image("C:/Users/filib/Documents/Praktika/Sempach/Montserrat/Montserrat_git/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData") # change last bit of the file name to most recent year


######################################################################################
#############  SIMPLE SUMMARY TABLE FOR MANUSCRIPT     ########################
######################################################################################
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")
summary<-COUNTDATA %>% group_by(Species,Year) %>%
  summarise(mean=mean(N, na.rm=T), sd=sd(N, na.rm=T)) %>%
  filter(Year!=2020) %>% ### no counts were done in 2020
  mutate(Species=fullnames[match(Species,SPECIES)])
summary


#write.table(summary,"Montserrat_raw_count_summaries2023.csv", row.names=F, sep=",")
#write.table(summary,"clipboard", row.names=F, sep="\t")




######################################################################################
#############  SIMPLE plot for the key species     ########################
######################################################################################



ggplot()+
geom_line(data=summary, aes(x=Year, y=mean), linewidth=1)+
  facet_wrap(~Species, ncol=2, scales="free_y")+
geom_point(data=summary, aes(x=Year, y=mean), size=1.5,col='black')+
#geom_ribbon(data=annestimates,aes(x=Year, ymin=lower95CI,ymax=upper95CI),alpha=0.2)+
geom_errorbar(data=summary,aes(x=Year, ymin=(mean-0.5*sd),ymax=(mean+0.5*sd)), colour="black", width=.1) +

## remove the legend
theme(legend.position="none")+
guides(fill=FALSE)+
theme(legend.title = element_blank())+
theme(legend.text = element_blank())+

## format axis ticks
scale_x_continuous(name="Year", breaks=seq(2011,2023,1), labels=as.character(seq(2011,2023,1)))+
#scale_y_continuous(name="Number of Birds at 67 Sampling Points", breaks=seq(0,4000,500), labels=as.character(seq(0,4000,500)))+
ylab("Number of Birds per Sampling Point") +

## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"),
	  axis.title.y=element_text(margin=margin(0,20,0,0)), 
        strip.background=element_rect(fill="white", colour="black")) 

ggsave("Montserrat_Forest_Birds2023.pdf", width=13, height=10)




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

 rmarkdown::render('C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Annual_abundance_report_word.Rmd',
                   output_file = "Montserrat_ForestBird_AnnualSummary2023.docx",
                   output_dir = 'C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring')

