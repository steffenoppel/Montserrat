#### Script for automaticall visualising and rendering the output from the Nmix models run in NIMBLE #####
#### this script is part of the automated workflow for an annual report for the Centre Hills Forest Bird Monitoring in Montserrat ####
#### Script written by Filibert Heim, filibert.heim@posteo.de, in Nov 2024 with important parts of the script token and adapted from Steffen Oppel 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Preparations  --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages 
library(tidyverse)
library(data.table) 
library(rmarkdown)
library(nimble)


# load data which has been prepared in the script 'Montserrat_forest_bird_monitoring_data_prep_for_Nmix_nimble.R'
load(file = 'data/Montserrat_forest_bird_monitoring_yearly_NIMBLE_model_data.RData')

# save YEAR as txt to make it available as environmental variable across steps
cat(YEAR, "YEAR.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Load in data from single models and create overview table with trends and annual estimates --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get number of sites 
nsites<-length(unique(siteCov$Point)) 

# save names of Birds 
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")

# collect all files that have been created with output from Nmix models and save data in annestimates and trendout
allout<-list.files(path = 'output/', pattern="trend_estimates.csv")
annestimates<-tibble()
trendout<-tibble()
for (f in allout){
  x<-fread(paste0('output/', f))
  trendout<-trendout %>% bind_rows(x %>% filter(parameter %in% c("trend","trend2")))
  annestimates<-annestimates %>% bind_rows(x %>% filter(substr(parameter,1,6)=="totalN"))
}

# export two files with collected data 
write.table(annestimates, "output/Annual_estimates.csv", row.names=F, sep=",")
write.table(trendout,"output/Trend_estimates.csv", row.names=F, sep=",")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Produce figures for population trend and detection prbability --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# add full species names to tables
annestimates$fullspec<-fullnames[match(annestimates$species, SPECIES)]
trendout$fullspec<-fullnames[match(trendout$species, SPECIES)]

# create a plot with abundance trends 
trendout<-trendout %>%
  mutate(col=ifelse(lcl<0,ifelse(ucl<0,"darkred","black"),ifelse(ucl>0,"forestgreen","black"))) #%>%
#mutate(col=ifelse(species=="CAEL","darkred",col))  ## this is a hack we should not need
annestimates %>% 
  mutate(Year=rep(seq(2011,YEAR), length(allout))) %>%  ## need to futureproof this by making max year dynamic
  filter(Year!=2020) %>%
  mutate(col = as.factor(trendout$col[match(species,trendout$species)])) %>% 
  
  
  ggplot()+
  geom_line(aes(x=Year, y=mean,col=col), linewidth=1)+
  facet_wrap(~fullspec, ncol=2, scales="free_y")+
  geom_point(aes(x=Year, y=mean,col=col), size=2)+
  #geom_ribbon(data=annestimates,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  geom_errorbar(aes(x=Year, ymin=lcl,ymax=ucl,col=col), width=.1) +
  
  ## remove the legend
  theme(legend.position="none")+
  guides(scale="none",fill=FALSE)+
  theme(legend.title = element_blank())+
  theme(legend.text = element_blank())+
  
  ## format axis ticks
  scale_x_continuous(name="Year", breaks=seq(2011,YEAR,2), labels=as.character(seq(2011,YEAR,2)))+   ## need to futureproof this by making max year dynamic
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

ggsave(sprintf("output/Montserrat_ForestBird_Trends%s.pdf",YEAR), width=13, height=16)   ## need to futureproof this by making max year dynamic


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Provide basic summary report for annual news article --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# all this output is created for the most recent YEAR
surveys<-obsCov %>% filter(year==YEAR) %>% 
  mutate(Point=as.integer(as.character(Point)))
birds2024<-countdata %>% filter(year==YEAR) %>% select(-Time,-Rain,-Wind,-day,-time,-activity,-Date,-Time,-VisitID, -Observers) %>%
  gather(key=Species, value=N, -year,-Point,-Count) %>%
  mutate(Point=as.integer(as.character(Point)))
summary<-surveys %>% select(year, Point, Count, Date) %>%
  left_join(birds2024, by=c('year','Point','Count')) %>%
  group_by(Count, Species) %>%
  summarise(N=sum(N, na.rm=T))

totals<-summary %>% group_by(Count) %>%
  summarise(N=sum(N), n_spec=length(unique(Species)))

table1<-summary %>% spread(key=Count, value=N, fill = 0) %>%
  filter(!is.na(Species)) %>%
  mutate(Species=tblSpecies$Species[match(Species,tblSpecies$SpeciesCode)]) %>%
  arrange(desc(`1`))

table2<-trendout %>%
  mutate(dir=ifelse(lcl<0,ifelse(ucl<0,"decrease","stable"),ifelse(ucl>0,"increase","stable"))) %>%
  mutate(dir=ifelse(species=="CAEL","(decrease)",dir)) %>%
  mutate(conf=paste(round(mean,3)," (",round(lcl,3)," - ",round(ucl,3),")", sep="")) %>%
  select(fullspec,dir,conf,BayesP,Rhat) %>%
  arrange(desc(BayesP))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Create HTML report for overall summary report --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




### create HTML report for overall summary report
#Sys.setenv(RSTUDIO_PANDOC="C:/Users/Inge Oppel/AppData/Local/Pandoc")
# Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")
#
# rmarkdown::render('Annual_abundance_report_rawdat.Rmd', ########## THIS DOCUMENT IS NOT UPDATED YET
#                  output_file = "Montserrat_ForestBird_AnnualSummary_rawdat.html")

rmarkdown::render('Annual_abundance_report_modelled.Rmd',
                  output_file = sprintf("output/Montserrat_ForestBird_AnnualReport_%s.html",YEAR)) ## need to futureproof this by making max year dynamic so we save annual reports


