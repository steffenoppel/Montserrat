#.....................................................................................
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  CREATION OF ANNUAL REPORT    ##########################################
#############  steffen.oppel@gmail.com      ##########################################
#.....................................................................................

### based on multi-year trend model of Kery et al. 2009
### REQUIRES ANALYIS TO BE COMPLETED IN  "MONTSERRAT_BIRD_MONITORING_Nimble.R"


#.....................................................................................
# 1.  Load required packages       ---------------------
#.....................................................................................

rm(list=ls())
library('rsconnect')

if (!require("pacman")) install.packages("pacman")
# for shiny app
pacman::p_load(tidyverse, #includes tidyr (for using gather(), rearranging data), dplyr, ggplot2
               rgdal,
               raster,
               scales,
               data.table,
               dplyr,
               dtplyr,
               lubridate,
               ggplot2,
               knitr,
               tmap,
               basemaps,
               sf,
               terra,
               shiny,
               rnaturalearth,
               rnaturalearthdata,
               shinydashboard,
               gridExtra,
               readr,
               readxl) #arrange ggplots
# library(tidyverse)
# library(data.table)
# library(dplyr)
# library(dtplyr)
# library(lubridate)
# library(ggplot2)
# library(knitr)
# library(rmarkdown)
# library(tmap)
# library(basemaps)  ### loads basemap in EPSG:3857 and therefore all other sf objects must be st_transform(3857)
# library(sf)
# library(terra)
# library("rnaturalearth")
# library("rnaturalearthdata")
# library(gridExtra)
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

YEAR<- if_else(month(Sys.time())<6,year(Sys.time())-1,year(Sys.time()))


#.....................................................................................
# 2. Set your working directory to find and load model outputs ---------
#.....................................................................................

#setwd("C:\\STEFFEN\\RSPB\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")
#setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")
#setwd("C:\\Users\\sop\\Documents\\Steffen\\RSPB\\Montserrat\\Montserrat")



#.....................................................................................
## 2.1 load the pre-prepared modelling INPUT dataset					--------
#.....................................................................................


file_url <- "https://github.com/steffenoppel/Montserrat/blob/main/data/MONTSERRAT_ANNUAL_DATA_INPUT.RData?raw=true"
load(url(file_url))


#load(sprintf("data/MONTSERRAT_ANNUAL_DATA_INPUT%s.RData",YEAR))
fullnames<-c("Montserrat Oriole", "Forest Thrush", "Bridled Quail-Dove", "Brown Trembler",
             "Antillean Crested Hummingbird","Purple-throated Carib",
             "Pearly-eyed Thrasher","Green-throated Carib","Scaly-breasted Thrasher","Scaly-naped Pigeon",
             "Caribbean Elaenia","Bananaquit")

fullyears<-seq(2011,YEAR,1)



#.....................................................................................
## 2.2 load the pre-prepared modelling OUTPUT estimates					----------
#.....................................................................................


## read in data from local repository
# annestimates<-read.table(sprintf("output/Annual_estimates%s.csv",YEAR), header=T, sep=",")
# trendout<-read.table(sprintf("output/Trend_estimates%s.csv",YEAR), header=T, sep=",")
# mapdata<-fread(sprintf("output/Annual_estimates%s_mapdata.csv",YEAR))

## read in data directly from GitHub (necessary for remote deployment)
urlfile<-sprintf("https://raw.githubusercontent.com/steffenoppel/Montserrat/refs/heads/main/output/Annual_estimates%s.csv",YEAR)
annestimates<-read_csv(url(urlfile),show_col_types = FALSE)
urlfile<-sprintf("https://raw.githubusercontent.com/steffenoppel/Montserrat/refs/heads/main/output/Trend_estimates%s.csv",YEAR)
trendout<-read_csv(url(urlfile),show_col_types = FALSE)
urlfile<-sprintf("https://raw.githubusercontent.com/steffenoppel/Montserrat/refs/heads/main/output/Annual_estimates%s_mapdata.csv",YEAR)
mapdata<-read_csv(url(urlfile),show_col_types = FALSE)

#.....................................................................................
## 2.3 manipulate the loaded data					----------
#.....................................................................................

## add full species names
annestimates$fullspec<-fullnames[match(annestimates$species, SPECIES)]
trendout$fullspec<-fullnames[match(trendout$species, SPECIES)]
mapdata$fullspec<-fullnames[match(mapdata$species, SPECIES)]

## add colour
trendout<-trendout %>%
  mutate(col=ifelse(lcl<0,ifelse(ucl<0,"darkred","black"),ifelse(ucl>0,"forestgreen","black"))) %>%
  mutate(col=ifelse(species=="CAEL","darkred",col))


## add numeric year  
annestimates$Year<-str_replace_all(annestimates$parameter,pattern="[^[:alnum:]]", replacement="")
annestimates$Year<-str_replace_all(annestimates$Year,pattern="totalN", replacement="")
annestimates$Year<-as.numeric(annestimates$Year)+2010

## define dimensions of arrays
nsites<-length(unique(siteCov$Point))
nyears<-length(unique(countdata$year))


## create simple feature

mapdata_sf<-mapdata %>%
  st_as_sf(coords = c("Eastings", "Northings"), crs=2004) %>%
  st_transform(4326)                        ## Montserrat is EPSG 4604 or 2004




#.....................................................................................
# 3. CREATE THE APP ---------
#.....................................................................................


#.....................................................................................
## 3.1 Pre-shiny Preparation -------					
#.....................................................................................


# code to make actionButtons also work by hitting Enter, left arrow and right arrow
# @Stéphane Laurent, https://stackoverflow.com/questions/56600232/detecting-arrow-key-cursor-key-in-shiny

# code to make actionButtons also work by hitting Enter (after clicking on them once)
jscode <- '$(document).keyup(function(event)) {
  if ((event.keyCode == 13)) {
  $("#button").click();}});'

# code formatting mouseover legend in dygraphs
valueFormatter <- "function(x) {
          var options = {weekday: 'short', year: 'numeric', month: '2-digit', day: '2-digit', hour: '2-digit', minute: '2-digit', 
          hour12: false, timeZone: 'UTC'};
          var dayX = new Date(x);
          return dayX.toLocaleString('en-SE', options);
        }"






#.....................................................................................
## 3.2 USER INTERFACE -------					
#.....................................................................................


ui <- fluidPage(
  titlePanel("Montserrat Forest Bird Abundance Estimates"),
  sidebarPanel(
      fluidRow(
        h4("Map of estimated abundance at all forest points"),
        selectInput("species", "Species:", choices = fullnames, selected = "Montserrat Oriole"),
        p("Select the year for which you want to show abundance:"),
        selectInput("year", "Year:", choices = fullyears, selected = max(fullyears)),
        plotOutput("Index")
    )),
    mainPanel(
      tmapOutput("map",height="100vh")
    )
  
)


#.....................................................................................
## 3.2 SERVER SIDE OF THE APP -------					
#.....................................................................................


server <- function(input, output,session) {

  output$Index <-
    renderPlot({
      
      ggplot(data= annestimates %>% 
          filter(Year!=2020, fullspec==input$species) %>%
          mutate(col = as.factor(trendout$col[match(species,trendout$species)])))+
        geom_line(aes(x=Year, y=mean,col=col), linewidth=1)+
        geom_point(aes(x=Year, y=mean,col=col), size=2)+
        geom_errorbar(aes(x=Year, ymin=lcl,ymax=ucl,col=col), width=.1) +
        
        ## remove the legend
        theme(legend.position="none")+
        guides(scale="none",fill=FALSE)+
        theme(legend.title = element_blank())+
        theme(legend.text = element_blank())+
        
        ## format axis ticks
        scale_x_continuous(name="Year", breaks=seq(2011,YEAR,2), labels=as.character(seq(2011,YEAR,2)))+
        ylab(sprintf("Number of birds at %i sampling points",nsites)) +
        scale_color_manual(values = c("black","darkred", "forestgreen"))+
        ## beautification of the axes
        theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text=element_text(size=18, color="black"),
              axis.title=element_text(size=18),
              strip.text.x=element_text(size=18, color="black"),
              axis.title.y=element_text(margin=margin(0,20,0,0)),
              strip.background=element_rect(fill="white", colour="black"))
    })
      
  output$map <- renderTmap({

    pointlabels<-mapdata_sf %>%
      st_drop_geometry() %>%
      filter(Year==max(fullyears)) %>%
      filter(median>0) %>%
      group_by(Point) %>%
      summarise(N_species=length(unique(species))) %>%
      select(Point,N_species)
    
    abundance_sf<-mapdata_sf %>%
      filter(Year==max(fullyears)) %>%
      group_by(Point, geometry) %>%
      summarise(N=sum(median),lcl=sum(lcl),ucl=sum(ucl)) %>%
      mutate(Total= paste0(N, " [",lcl,"-",ucl,"]")) %>%
      left_join(pointlabels, by="Point") %>%
      select(Point,N, Total,N_species)
    
    tm_basemap(c(StreetMap = "OpenStreetMap",
                   TopoMap = "OpenTopoMap")) +

        tm_shape(abundance_sf)  +
        tm_symbols(col="N", palette = c("lightgreen", "darkred"),
                   size=0.8, alpha=0.7,
                   title.col = "Estimated number", zindex=401,popup.vars=c("Total","N_species")) +
      tm_legend(show = TRUE) 
      
        }, env = parent.frame(), quoted = FALSE)
  
  
  
  observe({

    pointlabels<-mapdata_sf %>%
      st_drop_geometry() %>%
      filter(fullspec==input$species) %>%
      group_by(Point) %>%
      summarise(min=min(median),max=max(median),
                yearmin=min(Year[which(median==min(median))]),
                yearmax=max(Year[which(median==max(median))])) %>%
      mutate(History= paste("Min of",min,"in",yearmin,"; Max of",max,"in",yearmax)) %>%
      select(Point,History)
    
    abundance_sf<-mapdata_sf %>%
      filter(fullspec==input$species) %>%
      filter(Year==input$year) %>%
      mutate(Present= paste0(median, " [",lcl,"-",ucl,"]")) %>%
      left_join(pointlabels, by="Point") %>%
      select(Point,mean, Present,History)

    tmapProxy("map", session, {
      tm_remove_layer(401) +
      tm_basemap(c(StreetMap = "OpenStreetMap",
                 TopoMap = "OpenTopoMap")) +

      tm_shape(abundance_sf)  +
        tm_symbols(col="mean", palette = c("lightgreen", "darkred"),
                   size=0.8, alpha=0.7,
                   title.col = "Estimated number", zindex = 401,popup.vars=c("Present","History")) +
        tm_legend(show = TRUE) 

    })
  })
}	


shinyApp(ui, server)