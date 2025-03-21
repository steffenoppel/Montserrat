#.....................................................................................
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  CREATION OF ANNUAL REPORT    ##########################################
#############  steffen.oppel@gmail.com      ##########################################
#.....................................................................................

### based on multi-year trend model of Kery et al. 2009
### REQUIRES ANALYIS TO BE COMPLETED IN  "MONTSERRAT_BIRD_MONITORING_Nimble.R"



## THIS APP DOES NOT WORK DUE TO A PROBLEM IN tmap v4: https://github.com/r-tmap/tmap/issues/1024


#.....................................................................................
# 1.  Load required packages       ---------------------
#.....................................................................................

rm(list=ls())
library('rsconnect')

if (!require("pacman")) install.packages("pacman")
# for shiny app
pacman::p_load(tidyverse, #includes tidyr (for using gather(), rearranging data), dplyr, ggplot2
               data.table,
               dplyr,
               dtplyr,
               lubridate,
               ggplot2,
               tmap,
               basemaps,
               sf,
               terra,
               shiny,
               shinydashboard,
               gridExtra,
               readr) #arrange ggplots

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


file_url <- sprintf("https://github.com/steffenoppel/Montserrat/blob/main/data/MONTSERRAT_ANNUAL_DATA_INPUT%s.RData?raw=true",YEAR)
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
  #mutate(col=ifelse(species=="CAEL","darkred",col)) %>%
  group_by(species, fullspec) %>%
  summarise(col=first(col))

## add numeric year  
annestimates$Year<-str_replace_all(annestimates$parameter,pattern="[^[:alnum:]]", replacement="")
annestimates$Year<-str_replace_all(annestimates$Year,pattern="totalN", replacement="")
annestimates$Year<-as.numeric(annestimates$Year)+2010

## define dimensions of arrays
nsites<-length(unique(siteCov$Point))

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


  server <- function(input, output, session) {
    # Create mapdata_sf as a reactive object
    mapdata_sf <- reactive({
      req(mapdata)  # Ensure mapdata is available
      mapdata %>%
        st_as_sf(coords = c("Eastings", "Northings"), crs = 2004) %>%
        st_transform(4326)  # Transform to EPSG 4326
    })
    
    output$Index <-
      renderPlot({
        
        ggplot(data= annestimates %>% 
                 filter(Year!=2020, fullspec==input$species))+
          geom_line(aes(x=Year, y=mean),col=trendout$col[match(input$species,trendout$fullspec)], linewidth=1)+
          geom_point(aes(x=Year, y=mean),col=trendout$col[match(input$species,trendout$fullspec)], size=2)+
          geom_errorbar(aes(x=Year, ymin=lcl,ymax=ucl),col=trendout$col[match(input$species,trendout$fullspec)], width=.1) +
          
          ## remove the legend
          theme(legend.position="none")+
          guides(scale="none",fill="none")+
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
    
    # Render the tmap output
    output$map <- renderTmap({
      # Access the reactive mapdata_sf using mapdata_sf()
      abundance_sf <- mapdata_sf() %>%
        filter(Year == max(fullyears)) %>%
        group_by(Point, geometry) %>%
        summarise(N = sum(median), lcl = sum(lcl), ucl = sum(ucl)) %>%
        mutate(Total = paste0(N, " [", lcl, "-", ucl, "]")) %>%
        select(Point, N, Total)
      
      tm_basemap(c(StreetMap = "OpenStreetMap", TopoMap = "OpenTopoMap")) +
        tm_shape(abundance_sf) +
        tm_symbols(
          col= "darkgray",
          fill = "N",
          fill.scale=tm_scale(values=c("lightgreen", "darkred")),
          size = 0.8,
          fill_alpha = 0.7,
          fill.legend=tm_legend("Estimated number"),
          popup.vars = c("Total")
        )+
        tm_layout(legend.position=c("right", "top"))
    })
    
    
    observe({
      
      pointlabels<-mapdata_sf() %>%
        st_drop_geometry() %>%
        filter(fullspec==input$species) %>%
        group_by(Point) %>%
        summarise(min=min(median),max=max(median),
                  yearmin=min(Year[which(median==min(median))]),
                  yearmax=max(Year[which(median==max(median))])) %>%
        mutate(Lowest= paste("Min of",min,"in",yearmin)) %>%
        mutate(Highest= paste("Max of",max,"in",yearmax)) %>%
        select(Point,Lowest,Highest)
      
      abundance_sf<-mapdata_sf() %>%
        filter(fullspec==input$species) %>%
        filter(Year==input$year) %>%
        mutate(Selected_Year= paste0(median, " [",lcl,"-",ucl,"]")) %>%
        left_join(pointlabels, by="Point") %>%
        select(Point,mean, Selected_Year,Lowest,Highest)
      
      tmapProxy("map", session, {
        tm_basemap(c(StreetMap = "OpenStreetMap", TopoMap = "OpenTopoMap")) +
          tm_shape(abundance_sf) +
          tm_symbols(
            col= "darkgray",
            fill = "mean",
            fill.scale=tm_scale(values=c("lightgreen", "darkred")),
            size = 0.8,
            fill_alpha = 0.7,
            fill.legend=tm_legend("Estimated number"),
            popup.vars = c("Selected_Year","Lowest","Highest")
          ) +
          tm_layout(legend.position=c("right", "top"))
      })
    })
    
    
  }

shinyApp(ui, server)

