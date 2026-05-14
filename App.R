#.....................................................................................
#############  MONTSERRAT BIRD MONITORING   ##########################################
#############  MAP ANNUAL ABUNDANCE         ##########################################
#############  steffen.oppel@gmail.com      ##########################################
#.....................................................................................


#.....................................................................................
# 1.  Load required packages       ---------------------
#.....................................................................................

rm(list=ls())
library('rsconnect')
library(shiny)
library(tidyverse)
library(data.table)

library(lubridate)
library(ggplot2)
library(leaflet)
library(gridExtra)
library(readr)


select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

YEAR<- if_else(month(Sys.time())<5,year(Sys.time())-1,year(Sys.time()))


#.....................................................................................
# 2. LOAD AND MANIPULATE DATA  ---------
#.....................................................................................



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

## read in data directly from GitHub (necessary for remote deployment)
urlfile<-"https://raw.githubusercontent.com/steffenoppel/Montserrat/refs/heads/main/output/Annual_estimates.csv"
annestimates<-read_csv(url(urlfile),show_col_types = FALSE)
urlfile<-"https://raw.githubusercontent.com/steffenoppel/Montserrat/refs/heads/main/output/Trend_estimates.csv"
trendout<-read_csv(url(urlfile),show_col_types = FALSE)
urlfile<-"https://raw.githubusercontent.com/steffenoppel/Montserrat/refs/heads/main/output/Annual_estimates_mapdata.csv"
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


#.....................................................................................
# 3. CREATE THE SHINY APP ---------
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
    leafletOutput("map", height="100vh")
  )
  
)


#.....................................................................................
## 3.2 SERVER SIDE OF THE APP -------					
#.....................................................................................


server <- function(input, output, session) {

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
  
  # Render the leaflet output
  
  output$map <- renderLeaflet({
    
    # ----- recreate your summaries EXACTLY -----
    
    
    pointlabels <- mapdata %>%
      filter(fullspec == input$species) %>%
      group_by(Point) %>%
      summarise(
        min = min(median, na.rm = TRUE),
        max = max(median, na.rm = TRUE),
        yearmin = Year[which.min(median)],
        yearmax = Year[which.max(median)],
        .groups = "drop"
      ) %>%
      mutate(
        Lowest = paste("Min of", round(min, 2), "in", yearmin),
        Highest = paste("Max of", round(max, 2), "in", yearmax)
      ) %>%
      select(Point, Lowest, Highest)
    
    abundance_sf <- mapdata %>%
      filter(fullspec == input$species,
             Year == input$year) %>%
      mutate(
        Selected_Year = paste0(round(median, 2), " [", round(lcl,2), "-", round(ucl,2), "]")
      ) %>%
      left_join(pointlabels, by = "Point") %>%
      select(Point, mean, Selected_Year, Lowest, Highest, Eastings, Northings)
    
    # ----- color scale (replaces tm_symbols palette) -----
    pal <- colorNumeric(
      palette = c("cornflowerblue", "darkred"),
      domain = abundance_sf$mean,
      na.color = "transparent"
    )
    
    # ----- leaflet map -----
    leaflet(abundance_sf) %>%
      
      addProviderTiles(providers$OpenStreetMap, group = "OpenStreetMap") %>%
      addProviderTiles(providers$OpenTopoMap, group = "OpenTopoMap") %>%

      addCircleMarkers(
        lng = ~Eastings,
        lat = ~Northings,
        radius = 8,                     # fixed size like tm_symbols size=0.8
        
        stroke = TRUE,        #  turn on border
        color = "black",      #  border colour
        weight = 1,           #  border thickness
        
        fillColor = ~pal(mean),             # dynamic colouring 
        fillOpacity = 0.7,
        
        popup = ~paste0(
          "<b>Selected year:</b> ", Selected_Year, "<br>",
          "<b>", Lowest, "</b><br>",
          "<b>", Highest, "</b>"
        )
        
      ) %>%
      
      # ----- legend (replaces tm_legend) -----
    addLegend(
      "bottomright",
      pal = pal,
      values = ~mean,
      title = "Estimated number",
      opacity = 0.7
    ) %>%
    
    
    addLayersControl(
      baseGroups = c("OpenStreetMap", "OpenTopoMap"),
      options = layersControlOptions(collapsed = FALSE)
    )
    
  })


}

shinyApp(ui, server)


#.....................................................................................
# 5.  DEPLOY APP TO SERVER       ---------------------
#.....................................................................................


# rsconnect::deployApp(appDir="C:/STEFFEN/Vogelwarte/Shinyapps/Montserrat",
# 				appFiles="App.R",
# 				appName="Montserrat",
# 				appTitle="Montserrat")
