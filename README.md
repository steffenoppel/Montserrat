This repository holds data and code for the annual analysis of Forest Bird Monitoring in the Centre Hills, Montserrat, West Indies, conducted by the Department of Environment Montserrat and supported by the RSPB. 
The yearly analysis of the repeated bird counts is implemented for the twelve most abundant species by N-mixture models in NIMBLE (Code by Steffen Oppel) which automatically run every year at the beginning of June in git actions. 

The current state of this worflow is: [![yearly analysis](https://github.com/filibertmoritz/Montserrat/actions/workflows/yearly_analysis.yml/badge.svg)](https://github.com/filibertmoritz/Montserrat/actions/workflows/yearly_analysis.yml)

## ANNUAL STEPS TO TAKE TO CREATE REPORT

# 1. Manual download of data from ArcGIS Online

Date are collected with the smartphone app 'Survey123', and are stored online. From the ArcGIS Online servers, the data must be downloaded by an RSPB staff member. Go to the feature layer [Montserrat Forest Bird Survey](https://rspbeu.maps.arcgis.com/home/item.html?id=c14299a5c1214bc59658e2b1e06d2e2d) and press the button 'Export Data' -> 'Export to CSV file'. Ensure that the pop-up window specifies the column names to be the 'Display name' and press 'Export'. You will be redirected to a different site where the CSV collection of the feature layer is displayed. Press the button 'Download' in the top right corner. This will download a zip file called 'Montserrat_Forest_Bird_Survey.zip'.

# 2. Manual upload of data to GitHub

This zip file called 'Montserrat_Forest_Bird_Survey.zip' downloaded in Step 1 needs to be uploaded to the [data folder in this repository](https://github.com/steffenoppel/Montserrat/tree/main/data). Enter a simple commit message ("uploaded raw data for year X") and commit the zip file to the data folder.

# 3. Run annual data preparation and analysis

Once the zip file called 'Montserrat_Forest_Bird_Survey.zip' has been uploaded to the [data folder in this repository](https://github.com/steffenoppel/Montserrat/tree/main/data), the analysis can be launched via the automated [workflow 'ANNUAL BIRD TREND ANALYSIS'](https://github.com/steffenoppel/Montserrat/actions/workflows/ANNUAL_ANALYSIS_WORKFLOW.yml). Select the workflow and click the 'Run workflow' button to launch the analysis. This workflow will automatically extract the data from the added zip file, prepare the historic and current data for analysis, and run separate N-mixture models in NIMBLE for each species. It will automatically generate tables, figures, and an output report called 'Montserrat_ForestBird_AnnualReport_20XX.html' in the [output folder of this repository](https://github.com/steffenoppel/Montserrat/tree/main/output).
