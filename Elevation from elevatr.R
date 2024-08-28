#Elevation from coordinates using elevatr
#Created: 2024-03-18
#Last Edited: 2024-03-18

#### Load packages ####

library(tidyverse)
library(elevatr)

#### Clear workspace ####

rm(list=ls())

#### Test ####

#Create dataframe with only coordinates
df_1933_cd<-data.frame(x = df_1933_fl$Longitude,
                       y = df_1933_fl$Latitude)
  
#Get elevation from coordinates
crs_dd<-4326
df_1933_ev<-data.frame(Elevation = 0)
df_1933_ev<-get_elev_point(df_1933_cd, prj = crs_dd, src = "epqs")
df_1933_ev_aws<-get_elev_point(df_1933_cd, prj = crs_dd, src = "aws")








