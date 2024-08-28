#Coordinate cleaning using CoordinateCleaner
#Created: 2024-03-18
#Last Edited: 2024-03-25

#### Clear workspace ####

#rm(list=ls())

#### Load packages ####

library(tidyverse)
library(CoordinateCleaner)
library(elevatr)

#OBS. run scripts in random graphs before this one

#### Elevation ####

#Change the value of "aws_res" to determine resolution of elevation data
#Resolutions range from 1-15, with 15 being the most precise
#For info on resolutions, see https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution

aws_res <- 8

#Any resolution above 8 will yield a warning, due to large file sizes!
#To get around this, remove the number sign at row 51 to activate the line "override_size_check = TRUE"


#Run this to get elevation data according to resolution given above
#Output will be the dataframe "df_1933_ev"
{
#Create dataframe with only coordinates
df_1933_cd<-data.frame(x = df_1933_fl$Longitude,
                       y = df_1933_fl$Latitude)

#Get elevation from coordinates using elevatr
crs_dd<-4326
df_1933_ev_epqs<-data.frame(Elevation = 0)
df_1933_ev_aws<-data.frame(Elevation = 0)

#EPQS elevation
df_1933_ev_epqs<-get_elev_point(df_1933_cd,
                                prj = crs_dd,
                                src = c("epqs"))

#AWS elevation
df_1933_ev_aws<-get_elev_point(df_1933_cd, 
                               prj = crs_dd, 
                               src = "aws",
                               z=aws_res,
                               #override_size_check = TRUE
                               )  

#Compile into one dataframe
df_1933_ev <- data.frame(Elevation.EPQS = df_1933_ev_epqs$elevation,
                         Elevation.AWS = df_1933_ev_aws$elevation)

#Remove excess dataframes
rm(df_1933_cd, df_1933_ev_aws, df_1933_ev_epqs, crs_dd)
rm(df_1933_ev_eqps)
}

#### Coordinate Cleaning ####

#Run this to clean coordinates
#Output will be the dataframe "df_1933_clean"
{
df_1933_fl$Longitude<-as.numeric(df_1933_fl$Longitude)
df_1933_fl$Latitude<-as.numeric(df_1933_fl$Latitude)


df_1933_bfcl<-data.frame(decimalLongitude=df_1933_fl$Longitude,
                         decimalLatitude = df_1933_fl$Latitude,
                         species = "Oxyria digyna")



flags<-clean_coordinates(df_1933_bfcl)


#From this, the only issue seems to be records at sea, with a large majority of flagged specimens
#Let's first add the sea column to the dataframe, as well as elevation
df_1933_flagged <- cbind(df_1933_fl,
                         df_1933_ev,
                         flags = flags$.sea)

#Remove flagged specimens
df_1933_clean <- filter(df_1933_flagged, flags %in% c("TRUE"))

#Remove excess dataframes
rm(df_1933_bfcl, flags, df_1933_flagged)

}

#Visualize on map

ggplot() +
  coord_fixed() +
  borders("world", colour = "gray50", fill = "gray50") +
  geom_point(data = df_1933_clean,
             aes(x = Longitude, y = Latitude),
             colour = "red",
             size = 0.5,
             alpha = 0.5) +
  theme_bw()


#### Remove alpine specimens ####


