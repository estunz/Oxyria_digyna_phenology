#netCDF reading
#Created: 2024-04-09
#Last Edited: 2024-08-30

#### Clear workspace, load packages ####
### Run one section at a time ###
{
  #rm(list=ls())
  
  library(ncdf4)
  library(raster)
  library(ggplot2)
  library(CFtime)
  library(lattice)
  library(RColorBrewer)
  library(dplyr)
  library(tidyverse)
  
  setwd("~/Desktop/GOTHENBURG_materials/Ante_O_digyna_phenology_items/Ante_ms_to_publish")
}
#### Functions #### 
#------------------------------------------------------------------------------#
#           Writes local functions required for downstream analysis            #
#------------------------------------------------------------------------------#

write_temp <- function(array, years, folder_name){
                            
  tmp_df <- data.frame()
  min_month <- 1
  max_month <- years * 12
  while(max_month <= 1440){
    
    rm(tmp_df)
    for(i in min_month:max_month){
      
      #Pull temp. values for the month slice
      tmp_slice <- array[,,i]                                                     #By defining the third dimension of the array, this becomes a 2-dimensional array of cells for the month
      
      #Vector of "tmp" values
      tmp_vec <- as.vector(tmp_slice)                                                 #Reformats the cells of "tmp_slice" into a string vector
      
      #Bind matrix
      if (i == min_month){
          tmp_df <- (cbind(lonlat, tmp_vec, i))                                           #Creates dataframe out of the matrix and vector
    
      }else{
        tmp_df_else <- cbind(lonlat, tmp_vec, i)
        tmp_df <- rbind(tmp_df, tmp_df_else)
        rm(tmp_df_else)
      } #else close
    } # for close
    
    tmp_df <- na.omit(tmp_df)
    tmp_df <- as.data.frame(tmp_df)
    colnames(tmp_df) <- c("Longitude","Latitude", "Temperature (°C)", "Slice") #OBS: Would rather split "Slice" into two columns of "Year" and "Month"
    
    year_month <- t(mapply(conv.ym, tmp_df$Slice))
    colnames(year_month) <- c("year", "month")
    tmp_df$Year <- year_month[, "year"]
    tmp_df$Month <- year_month[, "month"]
    tmp_df <- subset(tmp_df, select = -c(Slice))
    
    if (file.exists(folder_name) == FALSE){
      dir.create(file.path(getwd(), folder_name))
    }
    
    year_span <- paste(min(tmp_df$Year), max(tmp_df$Year), sep = "-")
    write.csv2(tmp_df, file = paste(getwd(), "/cru/", "Temperature_monthly_", year_span, ".csv", sep = ""), row.names = FALSE)
    
    min_month <- min_month + (years * 12)
    max_month <- max_month + (years * 12)
  } #while close
}

conv.ym <- function(months_since_1901) {
  year <- floor((months_since_1901 - 1) / 12) + 1901
  month <- (months_since_1901 - 1) %% 12 + 1
  return(c(year = year, month = month))
}

#### Reformat_temperature ####
#------------------------------------------------------------------------------#
#        Monthly temperature data is available as a 0.5*0.5 deegrees           #
#       global grid in .nc format. This section aquires relevant data          #
#       from the .nc file, and also reformats the data as an array to          #
#                       later be written as .csv                               #
#------------------------------------------------------------------------------#

ncdata <- nc_open(file.path(getwd(), "CRU.nc"))
print(ncdata)

dname <- "tmp"

#Lat and Lon

Latitude <- ncvar_get(ncdata, "lat")
Longitude <- ncvar_get(ncdata, "lon")


nlat <- dim(Latitude)
nlon <- dim(Longitude)

print(c(nlon, nlat))

#Reformat the array "tmp_array" as a matrix
lonlat <- as.matrix(expand.grid(Longitude, Latitude))  

#Time

time <- ncvar_get(ncdata, "time")
time_units <- ncatt_get(ncdata, "time", "units")

ntime <- dim(time)

#Temperature#

tmp_array <- ncvar_get(ncdata, dname)
tmp_fillvalue <- ncatt_get(ncdata, dname, "_FillValue")
tmp_ln <- ncatt_get(ncdata, dname, "long_name")
tmp_units <- ncatt_get(ncdata, dname, "units")

dim(tmp_array)

#Close NetCDF file

nc_close(ncdata)

ls()

#### Prepare_temperature_array ####
#------------------------------------------------------------------------------#
#        Edits the units of time in the array, also removes grid               #
#               cells without values for faster output                         #
#------------------------------------------------------------------------------#

#Reformat time

cf <- CFtime(time_units$value, calendar = "proleptic_gregorian", time) #convert time to CFtime class

timestamps <- as_timestamp(cf) #get character-string times

time_cf <- parse_timestamps(cf, timestamps) #parse the string into date-components ##updated CFtimestamp
#to as_timestamp and CFparse to parse_timestamps

#Replace empty values with NA

tmp_array[tmp_array==tmp_fillvalue$value] <- NA
length(na.omit(as.vector(tmp_array[,,1])))
#length=54025

#### Write_csv ####
#------------------------------------------------------------------------------#
#        Results in a folder "cru" containing .csv files with                  #
#       temperature data spanning 10 years each, for a total of 12 files       #
#------------------------------------------------------------------------------#

write_temp(tmp_array, 10, "cru")  #takes time to run: started at 22.39, done sometime before 23.10

#------------------------------------------------------------------------------#