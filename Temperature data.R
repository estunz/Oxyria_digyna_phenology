#netCDF reading
#Created: 2024-04-09
#Last Edited: 2024-04-09

#### Clear workspace, load packages ####

rm(list=ls())

library(ncdf4)
library(raster)
library(ggplot2)
library(CFtime)
library(lattice)
library(RColorBrewer)
library(dplyr)
library(tidyverse)

#### Functions ####

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

#### Create netCDF file in R? ####

setwd("D:/R-stuff/Temperature/")

ncdata <- nc_open(file.path(getwd(), "output_001.nc"))
print(ncdata)

dname <- "tmp"

#### Lat and Lon ####

Latitude <- ncvar_get(ncdata, "lat")
Longitude <- ncvar_get(ncdata, "lon")


nlat <- dim(Latitude)
nlon <- dim(Longitude)

print(c(nlon, nlat))

lonlat <- as.matrix(expand.grid(Longitude, Latitude))  #Reformats the array "tmp_array" as a matrix

#### Time ####

time <- ncvar_get(ncdata, "time")
time_units <- ncatt_get(ncdata, "time", "units")

ntime <- dim(time)

#### Temperature ####

tmp_array <- ncvar_get(ncdata, dname)
tmp_fillvalue <- ncatt_get(ncdata, dname, "_FillValue")
tmp_ln <- ncatt_get(ncdata, dname, "long_name")
tmp_units <- ncatt_get(ncdata, dname, "units")

dim(tmp_array)

#### Global Attributes ####

# comment <- ncatt_get(ncdata, 0, "comment")
# history <- ncatt_get(ncdata, 0, "history")
# institution <- ncatt_get(ncdata, 0, "institution")
# licence <- ncatt_get(ncdata, 0, "licence")
# reference <- ncatt_get(ncdata, 0, "reference")
# source <- ncatt_get(ncdata, 0, "source")
# title <- ncatt_get(ncdata, 0, "title")
# version <- ncatt_get(ncdata, 0, "version")
# Conventions <- ncatt_get(ncdata, 0, "Conventions")

#### Close NetCDF file ####

nc_close(ncdata)

ls()

#### Change time from days since origin to date ####

cf <- CFtime(time_units$value, calendar = "proleptic_gregorian", time) #convert time to CFtime class

timestamps <- CFtimestamp(cf) #get character-string times

time_cf <- CFparse(cf, timestamps) #parse the string into date-components


#Replace fill values with NA

tmp_array[tmp_array==tmp_fillvalue$value] <- NA
length(na.omit(as.vector(tmp_array[,,1])))

#### FUNCTION ####


write_temp(tmp_array, 10, "cru")  



# tmp_slice <- tmp_array[,,m]
# 
# image(Longitude,Latitude,tmp_slice, col=rev(brewer.pal(10,"RdBu")))





View(tmp_df)



#Convert "Slice" to "year" and "month"


