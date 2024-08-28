#==============================================================================#

#Code for Oxyria analysis
#Created: 2024-03-07
#Last Edited: 2024-07-08

#==============================================================================#

#### Clear environment, load packages, set working directory ####
{
  # rm(list=ls())

  library(ggplot2)  
  library(CoordinateCleaner)
  library(scales)
  library(elevatr)
  library(lme4)
  library(nlme)
  library(sf)
  library(raster)
  library(dplyr)
  library(extrafont)
  library(geosphere)
  library(ncf)
  library(spdep)
  library(spatialreg)
  library(wiqid)
  library(tidyverse)
  loadfonts(device = "win")
}
setwd()
df_temp <- read_csv2(paste(getwd(), ("cru/Temperature_monthly_2001_2010.csv"), sep = "/"))
#==============================================================================#

#### Functions ####

#Function to get elevation from coordinates
get_elevation <- function(data, aws_res, osc = FALSE) {
  if(aws_res > 9 & osc == FALSE){
    message("WARNING: aws_res above 8 downloads very large files. Check 'osc = TRUE' to override")
    return(data)
  }else{
  # Create dataframe with only coordinates
  df_cd <- data.frame(x = data$Longitude, y = data$Latitude)
  
  # Get elevation from coordinates using elevatr
  
  crs_dd <- 4326
  
  #AWS elevation
  if (osc) {
    df_ev_aws <- get_elev_point(df_cd, prj = crs_dd, src = "aws", z = aws_res, override_size_check = TRUE)
  } else {
    df_ev_aws <- get_elev_point(df_cd, prj = crs_dd, src = "aws", z = aws_res)
  }
  
  # EPQS elevation
  df_ev_epqs <- get_elev_point(df_cd, prj = crs_dd, src = "epqs")
  
  # Compile into one dataframe
  df_ev <- data.frame(Elevation.EPQS = df_ev_epqs$elevation,
                      Elevation.AWS = df_ev_aws$elevation,
                      Difference = abs(df_ev_epqs$elevation - df_ev_aws$elevation))
  
  # Check the mean difference of elevation between EPQS and AWS methods
  mean_ev_diff <- mean(df_ev$Difference, na.rm = TRUE)
  
  # Remove excess dataframes
  
  
  return(cbind(data, data.frame(df_ev, mean_ev_diff)))
  
  rm(df_cd, df_ev_aws, df_ev_epqs, crs_dd, df_ev)
  }
}
#READ BELOW:
{
#EPQS uses an API, which does not have global coverage. High precision though
#AWS downloads a global raster
#AWS therefore has global coverage, but higher resolutions take a lot of time
  #Change the value of "aws_res" to set the resolution of
  #Resolutions range from 1-15, with 15 being the most precise
  #For info on resolutions, see https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution

#Any resolution above 8 will not let you proceed, due to large file sizes!
  #Set 'osc = TRUE' when using the function to override
}

#Function to clean and remove flagged specimens
#Note that there are columns referenced that are created through 'get_elevation'
clean_data <- function(data){
  data$Longitude<-as.numeric(data$Longitude)                        #Reform coordinates to be numeric
  data$Latitude<-as.numeric(data$Latitude)
  
  
  bfcl<-data.frame(decimalLongitude=data$Longitude,               #df_bfcl: bfcl stands for before cleaning
                      decimalLatitude = data$Latitude,               #df_bfcl is temporary
                      species = "Oxyria digyna")
  
  
  
  flags <- clean_coordinates(bfcl)                                        
  
  #From this, the only issue seems to be records at sea, with a large majority of flagged specimens
  #Let's first add the sea column to the dataframe, as well as elevation
  df_flagged <- cbind(data, flags = flags$.sea)
  
  #Remove flagged + alpine specimens
    cl <- df_flagged %>%
      filter(flags %in% c("TRUE")) #%>% 
      # subset(select = -c(flags, Environment))
  # Categorize alpine and arctic (arctic tundra > 700m)
  cl$Environment <- "Alpine"  # Default category
  
  for (i in 1:nrow(cl)) {
    if (cl$Latitude[i] >= 60) {
      cl$Environment[i] <- "Arctic tundra"
    }
  }
  for (i in 1:nrow(cl)) {
    if (cl$Elevation.AWS[i] < 700 & cl$Latitude[i] >= 60 & cl$Latitude[i] <= 68) {
      cl$Environment[i] <- "Boreal Arctic"
    }
  }
  
  
  
  #Remove excess dataframes                                                     #We keep df_flagged, since a lot of data is omitted between this one and df_clean
  
  return(cl)
  rm(df_bfcl, flags, cl)
}

#Function to find closest temperature point
find_closest_point <- function(coord1, temp_folder, radius){
  
  filenames <- list.files(temp_folder, pattern="*.csv", full.names=TRUE)
  pst <- data.frame()
  
  
  for(fln in 1:length(filenames)){
    temp <- read_csv2(paste(getwd(), filenames[fln], sep = "/"))
    
    min_year <- min(temp$Year)
    max_year <- max(temp$Year)
    year_span <- filter(coord1, Year >= min_year & Year <= max_year)
    
    if(nrow(year_span) == 0) next  # Skip if there are no rows in the year_span
    
    closest_points <- vector(mode = "list", length =  nrow(year_span))
    min_distances <- numeric(nrow(year_span))
    mean_temperature <- numeric(nrow(year_span))
    june_temperature <- numeric(nrow(year_span))
    
    for(i in 1:nrow(year_span)){
      radius_check = 0
      while(radius_check == 0){
        lon_upper_lim <- min(year_span$Longitude[i] + radius, max(temp$Longitude))
        lon_lower_lim <- max(year_span$Longitude[i] - radius, min(temp$Longitude))
        lat_upper_lim <- min(year_span$Latitude[i] + radius, max(temp$Latitude))
        lat_lower_lim <- max(year_span$Latitude[i] - radius, min(temp$Latitude))
        temp_local <- filter(temp, Year == year_span$Year[i] &
                               Month >= 6 & Month <= 8 &
                               Longitude <= lon_upper_lim &
                               Longitude >= lon_lower_lim &
                               Latitude <= lat_upper_lim &
                               Latitude >= lat_lower_lim )
        if(nrow(temp_local) > 0){
          radius_check <- 1
          distances <- distGeo(year_span[i, c("Longitude", "Latitude")], temp_local[, c("Longitude", "Latitude")])
          min_distances[i] <- min(distances)
          closest_point_index <- which.min(distances)
          
          closest_point_coords <- temp_local[closest_point_index, c("Longitude", "Latitude")]
          
          temp_local_summer <- filter(temp_local, Longitude == closest_point_coords$Longitude &
                                                  Latitude == closest_point_coords$Latitude)
          temp_local_june <- filter(temp_local, 
                                    Longitude == closest_point_coords$Longitude &
                                    Latitude == closest_point_coords$Latitude & 
                                    Month == 6)
          
          # #Filter for months, pull mean temperature between initiated month and june
          # if(year_span$Month[i] <= 5){
          #   # Include May occurrences when the current month is May
          #   temp_local_summer <- filter(temp_local, Month >= 5 & Month <= 6 &
          #                                 Longitude == closest_point_coords$Longitude &
          #                                 Latitude == closest_point_coords$Latitude)
          # } else {
          #   temp_local_summer <- filter(temp_local, Month >= 6 & Month <= year_span$Month[i] &
          #                                 Longitude == closest_point_coords$Longitude &
          #                                 Latitude == closest_point_coords$Latitude)
          # }
          # Calculate the mean temperature for the summer months
          mean_temp_summer <- mean(temp_local_summer$"Temperature (°C)", na.rm = TRUE)
          temp_june <- temp_local_june$"Temperature (°C)"
          
          #closest_points[[i]] <- temp_local[closest_point_index, c("Longitude", "Latitude")]
          june_temperature[i] <- temp_june
          mean_temperature[i] <- mean_temp_summer
        }else{
          radius <- radius + 5
        } 
      }
  }

  year_span$closest_distance <- min_distances
  # year_span$closest_point_lon <- sapply(closest_points, function(x) x[1])
  # year_span$closest_point_lat <- sapply(closest_points, function(x) x[2])
  year_span$June_temperature <- june_temperature
  year_span$Mean_temperature <- mean_temperature
  pst <- rbind(pst, year_span)
    
  }
  
  
  coord1 <- pst
}

#Function to find best neighborhood distance threshold
Neighborhood_generator=function(COOR) {
  models<<-list(
    nb2listw(knn2nb(knearneigh(COOR,1, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,2, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,3, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,4, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,5, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,6, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,7, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,8, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,9, longlat = T)),style="W",zero.policy =T),
    nb2listw(knn2nb(knearneigh(COOR,10, longlat = T)),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,250, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,500, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,750, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1000, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1250, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1500, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,2000, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,2500, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,3000, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,3500, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,4000, longlat = T),style="W",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,250, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,500, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,750, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1000, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1250, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1500, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,2000, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,2500, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,3000, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,3500, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,4000, longlat = T),style="U",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,250, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,500, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,750, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1000, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1250, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,1500, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,2000, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,2500, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,3000, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,3500, longlat = T),style="S",zero.policy =T),
    nb2listw(dnearneigh(COOR, 0,4000, longlat = T),style="S",zero.policy =T)
  )
}
#==============================================================================#

#### Load data ####

#Load ALL Oxyria data from GBIF
df_all <- read.csv2("O. digyna all occurrences.csv")
df_all$decimalLongitude<-as.numeric(df_all$decimalLongitude)
df_all$decimalLatitude<-as.numeric(df_all$decimalLatitude)

sum(is.nan(df_cl_temp))

#Load collected data, omit non-usable data
    #Non-usable data are e.g. non-flowering specimens and   
    #specimens outside the temporal range of flowering
{
df<-read.csv2("Oxyria_data.csv")
  
#Omit months 1-4 + 10-12, as well as uncertain dates
df<-df[df$Month >= 5 & df$Month <= 9 & df$Month != 7.8, ]

#Filter non-flowering specimens
df_fl<-df %>%                                                         #df_fl: fl stands for flowering
  filter(df$Flowers.present == "yes")
}


#Load temperature data

# df_temp_month <- read.csv2("Temperature_monthly.csv")
# 
# #Monthly temperature (mean for northern hemisphere)
# {
# # df_temp_meanmonth <- read_csv("HadCRUT_monthly.csv")
# # 
# # #Separate Year and Month into their own columns
# # df_temp_month <- df_temp_month %>%
# #   mutate(Year = str_extract(Time,"[0-9]{1,4}"),
# #          Month = str_extract(Time,"(?<=-)[0-9]{2}"),
# #          Month = str_remove(Month, "^0+")) %>%
# #   subset(select = -(Time)) %>%
# #   relocate(Year, Month)
# # 
# # colnames(df_temp_month) [c(3:5)] <- paste("month", colnames(df_temp_month)[c(3:5)], sep = "_")
# }

#==============================================================================#


#### Add elevation ####

df_fl <- get_elevation(df_fl, aws_res = 7, osc = TRUE)

#### Clean dataset ####

df_cl_unique <- distinct(df_cl, Occurrence.ID, .keep_all = TRUE)

merged_df <- left_join(df_fl, df_cl_unique[, c("Occurrence.ID", "Elevation.AWS")], by = "Occurrence.ID")

df_cl <- clean_data(merged_df)

length(which(df_cl$Environment == "Alpine"))

#Visualize on map

world_plot <- ggplot() +
  coord_fixed() +
  borders("world", colour = "black", fill = "gray50") +
  scale_x_continuous(breaks = seq(-180, 180, 30)) +
  scale_y_continuous(breaks = seq(-90, 90, 15)) +
  geom_point(data = df_all,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "pink",
             size = 1,
             alpha = 0.25) +
  geom_point(data = df_cl,
             aes(x = Longitude, y = Latitude),
             colour = "red",
             size = 1,
             alpha = 0.25) + 
  theme_bw() +
  labs(title = "Distribution of Oxyria digyna GBIF Occurrences",
       y = "Latitude",
       x = "Longitude") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave(path = "D:/R-stuff/Random graphs/Graphs",
       filename = "Distribution of O. digyna.png",
       plot = world_plot,
       units = "px",
       width = 1080,
       height = 540,
       scale = 3)
#==============================================================================#

#### Pull the first day of flowering from each year ####

#Convert MM-DD to days
            
df_dtc <- df_cl                                                                 #df_dtc: dtc stands for date-count
df_dtc$ymd <- paste(df_dtc$Year,
                    df_dtc$Month,
                    df_dtc$Day, 
                    sep = "-")
df_dtc$tot_days <- yday(ymd(df_dtc$ymd))


#Create datasets containing first day of flowering per year

{df_first <- df_dtc %>% 
  group_by(Year) %>%
  arrange(tot_days) %>%
  summarise(Month = first(Month),
            first_day = first(tot_days), 
            Latitude = first(Latitude),
            Longitude = first(Longitude),
            Elevation = first(Elevation.AWS),
            #Temperature = first(`Temperature (°C)`),
            #closest_distance = first(closest_distance),
            Environment = first(Environment)
            ) %>%
  mutate(
  # Lat_groups = case_when(Latitude < 60 ~ "< 60",
  #                               Latitude >= 60 & Latitude <= 68 ~ "60 - 68",
  #                               Latitude > 68 ~ "> 68"),
         Total_occurrences = 0)

#df_first$Lat_groups = as.factor(df_first$Lat_groups)

df_first <- find_closest_point(df_first, "cru", 1.5)
}

{df_first_alpine <- df_dtc %>% 
  filter(Environment == "Alpine") %>%
    group_by(Year) %>%
    arrange(tot_days) %>%
    summarise(Month = first(Month),
              first_day = first(tot_days), 
              Latitude = first(Latitude),
              Longitude = first(Longitude),
              Elevation = first(Elevation.AWS),
              #Temperature = first(`Temperature (°C)`),
              #closest_distance = first(closest_distance)
    ) %>%
    mutate(
      # Lat_groups = case_when(Latitude < 60 ~ "< 60",
      #                             Latitude >= 60 & Latitude <= 68 ~ "60 - 68",
      #                             Latitude > 68 ~ "> 68"),
           Total_occurrences = 0)
  
  #df_first_alpine$Lat_groups = as.factor(df_first_alpine$Lat_groups)
  
df_first_alpine <- find_closest_point(df_first_alpine, "cru", 1.5)
}

{df_first_arc_tun <- df_dtc %>% 
    filter(Environment == "Arctic tundra") %>%
    group_by(Year) %>%
    arrange(tot_days) %>%
    summarise(Month = first(Month),
              first_day = first(tot_days), 
              Latitude = first(Latitude),
              Longitude = first(Longitude),
              Elevation = first(Elevation.AWS),
              #Temperature = first(`Temperature (°C)`),
              #closest_distance = first(closest_distance)
    ) %>%
    mutate(
      # Lat_groups = case_when(Latitude < 60 ~ "< 60",
      #                             Latitude >= 60 & Latitude <= 68 ~ "60 - 68",
      #                             Latitude > 68 ~ "> 68"),
      Total_occurrences = 0)
  
  #df_first_arc_tun$Lat_groups = as.factor(df_first_arc_tun$Lat_groups)

  df_first_arc_tun <- find_closest_point(df_first_arc_tun, "cru", 1.5)
}

{df_first_bor_arc <- df_dtc %>% 
    filter(Environment == "Boreal Arctic") %>%
    group_by(Year) %>%
    arrange(tot_days) %>%
    summarise(Month = first(Month),
              first_day = first(tot_days), 
              Latitude = first(Latitude),
              Longitude = first(Longitude),
              Elevation = first(Elevation.AWS),
              #Temperature = first(`Temperature (°C)`),
              #closest_distance = first(closest_distance)
    ) %>%
    mutate(
      # Lat_groups = case_when(Latitude < 60 ~ "< 60",
      #                             Latitude >= 60 & Latitude <= 68 ~ "60 - 68",
      #                             Latitude > 68 ~ "> 68"),
      Total_occurrences = 0)
  
  #df_first_bor_arc$Lat_groups = as.factor(df_first_bor_arc$Lat_groups)
  
  df_first_bor_arc <- find_closest_point(df_first_bor_arc, "cru", 1.5)
}


#Add total number of occurrences per year into previously created column "Total_occurrences"

for (year_value in unique(df_cl$Year)) {
  total_rows <- sum(df_cl$Year == year_value)
  df_first$Total_occurrences[df_first$Year == year_value] <- total_rows
}

for (year_value in unique(df_cl$Year)) {
  total_rows <- sum(df_cl$Year == year_value & df_cl$Environment == "Alpine")
  df_first_alpine$Total_occurrences[df_first_alpine$Year == year_value] <- total_rows
}

for (year_value in unique(df_cl$Year)) {
  total_rows <- sum(df_cl$Year == year_value & df_cl$Environment == "Arctic tundra")
  df_first_arc_tun$Total_occurrences[df_first_arc_tun$Year == year_value] <- total_rows
}

for (year_value in unique(df_cl$Year)) {
  total_rows <- sum(df_cl$Year == year_value & df_cl$Environment == "Boreal Arctic")
  df_first_bor_arc$Total_occurrences[df_first_bor_arc$Year == year_value] <- total_rows
}

#==============================================================================#

#### Linear Regression ####


data_subset <- "_alpine" #"", "_arc_tun", "_bor_arc" or "_alpine"
#Filtering

df_cor <- filter(get(paste("df_first", data_subset, sep = "")), Total_occurrences >= 0)
df_cor <- na.omit(df_cor)
#Coefficient test

coefficient <- cor.test(df_cor$first_day, df_cor$Mean_temperature)

print(coefficient)


#Regression test
  regression <- lm(formula = first_day ~ Mean_temperature, data = df_cor)
  summary(regression)
  
  
#==============================================================================#

#### Plots LM ####
{
custom_colors <- c("All" = "black", "Arctic tundra" = "#0AB1D6", "Boreal Arctic" = "#418A00",  "Alpine" = "#C81A00")
legend_order <- c("Arctic tundra", "Boreal Arctic", "Alpine")  
#Flowering time over the years
plot1 <- df_first %>%
  mutate(Environment = factor(Environment, levels = legend_order)) %>%
  ggplot(aes(x = Year, y = first_day, color = Environment)) +
  scale_x_continuous(breaks = c(1901, 1920, 1940, 1960, 1980, 2000, 2020)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              linewidth = 1) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_bw() +
  labs(title = substitute(paste(italic("Oxyria digyna "), "flowering timing 1901 - 2020")),
       y = "First flowering occurrence (total days)",
       x = "Year") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))
print(plot1)
#==============================================================================#

#Flowering by Mean summer temperature
plot2 <- df_first %>%
  mutate(Environment = factor(Environment, levels = legend_order)) %>%
  ggplot(aes(x = Mean_temperature, y = first_day, color = Environment)) +
  scale_x_continuous(breaks = c(min(df_first$Mean_temperature),
                                5, 10, 15, max(df_first$Mean_temperature)),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              linewidth = 1) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_bw() +
  labs(title = substitute(paste(italic("Oxyria digyna "), "flowering timing vs mean summer temperature")),
       y = "First flowering occurrence (total days)",
       x = "Mean temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#

#Flowering by June temperature
plot3 <- df_first %>%
  mutate(Environment = factor(Environment, levels = legend_order)) %>%
  ggplot(aes(x = June_temperature, y = first_day, color = Environment)) +
  scale_x_continuous(breaks = c(min(df_first$June_temperature),
                                0, 5, 10, 15, max(df_first$June_temperature)),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              linewidth = 1) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_bw() +
  labs(title = substitute(paste(italic("Oxyria digyna "), "flowering timing vs mean June temperature")),
       y = "First flowering occurrence (total days)",
       x = "June temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#
# 
# #Flowering time over the years (facet)
# df_first %>%
#   mutate(Environment = factor(Environment, levels = legend_order)) %>%
#   ggplot(aes(x = Year, y = first_day, color = Environment)) +
#   scale_x_continuous(breaks = c(1901, 1920, 1940, 1960, 1980, 2000, 2020)) +
#   geom_point(size = 3,
#              alpha = 0.7) +
#   stat_smooth(method = "lm",
#               formula = y ~ x,
#               geom = "smooth",
#               color = custom_colors[1],
#               weight = 0.5,
#               alpha = 0.5,
#               size = 1) +
#   scale_color_manual(values = custom_colors) +  # Apply custom colors
#   theme_bw() +
#   labs(title = "Flowering by Year",
#        y = "First flowering occurrence (total days)",
#        x = "Year") +
#   theme(legend.position = "right",
#         text = element_text(family = "Times New Roman"),
#         plot.title = element_text(hjust = 0.5)) +
#   facet_grid(cols = vars(Environment))
# 
# #==============================================================================#
# 
# #Flowering by Mean_temperature (facet)
# df_first %>%
#   mutate(Environment = factor(Environment, levels = legend_order)) %>%
#   ggplot(aes(x = Mean_temperature, y = first_day, color = Environment)) +
#   scale_x_continuous(breaks = c(min(df_first$Mean_temperature), 0,
#                                 5, 10, 15, 16.6, max(df_first$Mean_temperature)),
#                      labels = label_number(accuracy = 0.1)) +
#   geom_point(size = 3,
#              alpha = 0.7) +
#   stat_smooth(method = "lm",
#               formula = y ~ x,
#               geom = "smooth",
#               color = custom_colors[1],
#               weight = 0.5,
#               alpha = 0.5,
#               size = 1) +
#   scale_color_manual(values = custom_colors) +  # Apply custom colors
#   theme_bw() +
#   labs(title = "Flowering by Mean temperature",
#        y = "First flowering occurrence (total days)",
#        x = "Mean temperature (°C)") +
#   theme(legend.position = "right",
#         text = element_text(family = "Times New Roman"),
#         plot.title = element_text(hjust = 0.5)) +
#   facet_grid(cols = vars(Environment))

#==============================================================================#
#ARCTIC TUNDRA
#Flowering by year
plot4 <- df_first_arc_tun %>%
  ggplot(aes(x = Year, y = first_day)) +
  scale_x_continuous(breaks = c(1901, 1920, 1940, 1960, 1980, 2000, 2020),
                     limits = c(1901, 2020)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(color = custom_colors[2],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              size = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Arctic Tundra ", paste(italic("Oxyria digyna"))), paste(" flowering timing 1901 - 2020"))),
       y = "First flowering occurrence (total days)",
       x = "Year") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#

#Flowering by Mean_temperature
plot5 <- df_first_arc_tun %>%
  ggplot(aes(x = Mean_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(min(df_first_arc_tun$Mean_temperature), 5, 10, 15, 
                                max(df_first_arc_tun$Mean_temperature)),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +  
  geom_point(color = custom_colors[2],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              size = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Arctic Tundra ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs mean summer temperature"))),
       y = "First flowering occurrence (total days)",
       x = "Mean temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#

#Flowering by June temperature
plot6 <- df_first_arc_tun %>%
  ggplot(aes(x = June_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(min(df_first_arc_tun$June_temperature), 0, 5, 10, 15, 
                                max(df_first_arc_tun$June_temperature)),
                     limits = c(-2, 14),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(color = custom_colors[2],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              size = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Arctic Tundra ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs mean June temperature"))),
       y = "First flowering occurrence (total days)",
       x = "June temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#
#BOREAL ARCTIC
#Flowering by year
plot7 <- df_first_bor_arc %>%
  ggplot(aes(x = Year, y = first_day)) +
  scale_x_continuous(breaks = c(1901, 1920, 1940, 1960, 1980, 2000, 2020),
                     limits = c(1901, 2020)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(color = custom_colors[3],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              size = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Boreal Arctic ", paste(italic("Oxyria digyna"))), paste(" flowering timing 1901 - 2020"))),
       y = "First flowering occurrence (total days)",
       x = "Year") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#

#Flowering by Mean_temperature
plot8 <- df_first_bor_arc %>%
  ggplot(aes(x = Mean_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(min(df_first_bor_arc$Mean_temperature), 0, 5, 10, 
                                max(df_first_bor_arc$Mean_temperature)),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(color = custom_colors[3],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              size = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Boreal Arctic ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs mean summer temperature"))),
       y = "First flowering occurrence (total days)",
       x = "Mean temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#

#Flowering by June temperature
plot9 <- df_first_bor_arc %>%
  ggplot(aes(x = June_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(min(df_first_bor_arc$June_temperature), 0, 5, 10, 15, 
                                max(df_first_bor_arc$June_temperature)),
                     limits = c(-2, 14),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(color = custom_colors[3],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              linewidth = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Boreal Arctic ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs mean June temperature"))),
       y = "First flowering occurrence (total days)",
       x = "June temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#
#ALPINE
#Flowering by year
plot10 <- df_first_alpine %>%
  ggplot(aes(x = Year, y = first_day)) +
  scale_x_continuous(breaks = c(1901, 1920, 1940, 1960, 1980, 2000, 2020),
                     limits = c(1901, 2020)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(color = custom_colors[4],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              size = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Alpine ", paste(italic("Oxyria digyna"))), paste(" flowering timing 1901 - 2020"))),
       y = "First flowering occurrence (total days)",
       x = "Year") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#

#Flowering by Mean_temperature
plot11 <- df_first_alpine %>%
  ggplot(aes(x = Mean_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(min(df_first_alpine$Mean_temperature), 0, 5, 10, 15, 
                                max(df_first_alpine$Mean_temperature)),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(color = custom_colors[4],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              size = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Alpine ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs mean summer temperature"))),
       y = "First flowering occurrence (total days)",
       x = "Mean temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

#==============================================================================#
#Flowering by June temperature
plot12 <- df_first_alpine %>%
  ggplot(aes(x = June_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(min(df_first_alpine$June_temperature), 0, 5, 10, 15, 
                                max(df_first_alpine$June_temperature)),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(140, 260, 20),
                     limits = c(140, 260)) +
  geom_point(color = custom_colors[4],
             size = 8,
             alpha = 0.7) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = custom_colors[1],
              weight = 0.5,
              alpha = 0.5,
              size = 1) +
  theme_bw() +
  labs(title = expression(atop(paste("Alpine ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs mean June temperature"))),
       y = "First flowering occurrence (total days)",
       x = "June temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))


#Save plots
predictor_names <- c("Year", "Mean_temperature", "June_temperature", "Total_occurrences")
  for(i in 1:4){
    data_subset <- data_names[i]
    for(j in 1:3){
      predictor_variable <- predictor_names[j]
      ggsave(path = "D:/R-stuff/Random graphs/Graphs",
             filename = paste("Flowering_by_", predictor_variable, data_subset, ".png", sep = ""),
             plot = get(paste("plot", j + (i-1)*3, sep = "")),
             units = "px",
             width = 1080,
             height = 800,
             scale = 3)
    }
  }
}
#==============================================================================#



#### SAR ####




world_map <- st_read("D:/R-stuff/Random graphs/ne_110m_land.shp")

data_names <- c("", "arc_tun", "bor_arc", "alpine")
model_list <- list()

#Print plots of all SAR-models
for(i in 1:4){
  data_subset <- data_names[i]
  
  if (data_subset == "") {
    used_data <- filter(get(paste("df_first", data_subset, sep = "")), Total_occurrences >= 0)
  } else {
    used_data <- filter(get(paste("df_first", data_subset, sep = "_")), Total_occurrences >= 0)
  }
  data_merged <- data.frame(cbind(used_data$first_day, 
                                  used_data$Latitude, 
                                  used_data$Longitude, 
                                  used_data$Mean_temperature, 
                                  used_data$June_temperature,
                                  used_data$Year,
                                  used_data$Total_occurrences))
  
  final_data <- data_merged[complete.cases(data_merged),]
  colnames(final_data) <- c("first_day", "Latitude", "Longitude", "Mean_temperature", "June_temperature", "Year", "Total_occurrences")
  
  predictor_names <- c("Year", "June_temperature", "Mean_temperature", "Total_occurrences")
  
  for(dt in 1:4){
    predictor_variable <- predictor_names[dt]
    coordinates <- cbind(final_data$Longitude, final_data$Latitude)
    coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS(projection(world_map)))
    coordinates_transformed = spTransform(coordinates_sp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    coordinates_df <- as.data.frame(coordinates_transformed)
    
    
    final_data$Longitude <- coordinates_df$coords.x1
    final_data$Latitude <- coordinates_df$coords.x2
    neighbourhood_models = Neighborhood_generator(cbind(final_data$Longitude, final_data$Latitude))
    
    AIC_LIST=numeric(length(neighbourhood_models))
    for (j in 1:length(neighbourhood_models)) {
      AIC_LIST[j] = AICc(errorsarlm(first_day ~ get(predictor_variable),
                                  data=final_data, listw = neighbourhood_models[[j]], tol.solve = 1e-12,
                                  zero.policy =T))
    }

    index_best_model = which(AIC_LIST==min(AIC_LIST))
    best_neighbour_model = neighbourhood_models[index_best_model]
    AIC_LIST
    model_name <- ifelse(data_subset == "", 
                         paste("sar_model", "all", predictor_variable, sep = "_"), 
                         paste("sar_model", data_subset, predictor_variable, sep = "_"))
    model_list[[model_name]] <- errorsarlm(first_day ~ get(predictor_variable), 
                                           data = final_data,
                                           listw = best_neighbour_model[[1]], tol.solve = 1e-12, 
                                           zero.policy = TRUE)
    model <- model_list[[((i-1)*4+dt)]]
    final_data$fitted_values <- predict(model)
    
    # Create scatterplot with SAR model fit
    p <- ggplot(final_data, aes(x = get(predictor_variable), y = first_day)) +
      scale_x_continuous(breaks = c(1901, 1920, 1940, 1960, 1980, 2000, 2020)) +
      geom_point(color = custom_colors[i], size = 2) +
      geom_line(aes(y = fitted_values), color = "black", size = 1) +
      stat_smooth(method = "lm",
                  formula = y ~ x,
                  geom = "smooth",
                  color = custom_colors[1],
                  weight = 0.5,
                  alpha = 0.5,
                  size = 1) +
      theme_bw() +
      labs(title = paste("Fitted Values by", predictor_variable, gsub("_", " ", data_subset)),
           y = "Fitted values",
           x = predictor_variable) +
      theme(legend.position = "right",
            text = element_text(family = "Times New Roman"),
            plot.title = element_text(hjust = 0.5))
    
    print(p)
  }
}

#Print summaries of all SAR-models (requires running lines 952 - 1028)
for (model_name in names(model_list)) {
  model <- model_list[[model_name]]
  cat("Model:", model_name, "\n")
  print(summary(model, Nagelkerke = TRUE))
}

# xweight = seq(range(final_data$Year)[1],
#               range(final_data$Year)[2], 0.1)
# yweight = predict(lm_model, list(Year = xweight),type="response")
# plot(final_data$Year, final_data$first_day, pch = 16, xlab =
#        "Year", ylab = "First flowering specimen per Year")
# lines(xweight, yweight,col='red')

# coordinates <- cbind(final_data$Longitude, final_data$Latitude)
# coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS(projection(world_map)))
# coordinates_transformed = spTransform(coordinates_sp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# coordinates_df <- as.data.frame(coordinates_transformed)
# 
# 
# 
# 
# 
# 
# 
# 
# autocorrelation_sar = correlog(final_data$Longitude, final_data$Latitude, sar_model$residuals,
#                                increment=1000, latlon=T, resamp=100)
# plot(autocorrelation_sar)
# 
# summary(sar_model, Nagelkerke = TRUE)
# summary(lm_model)


#### Miscellaneous (not used in thesis) ####

#barchart of monthly distribution

{
  df_monthcount <- data.frame(Month = unique(df_fl$Month), monthcount = 0)
  
  for (month_value in unique(df_fl$Month)) {
    total_rows <- sum(df_fl$Month == month_value)
    df_monthcount$monthcount[df_monthcount$Month == month_value] <- total_rows
  }
  }
df_monthcount %>% 
  ggplot(aes(x=Month,
             y=monthcount)) + 
  geom_col() +
  theme_bw() +
  labs(x="Month",
       y="Occurrences")

#barchart of fraction of no per year
    #Make a datafile containing fraction of yes/no per year
{
  df_yesno <- data.frame(Year = unique(df_clean$Year), yescount = 0, nocount = 0, fraction = 0)
  
  for (year_value in unique(df$Year)) {
    #Count occurrences of "yes" and "no" for each unique Year, as well as the fraction of no
    count_yes <- sum(df$Flowers.present[df$Year == year_value] == "yes")
    count_no <- sum(df$Flowers.present[df$Year == year_value] == "no")
    fraction_no <- count_no / (count_yes + count_no)
    
    #Update df_yesno with the count as well as fraction
    df_yesno$yescount[df_yesno$Year == year_value] <- count_yes
    df_yesno$nocount[df_yesno$Year == year_value] <- count_no
    df_yesno$fraction[df_yesno$Year == year_value] <- fraction_no
    
  }
  view(df_yesno)
}

df_yesno %>%
  ggplot(aes(x=Year, y=fraction)) +
  geom_col() +
  labs(x="Year",
       y="Fraction of non-flowering specimens")


#==============================================================================#

