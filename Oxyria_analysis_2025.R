#==============================================================================#
#Oxyria digyna flowering phenology
#Created: 2024-03-07
#Last Edited: 2025-08-04
#==============================================================================#
####Run Temperature_data_Jul_2025.R script first#######

{
  rm(list=ls())

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
  #loadfonts(device = "win")
}
setwd(getwd())

#==============================================================================#


#### Functions ####
#------------------------------------------------------------------------------#
#       Writes local functions required for temperature data                   #
#       filtering and downstream analyses                                      #
#------------------------------------------------------------------------------#

#Get_elevation_for_occurrences
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
###########NOTE:#################

#EPQS uses an API, which does not have global coverage. High precision though
#AWS downloads a global raster
#AWS therefore has global coverage, but higher resolutions take a lot of time
  #Change the value of "aws_res" to set the resolution
  #Resolutions range from 1-15, with 15 being the most precise
  #For info on resolutions, see 
#https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution

#Any resolution above 8 will not let you proceed, due to large file sizes!
  #Set 'osc = TRUE' when using the function to override

#Function to clean and remove flagged specimens
#Note that there are columns referenced that are created through 'get_elevation'
clean_data <- function(data){
  data$Longitude<-as.numeric(data$Longitude)
  data$Latitude<-as.numeric(data$Latitude)
  #Reform coordinates to be numeric
  
  bfcl<-data.frame(decimalLongitude=data$Longitude,
                      decimalLatitude = data$Latitude,               
                      species = "Oxyria digyna")
  #df_bfcl: bfcl stands for before cleaning  
  
  
  flags <- clean_coordinates(bfcl)                                        
  
  #From this, the only issue seems to be records at sea, with a large majority of flagged specimens
  #Let's first add the sea column to the dataframe, as well as elevation
  df_flagged <- cbind(data, flags = flags$.sea)
  
  #Remove flagged and group into Alpine, Boreal Arctic, and Arctic tundra categories
    cl <- df_flagged %>%
      filter(flags %in% c("TRUE"))
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
  

#Remove excess dataframes
  
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


#### Load_data ####
#------------------------------------------------------------------------------#
#         Load O. digyna occurrence data from GBIF and flowering               #
#          specimens retained for analysis to compare distributions            #
#------------------------------------------------------------------------------#

#Load ALL Oxyria data from GBIF
setwd("~/Desktop/GOTHENBURG_materials/Ante_O_digyna_phenology_items/Ante_ms_to_publish")
df_all <- read.csv2("O_digyna_all_occurrences.csv")
df_all$decimalLongitude<-as.numeric(df_all$decimalLongitude)
df_all$decimalLatitude<-as.numeric(df_all$decimalLatitude)

#Load collected data, omit non-usable data
    #Non-usable data are e.g. non-flowering specimens and   
    #specimens outside the temporal range of flowering
{
  df <- read.csv2("Oxyria_data.csv")
  df <- subset(df, select = -c(X, X.1))
  
  #Omit months 1-4 + 10-12, as well as uncertain dates
  df <- df[df$Month >= 5 & df$Month <= 9 & df$Month != 7.8, ]
  #Filter non-flowering specimens
  df_fl <- df %>%                                                               
    filter(df$Flowers.present == "yes")
}

#==============================================================================#


#### Add_elevation ####
#------------------------------------------------------------------------------#
#         get_elevation obtains elevation for occurrences                      #
#------------------------------------------------------------------------------#

df_fl <- get_elevation(df_fl, aws_res = 7, osc = TRUE) #629.5Mb download, takes ~50 min
#ERROR: "API returned an empty repsonse (e.g. location in ocean or not in U.S.). NA returned for elevation"
#EPQS uses API so this is expected, AWS elevations are found for all retained occurrences

write.csv(df_fl, file = "O_digyna_df_fl_elev_2025_2.csv")

#==============================================================================#


#### Clean_dataset ####
#------------------------------------------------------------------------------#
#          GBIF data often has issues regarding unreliable data,               #
#          so here we clean data and plot GBIF occurrences and                 #
#          retained flowering specimen locations                               #
#------------------------------------------------------------------------------#

df_fl_unique <- distinct(df_fl, Occurrence.ID, .keep_all = TRUE)

df_cl <- clean_data(df_fl_unique)
#Flagged 344 of 2582 records, EQ = 0.13.#

write.csv(df_cl, file = "O_digyna_df_fl_elev_clean_2025.csv")

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
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

print(world_plot)

ggsave(path = paste(getwd(), "Plots", sep = "/"),
       filename = "Fig1_Distribution_O_digyna_with_grid.png",
       plot = world_plot,
       units = "px",
       width = 1080,
       height = 540,
       scale = 3)

ggsave(path = paste(getwd(), "Plots", sep = "/"),
       filename = "Fig1_Distribution_O_digyna_no_grid.jpeg",
       plot = world_plot 
       + theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")),
       units = "px",
       width = 1080,
       height = 540,
       scale = 3)

#==============================================================================#


#### Flowering_first_day ####
#------------------------------------------------------------------------------#
#           This section pulls the first day of flowering per year             #
#           for each group. It also adds temperature data to all               #
#           occurrences. Make sure you've already run                          # 
#           Temperature_data_Jul_2025.R!!!                                     #
#------------------------------------------------------------------------------#

#Convert MM-DD to days
#df_dtc creates dataframe for date count            
df_dtc <- df_cl
df_dtc$ymd <- paste(df_dtc$Year,
                    df_dtc$Month,
                    df_dtc$Day, 
                    sep = "-")
df_dtc$tot_days <- yday(ymd(df_dtc$ymd))

#Create datasets containing first day of flowering per year

{
  df_first <- df_dtc %>% 
    group_by(Year) %>%
    arrange(tot_days) %>%
    summarise(Month = first(Month),
              first_day = first(tot_days), 
              Latitude = first(Latitude),
              Longitude = first(Longitude),
              Elevation = first(Elevation.AWS),
              Environment = first(Environment))
  
  df_first <- find_closest_point(df_first, "cru", 1.5)
}

{ 
  df_first_alpine <- df_dtc %>% 
    filter(Environment == "Alpine") %>%
    group_by(Year) %>%
    arrange(tot_days) %>%
    summarise(Month = first(Month),
              first_day = first(tot_days), 
              Latitude = first(Latitude),
              Longitude = first(Longitude),
              Elevation = first(Elevation.AWS))
  
  df_first_alpine <- find_closest_point(df_first_alpine, "cru", 1.5)
}

{
  df_first_arc_tun <- df_dtc %>% 
    filter(Environment == "Arctic tundra") %>%
    group_by(Year) %>%
    arrange(tot_days) %>%
    summarise(Month = first(Month),
              first_day = first(tot_days), 
              Latitude = first(Latitude),
              Longitude = first(Longitude),
              Elevation = first(Elevation.AWS))

  df_first_arc_tun <- find_closest_point(df_first_arc_tun, "cru", 1.5)
}

{
  df_first_bor_arc <- df_dtc %>% 
    filter(Environment == "Boreal Arctic") %>%
    group_by(Year) %>%
    arrange(tot_days) %>%
    summarise(Month = first(Month),
              first_day = first(tot_days), 
              Latitude = first(Latitude),
              Longitude = first(Longitude),
              Elevation = first(Elevation.AWS))
  
  df_first_bor_arc <- find_closest_point(df_first_bor_arc, "cru", 1.5)
}

#==============================================================================#


#### Linear Regression ####
#------------------------------------------------------------------------------#
#         Performs linear regressions on the specified data subset             #
#------------------------------------------------------------------------------#

data_subset <- "" # "", "_arc_tun", "_bor_arc" or "_alpine"
#Filtering

df_reg <- filter(get(paste("df_first", data_subset, sep = "")))
df_reg <- na.omit(df_reg)

#Regression test
  regression <- lm(formula = first_day ~ Mean_temperature, data = df_reg)
  summary(regression)
  
  
#==============================================================================#

  
#### Linear Regression Plots ####
#------------------------------------------------------------------------------#
#         Creates linear regression plots for each group                       #
#------------------------------------------------------------------------------#

{
  custom_colors <- c("All" = "black", "Arctic tundra" = "#0AB1D6", "Boreal Arctic" = "#418A00",  "Alpine" = "#C81A00")
  legend_order <- c("Arctic tundra", "Boreal Arctic", "Alpine")  
  data_names <- c("", "arc_tun", "bor_arc", "alpine")
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
    theme_bw() + theme(panel.border = element_blank(), 
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), 
                                      axis.line = element_line(colour = "black")) +
    labs(title = substitute(paste(italic("Oxyria digyna "), "flowering timing 1901 - 2020")),
          y = "First flowering occurrence (day of the year)",
          x = "Year") +
    theme(legend.position = "right",
          text = element_text(family = "Times New Roman", size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  print(plot1)
#==============================================================================#
    
#Flowering by Mean summer temperature
plot2 <- df_first %>%
  mutate(Environment = factor(Environment, levels = legend_order)) %>%
  ggplot(aes(x = Mean_temperature, y = first_day, color = Environment)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = substitute(paste(italic("Oxyria digyna "), "flowering timing vs. mean summer temperature")),
        y = "First flowering occurrence (day of the year)",
        x = "Mean temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
#==============================================================================#
    
#Flowering by June temperature
plot3 <- df_first %>%
  mutate(Environment = factor(Environment, levels = legend_order)) %>%
  ggplot(aes(x = June_temperature, y = first_day, color = Environment)) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = substitute(paste(italic("Oxyria digyna "), "flowering timing vs. mean June temperature")),
        y = "First flowering occurrence (day of the year)",
        x = "June temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Arctic Tundra ", paste(italic("Oxyria digyna"))), paste(" flowering timing 1901 - 2020"))),
        y = "First flowering occurrence (day of the year)",
        x = "Year") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
#==============================================================================#
    
#Flowering by Mean_temperature
plot5 <- df_first_arc_tun %>%
  ggplot(aes(x = Mean_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Arctic Tundra ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs. mean summer temperature"))),
        y = "First flowering occurrence (day of the year)",
        x = "Mean temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
#==============================================================================#
    
#Flowering by June temperature
plot6 <- df_first_arc_tun %>%
  ggplot(aes(x = June_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8, 10, 12, 14),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Arctic Tundra ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs. mean June temperature"))),
        y = "First flowering occurrence (day of the year)",
        x = "June temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Boreal Arctic ", paste(italic("Oxyria digyna"))), paste(" flowering timing 1901 - 2020"))),
        y = "First flowering occurrence (day of the year)",
        x = "Year") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
#==============================================================================#
    
#Flowering by Mean_temperature
plot8 <- df_first_bor_arc %>%
  ggplot(aes(x = Mean_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Boreal Arctic ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs mean summer temperature"))),
        y = "First flowering occurrence (day of the year)",
        x = "Mean temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
#==============================================================================#
    
#Flowering by June temperature
plot9 <- df_first_bor_arc %>%
  ggplot(aes(x = June_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8, 10, 12, 14, 16),
                      limits = c(-2, 16),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Boreal Arctic ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs. mean June temperature"))),
        y = "First flowering occurrence (day of the year)",
        x = "June temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Alpine ", paste(italic("Oxyria digyna"))), paste(" flowering timing 1901 - 2020"))),
        y = "First flowering occurrence (day of the year)",
        x = "Year") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
#==============================================================================#
    
#Flowering by Mean_temperature
plot11 <- df_first_alpine %>%
  ggplot(aes(x = Mean_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Alpine ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs. mean summer temperature"))),
        y = "First flowering occurrence (day of the year)",
        x = "Mean temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
#==============================================================================#
    
#Flowering by June temperature
plot12 <- df_first_alpine %>%
  ggplot(aes(x = June_temperature, y = first_day)) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),
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
  theme_bw() + theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black")) +
  labs(title = expression(atop(paste("Alpine ", paste(italic("Oxyria digyna"))), paste(" flowering timing vs. mean June temperature"))),
        y = "First flowering occurrence (day of the year)",
        x = "June temperature (°C)") +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman", size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
    
#==============================================================================#
    
#Save plots
predictor_names <- c("Year", "Mean_temperature", "June_temperature")
for(i in 1:4){
  data_subset <- data_names[i]
  if (data_subset == "") {
  } else {
    data_subset <- paste("_", data_subset, sep = "")
  }
  for(j in 1:3){
    predictor_variable <- predictor_names[j]
    ggsave(path = paste(getwd(), "Plots", sep = "/"),
            filename = paste("Fig", j + 1, data_subset, " Flowering_by_", predictor_variable, ".png", sep = ""),
            plot = get(paste("plot", j + (i-1)*3, sep = "")),
            units = "px",
            width = 1080,
            height = 800,
            scale = 3)
  }
 }
}

#==============================================================================#

#### Spatial Autoregressive (SAR) models ####
#------------------------------------------------------------------------------#
#         Creates SAR models, plots, and summaries for each group              #
#------------------------------------------------------------------------------#
world_map <- st_read(paste(getwd(), "ne_110m_land.shp", sep="/"))
    
model_list <- list()
    
#Create SAR-models
for(i in 1:4){
  data_subset <- data_names[i]
      
  if (data_subset == "") {
    used_data <- filter(get(paste("df_first", data_subset, sep = "")))
    } else {
        used_data <- filter(get(paste("df_first", data_subset, sep = "_")))
    }
    data_merged <- data.frame(cbind(used_data$first_day, 
                                  used_data$Latitude, 
                                  used_data$Longitude, 
                                  used_data$Mean_temperature, 
                                  used_data$June_temperature,
                                  used_data$Year))
      
  final_data <- data_merged[complete.cases(data_merged),]
  colnames(final_data) <- c("first_day", "Latitude", "Longitude", "Mean_temperature", "June_temperature", "Year")
      
  predictor_names <- c("Year", "June_temperature", "Mean_temperature")
      
  for(dt in 1:3){
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
    model <- model_list[[((i-1)*3+dt)]]
    final_data$fitted_values <- predict(model)
        
    # Create scatterplot with SAR model fit
    p <- ggplot(final_data, aes(x = get(predictor_variable), y = first_day)) +
      scale_x_continuous(breaks = c(1901, 1920, 1940, 1960, 1980, 2000, 2020)) +
      #scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22)) +
      geom_point(color = custom_colors[i], size = 2) +
      geom_line(aes(y = fitted_values), color = "black", size = 1) +
      stat_smooth(method = "lm",
                  formula = y ~ x,
                  geom = "smooth",
                  color = custom_colors[1],
                  weight = 0.5,
                  alpha = 0.5,
                  size = 1) +
          theme_bw() + theme(panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), 
                             axis.line = element_line(colour = "black")) +
          labs(title = paste("SAR model for", predictor_variable, gsub("_", " ", data_subset)),
               y = "First flowering occurrence (day of the year)",
               x = predictor_variable) +
          theme(legend.position = "right",
                text = element_text(family = "Times New Roman"),
                plot.title = element_text(hjust = 0.5))
        
    print(p)
  }
}

#Print summaries of all SAR-models
for (model_name in names(model_list)) {
  model <- model_list[[model_name]]
  cat("Model:", model_name, "\n")
  print(summary(model, Nagelkerke = TRUE))
}