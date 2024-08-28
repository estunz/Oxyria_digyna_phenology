#Pull first date of flowering
#Created: 2024-04-02
#Last Edited: 2024-04-02

#### Clear workspace, load libraries ####
{
  #rm(list=ls())
  library(ggplot2)
  library(CoordinateCleaner)
  library(elevatr)
  library(dplyr)
  library(lubridate)
  library(lme4)
  library(tidyverse)
  setwd("C:/Users/zodas/OneDrive/Skrivbord/R-stuff/Bachelor's project/Random graphs")
}


#### Convert MM-DD to days ####

df_1933_dtc <- df_1933_clean
df_1933_dtc$ymd <- paste(df_1933_dtc$Year,
                         df_1933_dtc$Month,
                         df_1933_dtc$Day, 
                         sep = "-")
df_1933_dtc$tot_days <- yday(ymd(df_1933_dtc$ymd))


#### Pull first date of flowering per year ####

df_1933_first <- df_1933_dtc %>%
  group_by(Year) %>%
  arrange(tot_days) %>%
  summarise(first_day = first(tot_days), 
            Latitude = first(Latitude)) %>%
  mutate(Lat_groups = case_when(Latitude <= 60 ~ "<60",
                                Latitude >= 60 & Latitude <= 70 ~ "60-70",
                                Latitude >= 70 ~ ">70"),
         Total_occurrences = 0)

for (year_value in unique(df_1933_clean$Year)) {
  total_rows <- sum(df_1933$Year == year_value)
  df_1933_first$Total_occurrences[df_rowcount$Year == year_value] <- total_rows
}
#### Bin ####


  

#### Graph ####

df_1933_first %>%
  ggplot(aes(x = Year, y = first_day, color = Lat_groups)) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x,
              geom = "smooth",
              alpha = 0,
              lty = "dashed") +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth",
              color = "black",
              alpha = 0.5) +
  theme_bw() +
  labs(y = "First day of flowering",
       x = "Year",
       legend = "Latitude")


regression2 <- lm(formula = first_day ~ Year + Lat_groups, data = df_1933_first)
summary(regression2)

#R^2 = 0.09135 (adjusted). Rather low, don't know what would be reasonable here though
#Residual SE = 15.38
#First quartile = -10.280, median = -2.190, third quartile = 5.121





