# Project: IWG_NAM_Introduction
# Script 1 - Calculating Growing Degree Days
# Author: Kayla R. Altendorf
# Date: 08/19/2020

# load required packages
library("gridExtra")
library("tidyr")
library("tibble")
library("ggplot2")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("/Growing Degree Day Calculation")


#### Step 1: Prepare Weather Data ####
# Location: STP
# Note: this data was downloaded from NOAA, from Station ID: USC00218450
stp_dat <- read.csv(paste(dir, script, "/data/", "STP Weather Data 2017 and 2018.csv", sep = ""), na.strings = c("."))

stp_dat <- stp_dat %>% 
  separate(DATE, sep = "/", into = c("month", "day", "year")) %>%
  mutate(loc = "STP", 
         year = as.factor(paste("20", year, sep = ""))) %>%
  select(loc, year, month, day, TMAX, TMIN)

# Location: TLI
# Note: this data was obtained from The Land Institute Weather Station
tli16 <- read.csv(paste(dir, script, "/data/", "TLI Weather Data 2016.csv", sep = ""), na.strings = c("."))
tli17 <- read.csv(paste(dir, script, "/data/", "TLI Weather Data 2017.csv", sep = ""), na.strings = c("."))
tli18 <- read.csv(paste(dir, script, "/data/", "TLI Weather Data 2018.csv", sep = ""), na.strings = c("."))
tli_dat <- bind_rows(tli16, tli17, tli18)

# group by date and get min and max per date
tli_dat <- tli_dat %>% 
  separate(col = date_time, into = c("date", "time"), sep = " ") %>%
  group_by(date) %>%
  summarise(TMAX = max(temp_f, na.rm = T), TMIN = min(temp_f, na.rm = T)) %>%
  separate(col = date, into = c("month", "day", "year"), sep = "/") %>%
  mutate(month = as.numeric(as.character(month)),
         day = as.numeric(as.character(day)), 
         loc = "TLI", 
         year = paste("20", year, sep = "")) %>%
  select(loc, year, month, day, TMAX, TMIN) %>%
  filter(year %in% c("2017", "2018"))

# combine data together
all_dat <- rbind(stp_dat, tli_dat) %>% 
  mutate(TMIN = ((as.numeric(TMIN) - 32) / 1.8),
         TMAX = ((as.numeric(TMAX) - 32) / 1.8))

#### Step 2: Calculate GDDs and Format Data for the Figure ####
# set locs and years to iterate through each environment
locs <- c("STP", "TLI")
years <- c("2017", "2018")

# create lists for output
year_output <- list()
loc_output <- list()

# iterate through each environment
for (i in 1:length(locs)) {
  for (j in 1:length(years)) {
    dat1 <- all_dat[all_dat$loc == locs[i],]
    dat2 <- dat1[dat1$year == years[j],]
  
    # if the temperature exceeds 37, change it to 37 as it would not be accumulating GDD beyond that
    for (k in 1:nrow(dat2)) {
      if (! is.na(dat2$TMAX[k]) & dat2$TMAX[k] > 37.0) {dat2$TMAX[k] <- 37.0}
    }
    
    # arrange by month and day, add a day of year and caculate GDD
    dat3 <- dat2 %>% 
      mutate(month = as.numeric(month), 
             day = as.numeric(day)) %>%
      arrange(month, day) %>% 
      mutate(day_of_year = row_number(),
             gdd = (TMAX + TMIN) / 2)
    
    # iterating through each environment separately and filtering that data
    # specifically by when that environment began and stopped accumulating GDDs
    # these values were determined by visual inspection of the data 
    if (locs[i] == "STP" & years[j] == "2017") {
      dat4 <- dat3 %>% filter(day_of_year > 88 & day_of_year < 305)
      for (p in 1:nrow(dat4)) {
        if (dat4$day_of_year[p] < 88) {dat4$gdd[p] <- 0} # before day 88, there were no GDDs accumulated
        else if (dat4$day_of_year[p] > 315) {dat4$gdd[p] <- 0} # same with after day 315
      }
    }
    
    if (locs[i] == "STP" & years[j] == "2018") {
      dat4 <- dat3 %>% filter(day_of_year > 113 & day_of_year < 314)
      for (p in 1:nrow(dat4)) {
        if (dat4$day_of_year[p] < 113) {dat4$gdd[p] <- 0}
        else if (dat4$day_of_year[p] > 316) {dat4$gdd[p] <- 0}
      }
    }
    
    if (locs[i] == "TLI" & years[j] == "2017") {
      dat4 <- dat3 %>% filter(day_of_year > 21 & day_of_year < 320)
      for (p in 1:nrow(dat4)) {
        if (dat4$day_of_year[p] < 21) {dat4$gdd[p] <- 0}
        else if (dat4$day_of_year[p] > 314) {dat4$gdd[p] <- 0}
      }
    }
    
    if (locs[i] == "TLI" & years[j] == "2018") {
      dat4 <- dat3 %>% filter(day_of_year > 32 & day_of_year < 321)
      for (p in 1:nrow(dat4)) {
        if (dat4$day_of_year[p] < 23) {dat4$gdd[p] <- 0}
        else if (dat4$day_of_year[p] > 317) {dat4$gdd[p] <- 0}
      }
    }
    
    year_output[[j]] <- dat4 %>% mutate(GDD_acc = cumsum(replace_na(gdd, 0))) # combining years 
  }
  loc_output[[i]] <- do.call("rbind", year_output) # combining locs
}

output <- do.call("rbind", loc_output) # combining everything



#### Step 3: Create a Dataframe of Events ####
# STP
event <- c("50% Anthesis", "50% Anthesis", "Harvest", "Harvest", "50% Spike Emergence", "50% Spike Emergence")
year <- as.factor(c(2017, 2018, 2017, 2018, 2017, 2018))
date <- c("2017-06-26", "2018-06-21", "2017-08-15", "2018-08-02", "2017-06-08", "2018-06-07")
day_of_year <- c(177, 172, 227, 214, 159, 158)
loc <- rep("STP", 6)

stp_events <- data.frame(event, year, date, day_of_year, loc)

# TLI
event <- c("50% Anthesis", "50% Anthesis", "Harvest", "Harvest", "50% Spike Emergence", "50% Spike Emergence")
year <- as.factor(c(2017, 2018, 2017, 2018, 2017, 2018))
date <- c("2017-06-10", "2018-06-06", "2017-07-25", "2018-07-24", "2017-06-01", "2018-06-07")
day_of_year <- c(161, 157, 206, 205, 152, 158)
loc <- rep("TLI", 6)

tli_events <- data.frame(event, year, date, day_of_year, loc)

# combine events together
events <- rbind(stp_events, tli_events)

# make a column for their GDD_acc
events$GDD_acc <- NA

# determine the GDD accumulation when the event occured, add this to the dataframe 
for (i in 1:nrow(events)) {
  doy <- events$day_of_year[i]
  loc <- events$loc[i]
  year <- events$year[i]
  result <- output %>% filter(day_of_year == doy & loc == loc & year == year)
  gdd_acc <- result$GDD_acc
  events$GDD_acc[i] <- gdd_acc
}


#### Step 3: Make the Figure ####
gdd_plot <- ggplot(data = output, aes(x = day_of_year, y = GDD_acc)) + 
  geom_line(aes(color = year), size = 2) + 
  facet_grid(cols = vars(loc)) +
  geom_segment(data=events, aes(x = 0, y = GDD_acc, xend=day_of_year, yend=GDD_acc, linetype = event, color = year), inherit.aes=FALSE) + # horizontal
  geom_segment(data=events, aes(x = day_of_year, y = 0, xend=day_of_year, yend=GDD_acc, linetype = event, color = year), inherit.aes=FALSE) +  # vertical
  scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
  scale_color_manual(values = c("#000000", "#A9A9A9")) + 
  scale_y_continuous(name = "Accumulated Growing Degree Days (°C)", expand = c(0,0), limits = c(0, 5250)) + 
  scale_x_continuous(name = "Day of Year") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20, margin = margin(t=0, r = 15)), 
        axis.title.y = element_text(size = 20, margin = margin(t=0, r = 15))) 

# save the figure
ggsave("gdd_plot.tiff", 
       plot = last_plot(), 
       path = paste(dir, script, "/output", sep = ""),
       device = "tiff",
       scale = 1, 
       width = 15, 
       height = 10, 
       units = "in")
