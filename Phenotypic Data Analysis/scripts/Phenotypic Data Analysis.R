# Project: IWG_NAM_Introduction
# Analysis - Prepare and Analyze Phenotypic Data
# Author: Kayla R. Altendorf
# Date: 09/14/2020

# load required packages
library("dplyr") 
library("tidyr")
library("lme4")
library("lmerTest")
library("emmeans")
library("stringr")
library("reshape")
library("plyr")
library("multcompView")
library("scales")
library("cowplot")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("/Phenotypic Data Analysis")


# declare where you want the output to go 
out_path <- paste(dir, script, "/output", sep = "")
in_path <- paste(dir, script, "/data", sep = "")

# read in NAM data from the IWG database
dat <- read.table(paste(dir, "Phenotypic Data Analysis/data/NAM_Data_abridged.txt", sep = ""), header = T, sep = "\t")

#### Step 1: Format the Data ####
# create two dataframes for merging purposes 
trait_names <- data.frame(trait_id = c("HDEMPER", "ZDK") , trait_id_full = c("emergence_percent", "anthesis_score") )
feekes <- data.frame(anthesis_score = as.numeric(c("49", "51", "53", "55", "57", "59", "61", "65", "69", "71")), coded = 1:10) 

dat1 <- dat %>% filter(trait_id %in% c("HDEMPER", "ZDK")) %>% # filter by the traits we're using
  left_join(trait_names, by = "trait_id") %>% # rename the traits to their logical names
  dplyr::rename(year = phenotype_year, # rename cols to match personal preference
                famID = family_name) %>%
  mutate(loc = substr(experiment_id, 4, 6), # extract location
         parent = substr(germplasm_id, 6, 6)) %>% # extract parent (C for common, D for donor, P for parent)
  dplyr::select(famID, germplasm_id, parent, loc, year, rep, trait_id_full, phenotype_value, plant_id) %>%
  pivot_wider(names_from = trait_id_full, values_from = phenotype_value) %>% # make data wide
  dplyr::select(-plant_id) %>% 
  mutate(loc = str_replace(loc, "SAL", "TLI")) %>% # replace SAL (Salina) with TLI
  left_join(feekes, by = "anthesis_score") %>% # code feekes for easier data analysis
  dplyr::select(-anthesis_score) %>%
  dplyr::rename(anthesis_score = coded) %>%
  mutate(emergence_percent = as.numeric(emergence_percent))



#### Step 2: Calculate Parental Means ####

# calculate the family summaries for each parent
fam_sum <- dat1 %>% filter(parent == "P") %>%
    group_by(loc, year, famID) %>%
    summarise_at(c("emergence_percent", "anthesis_score"), mean, na.rm = TRUE)

# separate out common and donor parent means
common_parent_means <- fam_sum %>% filter(famID == "WGN59") 
donor_parent_means <- fam_sum %>% filter(famID != "WGN59") %>% mutate(emergence_percent_diff = NA, 
                                                                     anthesis_score_diff = NA)
# calculate the difference between each within each environment
for (i in 1:nrow(donor_parent_means)) {
  if (donor_parent_means$loc[i] == "STP" & donor_parent_means$year[i] == "2017") {
    donor_parent_means$emergence_percent_diff[i] <- abs(donor_parent_means$emergence_percent[i] - common_parent_means$emergence_percent[1])
    donor_parent_means$anthesis_score_diff[i] <- abs(donor_parent_means$anthesis_score[i] - common_parent_means$anthesis_score[1])
  }
  if (donor_parent_means$loc[i] == "STP" & donor_parent_means$year[i] == "2018") {
    donor_parent_means$emergence_percent_diff[i] <- abs(donor_parent_means$emergence_percent[i] - common_parent_means$emergence_percent[2])
    donor_parent_means$anthesis_score_diff[i] <- abs(donor_parent_means$anthesis_score[i] - common_parent_means$anthesis_score[2])
  }
  if (donor_parent_means$loc[i] == "TLI" & donor_parent_means$year[i] == "2017") {
    donor_parent_means$emergence_percent_diff[i] <- abs(donor_parent_means$emergence_percent[i] - common_parent_means$emergence_percent[3])
    donor_parent_means$anthesis_score_diff[i] <- abs(donor_parent_means$anthesis_score[i] - common_parent_means$anthesis_score[3])
  }
  if (donor_parent_means$loc[i] == "TLI" & donor_parent_means$year[i] == "2018") {
    donor_parent_means$emergence_percent_diff[i] <- abs(donor_parent_means$emergence_percent[i] - common_parent_means$emergence_percent[4])
    donor_parent_means$anthesis_score_diff[i] <- abs(donor_parent_means$anthesis_score[i] - common_parent_means$anthesis_score[4])
  }
}

#### Step 3: Variable Parents Predicting Variable Progeny ####
# read in emmeans from previous analysis
traits <- c("emergence_percent", "anthesis_score")
year <- c("2017", "2018", "2017", "2018")
loc <- c("STP", "STP", "TLI", "TLI")

# add columns to this dataframe for output
donor_parent_means <- donor_parent_means %>% mutate(emergence_percent_sd = 0, 
                                                    anthesis_score_sd = 0) 

for (i in 1:length(traits)) {
  emmeans <- list() # empty this list before a new trait
  files <- list.files(paste(in_path, traits[i], sep = "/"), "emmeans_genet", full.names = T)
  
  for (j in 1:length(files)) {
  emmeans[[j]] <- read.table(files[j], header = T) %>% mutate(year = year[j], loc = loc[j])
  }
  
  emmeans <- do.call("rbind", emmeans)

  for (k in 1:nrow(donor_parent_means)) {
    emmeans1 <- emmeans[emmeans$loc == donor_parent_means$loc[k] ,]
    emmeans2 <- emmeans1[emmeans1$year == donor_parent_means$year[k] ,]
    emmeans3 <- emmeans2[emmeans2$famID == donor_parent_means$famID[k] ,]
    donor_parent_means[k , paste(traits[i], "_sd", sep = "")] <- sd(emmeans3$emmean)
    }
}

# create dataframes for output
cor_df <- data.frame(loc = loc, year = year, p = NA, est = NA)
cor_df <- replicate(2, cor_df, simplify = FALSE)

# test to see if there's a correlation between them
for (i in 1:length(traits)) {
  for (j in 1:4) {
    cor_dat <- donor_parent_means[donor_parent_means$loc == loc[j] ,]
    cor_dat <- as.data.frame(cor_dat[cor_dat$year == year[j] ,])
    
    cor_dat <- cor_dat[,grep(traits[i], colnames(cor_dat))][,-1]
    colnames(cor_dat) <- c("diff", "sd")
    
    cor_out <- cor.test(cor_dat$diff, cor_dat$sd)
  
    cor_df[[i]][j, "p"]  <- cor_out$p.value
    cor_df[[i]][j, "est"] <- cor_out$estimate
    cor_df[[i]]$trait <- traits[i]
  }
}

cor_df <- do.call("rbind", cor_df)
cor_df <- pivot_wider(cor_df, names_from = c(trait), values_from = c(p, est)) %>% dplyr::select(loc, year, est_emergence_percent, p_emergence_percent, 
                                                                                         est_anthesis_score, p_anthesis_score)
write.csv(cor_df, paste(out_path, "variable_parents_progeny.csv", sep = "/"), row.names = F)


#### Step 4: T-test Between Progeny Derived from Common or Donor Parents ####
# merge with the backbone IDs to get their maternal parent origin

# prepare backbone
backbone <- read.csv(paste(in_path, "/backbone.csv", sep = ""), header = T)
backbone1 <- backbone %>% select(famID, parent, plantID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  distinct() %>%
  select(-plantID3, -famID)

# create a vector of families
families <- unique(emmeans2$famID)


for (i in 1:length(traits)) {
  emmeans <- list() # empty this list before a new trait
  files <- list.files(paste(in_path, traits[i], sep = "/"), "emmeans_genet", full.names = T)
  
  for (j in 1:length(files)) {
    emmeans[[j]] <- read.table(files[j], header = T) %>% mutate(year = year[j], loc = loc[j])
  }
  
  emmeans <- do.call("rbind", emmeans)
  emmeans1 <- emmeans %>% mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) # create a column for joining
  emmeans2 <- left_join(emmeans1, backbone1, by = "famID_plantID3")
  
  # create output dataframes
  c_df <- data.frame(c_stp17 = 1:10, c_stp18 = NA, c_tli17 = NA, c_tli18 = NA)
  d_df <- data.frame(d_stp17 = 1:10, d_stp18 = NA, d_tli17 = NA, d_tli18 = NA)
  p_df <- data.frame(p_stp17 = 1:10, p_stp18 = NA, p_tli17 = NA, p_tli18 = NA)
  
  for (k in 1:length(loc)) {
    
    emmeans_loc <- emmeans2[emmeans2$loc == loc[k], ]
    emmeans_loc_year <- emmeans_loc[emmeans_loc$year == year[k], ]
    
    for (l in 1:length(families)) {
      family <- emmeans_loc_year[emmeans_loc_year$famID == families[l], ]
      c <- family[family$parent == "C", ]
      d <- family[family$parent == "D", ]
    
      c_df[l, k] <- mean(c$emmean, na.rm = T)
      d_df[l, k] <- mean(d$emmean, na.rm = T)
    
      t <- t.test(c$emmean, d$emmean)
      stat <- t$statistic
      p <- t$p.value
      
      if (p > 0.05) {
        p_df[l, k] <- "NS"
      }
      else {p_df[l, k] <- p}
    }
  }
  
  # bring all results together
  t_test_out <- cbind(as.data.frame(families), c_df, d_df, p_df) %>% select(families, c_stp17, d_stp17, p_stp17,
                                                   c_stp18, d_stp18, p_stp18,
                                                   c_tli17, d_tli17, p_tli17,
                                                   c_tli18, d_tli18, p_tli18) %>%
    mutate_if(is.numeric, round, 2)
  
  # write out table
  write.csv(t_test_out, paste(out_path,"/", traits[i], "_t_test_between_parents.csv", sep = ""), row.names = F)
}



#### Step 5: Boxplots of Variation within Families ####

# emergence_percent 
# read in data
emmeans <- list() # empty this list before a new trait
files <- list.files(paste(in_path, "emergence_percent", sep = "/"), "emmeans_genet", full.names = T)

for (j in 1:length(files)) {
  emmeans[[j]] <- read.table(files[j], header = T) %>% mutate(year = year[j], loc = loc[j])
}

emmeans_emergence_percent <- do.call("rbind", emmeans)
donor_parent_means <- fam_sum %>% 
  select(-anthesis_score) %>% 
  dplyr::rename(emmean = emergence_percent) %>%
  filter(famID != "WGN59")

common_parent_mean <- fam_sum %>% 
  select(-anthesis_score) %>% 
  dplyr::rename(emmean = emergence_percent) %>%
  filter(famID == "WGN59")

# get the order of means for the first environment
emmeans_emergence_percent %>% filter(loc == "STP", year == "2017") %>% group_by(famID) %>% summarise(avg = mean(emmean)) %>% arrange(avg)
emmeans_emergence_percent$famID <- factor(emmeans_emergence_percent$famID, levels = c("WGN39", "WGN63", "WGN26", "WGN36", 
                                                                                      "WGN38", "WGN15", "WGN07", "WGN45", 
                                                                                      "WGN46", "WGN55"))

colors <- c("#c1d42f", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#770026",  "#68bbaf", "#D55E00", "#0b3d4c", "#CC79A7")



# plot
emp <- ggplot(emmeans_emergence_percent, aes(x = famID, y = emmean)) + 
  geom_boxplot(size = 0.75, aes(middle = mean(emmean)), outlier.size = 0.75) + # change default median line to "mean"
  facet_grid(loc~year, scales = "free_y") +
  geom_point(data = donor_parent_means, aes(color = famID), size = 5) + # donor parent dots
  geom_hline(data = common_parent_mean, aes(yintercept = emmean), linetype = "dotted", size =2, color = "#616365") + # common parent line
  theme_bw() + 
  ylab("Emergence %") +
  xlab("Family ID") +
  scale_color_manual(values = c("#c1d42f", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#770026",  "#68bbaf", "#D55E00", "#0b3d4c", "#CC79A7"))+ 
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size=30, face = "bold"), 
        strip.text.x = element_text(size = 25),
        strip.background.x = element_rect(fill = "white", colour = NA), 
        strip.text.y = element_text(size = 25),
        strip.background.y = element_rect(fill = "white", colour = NA), 
        legend.position = "none")

# export with tiff 1500 x 1000


# anthesis_score
# read in data
emmeans <- list() # empty this list before a new trait
files <- list.files(paste(in_path, "anthesis_score", sep = "/"), "emmeans_genet", full.names = T)

for (j in 1:length(files)) {
  emmeans[[j]] <- read.table(files[j], header = T) %>% mutate(year = year[j], loc = loc[j])
}

emmeans_anthesis_score <- do.call("rbind", emmeans)
donor_parent_means <- fam_sum %>% 
  select(-emergence_percent) %>% 
  dplyr::rename(emmean = anthesis_score) %>%
  filter(famID != "WGN59")

common_parent_mean <- fam_sum %>% 
  select(-emergence_percent) %>% 
  dplyr::rename(emmean = anthesis_score) %>%
  filter(famID == "WGN59")

# get the order of means for the first environment
emmeans_anthesis_score %>% filter(loc == "STP", year == "2017") %>% group_by(famID) %>% summarise(avg = mean(emmean)) %>% arrange(avg)
emmeans_anthesis_score$famID <- factor(emmeans_anthesis_score$famID, levels = c("WGN36", "WGN39", "WGN38", "WGN26", 
                                                                                      "WGN63", "WGN15", "WGN07", "WGN45", 
                                                                                      "WGN46", "WGN55"))
# set function for y axis ticks
integer_breaks <- function(n = 5, ...) {
  breaker <- pretty_breaks(n, ...)
  function(x) {
    breaks <- breaker(x)
    breaks[breaks == floor(breaks)]
  }
}

# plot
ant <- ggplot(emmeans_anthesis_score, aes(x = famID, y = emmean)) + 
  geom_boxplot(size = 0.75, aes(middle = mean(emmean)), outlier.size = 0.75) + # change default median line to "mean"
  facet_grid(loc~year, scales = "free_y") +
  scale_y_continuous(breaks = integer_breaks())+ 
  geom_point(data = donor_parent_means, aes(color = famID), size = 5) + # donor parent dots
  geom_hline(data = common_parent_mean, aes(yintercept = emmean), linetype = "dotted", size =2, color = "#616365") + # common parent line
  theme_bw() + 
  ylab("Anthesis Score") +
  xlab("Family ID") +
  scale_color_manual(values = c("#c1d42f", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#770026",  "#68bbaf", "#D55E00", "#0b3d4c", "#CC79A7"))+ 
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size=30, face = "bold"), 
        strip.text.x = element_text(size = 25),
        strip.background.x = element_rect(fill = "white", colour = NA), 
        strip.text.y = element_text(size = 25),
        strip.background.y = element_rect(fill = "white", colour = NA), 
        legend.position = "none")

# export with tiff 1500 x 1000

boxplots <- plot_grid(emp, ant, labels = c('A', 'B'), label_size = 20)

ggsave("boxplots.tiff", 
       plot = last_plot(), 
       path = out_path,
       device = "tiff",
       scale = 1, 
       width = 20, 
       height = 10, 
       units = "in")

#### Step 6: Calculate Relationship Between Traits ####
emmeans_emergence_percent <- emmeans_emergence_percent %>% mutate(merge_col = paste(famID, plantID3, loc, year, sep = "_")) %>%
  select(merge_col, emmean, famID, loc, year) %>%
  dplyr::rename(emergence_percent = emmean)

emmeans_anthesis_score <- emmeans_anthesis_score %>% mutate(merge_col = paste(famID, plantID3, loc, year, sep = "_")) %>%
  select(merge_col, emmean) %>%
  dplyr::rename(anthesis_score = emmean)

cor_join <- left_join(emmeans_anthesis_score, emmeans_emergence_percent, by = "merge_col") %>% filter(! is.na(loc))
head(cor_join)
# filter by env
stp17 <- filter(cor_join, loc == "STP", year == "2017")
stp18 <- filter(cor_join, loc == "STP", year == "2018")
tli17 <- filter(cor_join, loc == "TLI", year == "2017")
tli18 <- filter(cor_join, loc == "TLI", year == "2018")

# STP 2017
fit.stp17 <- lm(emergence_percent ~ anthesis_score, data = stp17)
fit2.stp17 <- lm(emergence_percent~anthesis_score+I(anthesis_score^2), data = stp17)
summary(fit.stp17)

df <- summary(fit.stp17)
r2.stp17 <- df$r.squared
b1.stp17 <- df$coefficients[1]
b2.stp17 <- "NS"
summary(fit2.stp17) # quadratic term is not significant, indicating linear relationship 

# STP 2018
fit.stp18 <- lm(emergence_percent ~ anthesis_score, data = stp18)
fit2.stp18 <- lm(emergence_percent~anthesis_score+I(anthesis_score^2), data = stp18)
summary(fit.stp18)

df <- summary(fit.stp18)
r2.stp18 <- df$r.squared
b1.stp18 <- df$coefficients[2]
b2.stp18 <- "NS"
summary(fit2.stp18) # quadratic term is not significant, indicating linear relationship 


# TLI 2017
fit.tli17 <- lm(emergence_percent ~ anthesis_score, data = tli17)
fit2.tli17 <- lm(emergence_percent~anthesis_score+I(anthesis_score^2), data = tli17)
summary(fit.tli17)
summary(fit2.tli17) # quadratic term is significant, indicating quadratic relationship 

df <- summary(fit2.tli17)
r2.tli17 <- df$r.squared
b1.tli17 <- df$coefficients[2]
b2.tli17 <- df$coefficients[3]

## TLI 2018
fit.tli18 <- lm(emergence_percent ~ anthesis_score, data = tli18)
fit2.tli18 <- lm(emergence_percent~anthesis_score+I(anthesis_score^2), data = tli18)
summary(fit.tli18)
summary(fit2.tli18) # quadratic term is significant, indicating quadratic relationship 

df <- summary(fit2.tli18)
r2.tli18 <- df$r.squared
b1.tli18 <- df$coefficients[2]
b2.tli18 <- df$coefficients[3]


# create dataframe for r2 
r2 <- data.frame(loc = c("STP", "STP", "TLI", "TLI"), 
                 year = c("2017", "2018", "2017", "2018"), 
                 r2 = c(r2.stp17, r2.stp18, r2.tli17, r2.tli18),
                 b1 = c(b1.stp17, b1.stp18, b1.tli17, b1.tli18),
                 b2 = c(b2.stp17, b2.stp18, b2.tli17, b2.tli18))

# create model output for graphing the relationships
# STP 2017
fit.dat_stp17 <- data.frame(anthesis_score = seq(min(stp17$anthesis_score, na.rm=T), max(stp17$anthesis_score, na.rm=T), length.out = 100),
                            loc = rep("STP"),
                            year = rep(2017))
fit.dat_stp17$emergence_percent<- predict(fit.stp17, newdata = fit.dat_stp17)

# STP 2018
fit.dat_stp18<- data.frame(anthesis_score = seq(min(stp18$anthesis_score, na.rm=T), max(stp18$anthesis_score, na.rm=T), length.out = 100),
                           loc = rep("STP"),
                           year = rep(2018))
fit.dat_stp18$emergence_percent<- predict(fit.stp18, newdata = fit.dat_stp18)

# TLI 2017
fit.dat_tli17<- data.frame(anthesis_score = seq(min(tli17$anthesis_score, na.rm=T), max(tli17$anthesis_score, na.rm=T), length.out = 100),
                           loc = rep("TLI"),
                           year = rep(2017))
fit.dat_tli17$emergence_percent<- predict(fit2.tli17, newdata = fit.dat_tli17)

# TLI 2018
fit.dat_tli18<- data.frame(anthesis_score = seq(min(tli18$anthesis_score, na.rm=T), max(tli18$anthesis_score, na.rm=T), length.out = 100),
                           loc = rep("TLI"),
                           year = rep(2018))
fit.dat_tli18$emergence_percent<- predict(fit2.tli18, newdata = fit.dat_tli18)


# rbind them together
model_fits <- rbind(fit.dat_stp17, fit.dat_stp18, fit.dat_tli17, fit.dat_tli18)

# calculate the mean emergence percent per environment 
mean_emergence_percent <- cor_join %>% group_by(year, loc) %>% summarise(avg_emergence = mean(emergence_percent))


# make the correlation plot 
ggplot(cor_join, aes(anthesis_score, emergence_percent)) + 
  geom_point(size = 2, color = "grey20") + 
  geom_line(data = model_fits, colour = "#0072B2", size = 2) + 
  facet_grid(loc ~ year, scales = "fixed") +
  geom_hline(data = mean_emergence_percent, aes(yintercept = avg_emergence), linetype = "dotted", size = 2, color = "#770026") + # horizontal line will be the mean for each env
  geom_hline(yintercept = 0.5, size = 2, color = "#770026") + # horizontal line will be the mean for each env
  theme_bw() + 
  ylab("Emergence %") +
  xlab("Anthesis") +
  geom_text(data = r2,  size = 8, mapping = aes(x = -Inf, y = 1.65, label = paste("r^2 = ", substr(r2, 1, 4), sep= ""), hjust = -0.1, vjust = -1)) + # add correlation results
  theme(axis.text.x = element_text(size = 25, hjust = 1), 
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size=30, face = "bold"), 
        strip.text.x = element_text(size = 25),
        strip.background.x = element_rect(fill = "white", colour = NA), 
        strip.text.y = element_text(size = 25),
        strip.background.y = element_rect(fill = "white", colour = NA), 
        legend.position = "none")

# save the figure
ggsave("relationship_between_traits.tiff", 
       plot = last_plot(), 
       path = out_path,
       device = "tiff",
       scale = 1, 
       width = 15, 
       height = 10, 
       units = "in")
