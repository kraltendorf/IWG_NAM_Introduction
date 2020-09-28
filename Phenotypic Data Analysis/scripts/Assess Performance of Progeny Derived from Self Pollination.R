# Project: IWG_NAM_Introduction
# Analysis - Assess Performance of Progeny Derived from Self Pollination 
# Author: Kayla Altendorf 
# Date: 9/17/2020

# load required packages
library("cowplot")
library("dplyr")
library("gtools")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("/Phenotypic Data Analysis")


#### Step 1: Prepare and Format the Data ####
dat <- read.table(paste(dir, script, "/data/NAM_Data_abridged.txt", sep = ""), header = T, sep = "\t")
selfs <- read.csv(paste(dir, script, "/data/selfs.csv", sep = ""), header = T) %>% rename(germplasm_id = longID)

# create two dataframes for merging purposes 
trait_names <- data.frame(trait_id = c("HDEMPER", "ZDK") , trait_id_full = c("emergence_percent", "anthesis_score") )
feekes <- data.frame(anthesis_score = as.numeric(c("49", "51", "53", "55", "57", "59", "61", "65", "69", "71")), coded = 1:10) 

dat1 <- dat %>% filter(trait_id %in% c("HDEMPER", "ZDK")) %>% # filter by the traits we're using
  left_join(trait_names, by = "trait_id") %>% # rename the traits to their logical names
  dplyr::rename(year = phenotype_year, # rename cols to match personal preference
                famID = family_name) %>%
  mutate(loc = substr(experiment_id, 4, 6), # extract location
         parent = substr(germplasm_id, 6, 6)) %>% # extract parent (C for common, D for donor, P for parent)
  select(famID, germplasm_id, parent, loc, year, rep, trait_id_full, phenotype_value, plant_id) %>%
  pivot_wider(names_from = trait_id_full, values_from = phenotype_value) %>% # make data wide
  select(-plant_id) %>% 
  mutate(loc = str_replace(loc, "SAL", "TLI")) %>% # replace SAL (Salina) with TLI
  left_join(feekes, by = "anthesis_score") %>% # code feekes for easier data analysis
  select(-anthesis_score) %>%
  dplyr::rename(anthesis_score = coded) %>%
  mutate(emergence_percent = as.numeric(emergence_percent))



dat <- dat1 %>% left_join(selfs, "germplasm_id")


# make a figure showing the distribution of selfs and outcrosses in the population and whether they were detectable 
dat <- dat %>% filter(loc %in% c("STP", "TLI") & year %in% c("2017", "2018"))
dat$year <- as.factor(dat$year)

dat <- dat %>% mutate(loc_year = paste(loc, year, sep = "_"))

trait <- c("emergence_percent", "anthesis_score")

# environments
loc <- c("STP", "STP", "TLI", "TLI")
year <- c("2017", "2018", "2017", "2018")
envs <- c("stp17", "stp18", "tli17", "tli18")

# conduct t-tests for selfs
selfs_test <-  data.frame(matrix(NA, nrow = 2, ncol = 4))
colnames(selfs_test) <- envs
rownames(selfs_test) <- trait

# filter selfs
selfs <- filter(dat, self == "self")
f1s <- filter(dat, self != "self")

# iterate through and calculate t-test
for (i in 1:length(trait)) {
  
  for (j in 1:length(envs)){
    
    # filter out by year and env
    selfs_loc <- selfs[selfs$loc == loc[j],]
    selfs_loc_year <- selfs_loc[selfs_loc$year == year[j],]
    
    f1s_loc <- f1s[f1s$loc == loc[j],]
    f1s_loc_year <- f1s_loc[f1s_loc$year == year[j],]
    
    # rename vectors
    self_trait <- selfs_loc_year[, trait[i]]
    f1s_trait <- f1s_loc_year[, trait[i]]
    
    # conduct t-test
    t.test_result <- t.test(self_trait, f1s_trait, var.equal = FALSE)
    selfs_test[i, j] <- t.test_result$p.value
  }
}

# set trait and env
selfs_test_ed <- selfs_test

for (i in 1:nrow(selfs_test)) {
  for (j in 1:ncol(selfs_test)) {
    p <- as.numeric(selfs_test[i,j])
    
    if (p > 0.05) {selfs_test_ed[i,j] <- "NS" }
    else {selfs_test_ed[i, j] <- stars.pval(as.numeric(p))
    }
  }
}


selfs_test_final <- selfs_test_ed %>% rownames_to_column("trait") %>%
  pivot_longer(cols = stp17:tli18, names_to = "envs", values_to = "sig")

selfs_test_final <- selfs_test_final %>% 
  mutate(loc = toupper(substr(envs, 1, 3)), 
         year = paste("20", substr(envs, 4, 5), sep = ""))


#### Step 2: Create Plots ####
# set label positions
ymax <- max(dat[,"emergence_percent"], na.rm = T) 
dat_fig <- dat %>% dplyr::select(loc_year, self, emergence_percent)
colnames(dat_fig)[3] <- c("trait_name")

# change trait and fig number
fig1 <- ggplot(dat_fig, aes(x = loc_year, y = trait_name)) + 
  geom_violin()  +
  geom_jitter(data = filter(dat_fig, ! is.na(self)), size = 1.5, aes(color = self), width = 0.1) + 
  annotate("text", label = paste(selfs_test_ed[i,j], sep = ""), x =1, y = ymax, size = 6, colour = "black") + 
  annotate("text", label = paste(selfs_test_ed[i,(j + 1)], sep = ""), x =2, y = ymax, size = 6, colour = "black") + 
  annotate("text", label = paste(selfs_test_ed[i,(j + 2)], sep = ""), x =3, y = ymax, size = 6, colour = "black") + 
  annotate("text", label = paste(selfs_test_ed[i,(j + 3)], sep = ""), x =4, y = ymax, size = 6, colour = "black") + 
  theme_bw() + 
  ylab("Emergence Percent") + 
  xlab("") +
  scale_color_manual(values = c("black", "gray50")) +  
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size=15, face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = 'none')

# set label position and select data
ymax <- max(dat[,"anthesis_score"], na.rm = T) 
dat_fig <- dat %>% dplyr::select(loc_year, self, anthesis_score)
colnames(dat_fig)[3] <- c("trait_name")

# change trait and fig number
fig2 <- ggplot(dat_fig, aes(x = loc_year, y = trait_name)) + 
  geom_violin()  +
  geom_jitter(data = filter(dat_fig, ! is.na(self)), size = 1.5, aes(color = self), width = 0.1) + 
  annotate("text", label = paste(selfs_test_ed[i,j], sep = ""), x =1, y = ymax, size = 6, colour = "black") + 
  annotate("text", label = paste(selfs_test_ed[i,(j + 1)], sep = ""), x =2, y = ymax, size = 6, colour = "black") + 
  annotate("text", label = paste(selfs_test_ed[i,(j + 2)], sep = ""), x =3, y = ymax, size = 6, colour = "black") + 
  annotate("text", label = paste(selfs_test_ed[i,(j + 3)], sep = ""), x =4, y = ymax, size = 6, colour = "black") + 
  theme_bw() + 
  ylab("Anthesis Score") + 
  xlab("") +
  scale_color_manual(values = c("black", "gray50")) +  
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size=15, face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = 'none')

# cowplot them together
plot_grid(fig1, fig2)

# save the figure
ggsave("self_plot.tiff", 
       plot = last_plot(), 
       path = paste(dir, script, "/output", sep = ""),
       device = "tiff",
       scale = 1, 
       width = 15, 
       height = 10, 
       units = "in")

