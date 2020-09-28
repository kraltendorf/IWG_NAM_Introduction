# Project: IWG_NAM_Introduction
# Analysis - Marker Summary
# Author: Kayla Altendorf 
# Date: 6/16/2020

# load required packages
library("tidyverse")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("Marker Summary")

# set vectors for iteration 
famID <- c("WGN07", "WGN15", "WGN26", "WGN36", "WGN38", "WGN39", "WGN45", "WGN46", "WGN55", "WGN63")


#### Step 1: Count Markers that were Considered for Map Creation ####

# create an empty dataframe for results
marker_summary <- data.frame(matrix(NA, nrow = 21, ncol = 11))
colnames(marker_summary)[1] <- "chr"
marker_summary$chr <- 1:21
colnames(marker_summary)[2:11] <- famID

# replicate this dataframe eight times
marker_summary_list <- replicate(8, marker_summary, simplify = FALSE)

# reading in loc files to determine how many markers were considered for map creation
for (i in 1:21) {
  for (j in 1:length(famID)) {
    file <- read.table(paste(dir, "/Additional Filtering and JoinMap File Prep/output/", famID[j], "_chr", i, ".txt", sep = ""), skip = 5, fill = T)
    markers <- file %>% filter(grepl("Chr", V1)) %>% dplyr::select(V1, V2)
    tally <- markers %>% group_by(V2) %>% tally()
    tally <- as.data.frame(tally) %>% column_to_rownames("V2")
    
    # add to marker summary list
    marker_summary_list[[1]][i, (j+1)] <- length(markers$V1)
    marker_summary_list[[2]][i, (j+1)] <- tally["<hkxhk>",]
    marker_summary_list[[3]][i, (j+1)] <- tally["<lmxll>",]
    marker_summary_list[[4]][i, (j+1)] <- tally["<nnxnp>",]
    
    # read in map and determine which ones were on the map
    map <- read.table(paste(dir, "/JoinMap/output/Map Files for MapQTL/Formatted for R/no_group.map", sep = "")) 
    markers_on_map <- markers %>% filter(V1 %in% map$V1)
    
    tally_map <- markers_on_map %>% group_by(V2) %>% tally()
    tally_map <- as.data.frame(tally_map) %>% column_to_rownames("V2")
  
    marker_summary_list[[5]][i, (j+1)] <- length(markers_on_map$V1)
    marker_summary_list[[6]][i, (j+1)] <- tally_map["<hkxhk>",]
    marker_summary_list[[7]][i, (j+1)] <- tally_map["<lmxll>",]
    marker_summary_list[[8]][i, (j+1)] <- tally_map["<nnxnp>",]
  }
}

# create a column for marker type and map type (creation = map making; mapping = for linkage mapping)
marker_summary_list[[1]] <- marker_summary_list[[1]] %>% mutate(marker_type = "nloc") %>% mutate(type = "creation")
marker_summary_list[[2]] <- marker_summary_list[[2]] %>% mutate(marker_type = "hkhk") %>% mutate(type = "creation")
marker_summary_list[[3]] <- marker_summary_list[[3]] %>% mutate(marker_type = "lmll") %>% mutate(type = "creation")
marker_summary_list[[4]] <- marker_summary_list[[4]] %>% mutate(marker_type = "nnnp") %>% mutate(type = "creation")

marker_summary_list[[5]] <- marker_summary_list[[5]] %>% mutate(marker_type = "nloc") %>% mutate(type = "mapping")
marker_summary_list[[6]] <- marker_summary_list[[6]] %>% mutate(marker_type = "hkhk") %>% mutate(type = "mapping")
marker_summary_list[[7]] <- marker_summary_list[[7]] %>% mutate(marker_type = "lmll") %>% mutate(type = "mapping")
marker_summary_list[[8]] <- marker_summary_list[[8]] %>% mutate(marker_type = "nnnp") %>% mutate(type = "mapping")

# rbind results
marker_summary_all <- do.call("rbind", marker_summary_list)

marker_summary_long <- pivot_longer(marker_summary_all, WGN07:WGN63, names_to = "family") %>% 
  arrange(marker_type) %>% 
  mutate(facet_var = paste(chr, toupper(substr(type, 1, 1)), sep = "_")) %>%
  mutate(family_type = paste(substr(family, 4, 5), toupper(substr(type, 1,1)), sep = "_"))

marker_summary_long$chr <- as.factor(marker_summary_long$chr)
marker_summary_long$value <- as.numeric(as.character(marker_summary_long$value))

# set the order for facets
marker_summary_long$chr_order = factor(marker_summary_long$chr, levels=c("1", "4", "7", "10", "13", "16", "19", "2", "5", "8", "11", "14", "17", "20", "3", "6", "9", "12", "15", "18", "21"))
head(marker_summary_long)

marker_summary_long$marker_type_ed <- NA

for (i in 1:nrow(marker_summary_long)) {
  if (marker_summary_long$marker_type[i] == "hkhk") {marker_summary_long$marker_type_ed[i] <- "hkxhk"}
  else if (marker_summary_long$marker_type[i] == "lmll") {marker_summary_long$marker_type_ed[i] <- "lmxll"}
  else if (marker_summary_long$marker_type[i] == "nnnp") {marker_summary_long$marker_type_ed[i] <- "nnxnp"}
}

#### Step 2: Make a Figure Displaying Marker Types #### 

ggplot(filter(marker_summary_long, marker_type != "nloc" & type == "mapping" & marker_type != "hkhk"), aes(y=value, x=family)) + 
  geom_bar(aes(fill = marker_type_ed), stat="identity",  width = 0.9) +
  facet_wrap(~chr_order, nrow = 3) +
  theme_bw() +
  ylab("Markers (ct)") + 
  xlab("Family") + 
  #scale_y_continuous(limits = c(0, 300)) +
  scale_fill_manual(values = c("gray40", "gray1"), name = "Marker Type") + #"gray70" for hkhk
  theme(axis.text.x = element_text(size = 11, angle = 90, hjust = 1), 
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size=20, face = "bold"), 
        strip.text.x = element_text(size = 15),
        strip.background.x = element_rect(fill = "white", colour = NA), 
        strip.text.y = element_text(size = 15),
        strip.background.y = element_rect(fill = "white", colour = NA), 
        legend.text=element_text(size=18),
        legend.title=element_text(size=15),
        legend.position="bottom")

# save the figure
ggsave("marker_types.tiff", 
       plot = last_plot(), 
       path = paste(dir, script, "/output", sep = ""),
       device = "tiff",
       scale = 1, 
       width = 15, 
       height = 10, 
       units = "in")

# average marker cts before mapping
marker_summary_long %>% group_by(type, marker_type) %>% summarise(mean(value, na.rm = T))

# average 143 markers per lg input for map creation, 43.5 were hkhk, 47.1 lmll, 52 were nnnp
# after filtering and map creation, total 78 markers per lg per family, 18.4 hkhk, 29.4 lmll, 33.1 nnnp

cts <- marker_summary_long %>% 
  filter(marker_type != "nloc") %>%
  group_by(type, family) %>% summarise(sum(value, na.rm = T))

cts %>% dplyr::filter(type == "creation") %>% summarise(mean(`sum(value, na.rm = T)`)) # 3003 markers input per family
cts %>% filter(type == "mapping") %>% summarise(mean(`sum(value, na.rm = T)`)) # 1652 markers input per family

# but if you exclude the hkhks
marker_sum_exclude_hkhk <- marker_summary_long %>% 
  filter(marker_type != "nloc" & marker_type != "hkhk") %>%
  filter(type == "mapping") %>% 
  group_by(chr, family) %>%
  summarise(sum(value, na.rm = T))

mean(marker_sum_exclude_hkhk$`sum(value, na.rm = T)`)
head(marker_summary_long)

# which family had the most segregation distortion? i.e. the most markers removed
c <- marker_summary_long %>% 
  group_by(family, type) %>%
  filter(type == "creation") %>%
  summarise(sum(value, na.rm = T))

m <- marker_summary_long %>% 
  group_by(family, type) %>%
  filter(type == "mapping") %>%
  summarise(sum(value, na.rm = T))

cm <- left_join(c, m, by = "family")

colnames(cm)[3] <- c("creation")
colnames(cm)[5] <- c("mapping")

cm_sum <- cm %>% 
  mutate(dif =  creation - mapping) %>%
  mutate(percent_lost = dif / creation)

mean(cm_sum$percent_lost) # 45% of markers were lost to seg distortion

colnames(map) <- c("marker", "cM")

map_length <- map %>%  # we used map R2
  filter(cM >= 0.00) %>% 
  mutate(chr = substr(marker, 1, 5)) %>%
  group_by(chr) %>%
  summarise(max(cM))

sum(map_length$`max(cM)`) # 3385.849 cM total length
mean(map_length$`max(cM)`) # 161 cM length

# marker density
marker_density <- map %>% filter(cM >= 0.00) %>% 
  mutate(chr = substr(marker, 1, 5)) %>%
  group_by(chr) %>%
  tally()

marker_legnth_density <- left_join(marker_density, map_length, by = "chr")
mean(marker_legnth_density$n / marker_legnth_density$`max(cM)`)
# 0.93

#### Step 3: Calculating Number of Markers Shared with Existing Map ####

# read in consensus map 
con <- read.csv(paste(dir, script, "/data/ConsensusMap.csv", sep = ""), header = T) %>%
  dplyr::select(GBSV2RS, cM) %>% 
  dplyr::rename(CON_cM = cM)

map1 <- map %>% 
  mutate(GBSV2RS = str_replace(marker, "Chr", "S")) %>%
  dplyr::select(-marker)

# find markers in common
common <- inner_join(map1, con, by = "GBSV2RS")
per_lg <- common %>% 
  mutate(chr = as.numeric(substr(GBSV2RS, 2, 3))) %>%
  group_by(chr) %>% 
  tally()

min(per_lg$n)
max(per_lg$n)


#### Step 4: Calculate Segregation Distortion Across Linkage Groups ####
# first create a file for each family individually, by rbinding the original loc files

# create output
file_list <- list()

# loop
for (j in 1:length(famID)) {
  for (i in 1:21) {
    file <- read.table(paste(dir, "/Additional Filtering and JoinMap File Prep/output/", famID[j], "_chr", i, ".txt", sep = ""), skip = 5, fill = T)

    # remove individuals
    line <- which(grepl("individual", file$V1)) 
    file1 <- slice(file, 1:line-1)
    
    # save their names
    names <- slice(file, (line+1):nrow(file)) %>% dplyr::select(V1) # extract them 
    
    # put the file into a list
    file_list[[i]] <- file1
  }
  
  # do.call the family's chromosomes
  file_bind <- do.call("rbind", file_list)
  
  # write out the file
  file <- paste(dir, script, "/output/", famID[j], "_sd.loc", sep = "")
  write.table(paste("name = ", famID[j], sep = "" ), file = file, row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(paste("popt = CP"), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  write.table(paste("nloc = ", length(file_bind$V1), sep = ""), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  write.table(paste("nind = ", length(names$V1), sep = ""), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  cat("\n", file = file, append=TRUE)
  write.table(file_bind, file = file, quote = F, row.names = F, col.names = F, append = T)
  cat("\n", file = file, append=TRUE)
  write.table("individual names:", file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  write.table(names, file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
}

    
# load data into join map, calculate segregation distortion, export
# done

# set directory
dir_seg <- paste(dir, script, "/output/seg_distortion/", sep = "")
  
# output
fam_list <- list()
fam_sum_list <- list()

# loop
for (j in 1:length(famID)) {
  fam <- read.table(paste(dir_seg, famID[j], ".txt", sep = ""), na.strings = c("NA","."," "), header = T, fill = TRUE)
  
  # calculate pvalue because joinmap output sucks & create chromosome column
  fam <- fam %>% 
    dplyr::select(Locus, Segregation, X2, Df) %>% 
    mutate(pval = pchisq(X2, df = Df, lower.tail=FALSE),
           log_pval = -log10(pval),
           chr = as.numeric(substr(Locus, 4, 5)),
           pos = substr(Locus, 7, nchar(as.character(Locus))),
           sig = NA, 
           famID = famID[j], 
           Segregation = substr(Segregation, 2, 6))
           
  
  for (i in 1:nrow(fam)) {
    if (fam$pval[i] < 0.1) {fam$sig[i] <- "yes"}
    else if (fam$pval[i] > 0.1) {fam$sig[i] <- "no"}
    else {fam$sig[i] <- NA}
  }
  
  # output marker dataframe
  fam_list[[j]] <- fam
  
  # calculate marker summary
  fam_sum <- fam %>% group_by(Segregation, sig, chr) %>% tally()
  fam_sum1 <- fam_sum %>% mutate(famID = famID[j])
  fam_sum_list[[j]] <- fam_sum1
}

# marker summary
marker_summary <- do.call("rbind", fam_sum_list)
total <- marker_summary %>% group_by(famID) %>% summarise(sum(n))
sig_dist <- marker_summary %>% filter(sig == "yes") %>% group_by(famID) %>% summarise(sum(n))
no_sig_dist <- marker_summary %>% filter(sig == "no") %>% group_by(famID) %>% summarise(sum(n)) 

all <- cbind(total, sig_dist[,2], no_sig_dist[,2])
colnames(all) <- c("famID", "total", "sig", "not_sig")
mean(all$not_sig)

# average percent distorted
percent_distorted <- all %>% mutate(percent_distorted = sig / total) 
mean(percent_distorted$percent_distorted) # 0.40 average percent distorted across families
min(percent_distorted$percent_distorted) # 0.36 range of percent distorted
max(percent_distorted$percent_distorted) # 0.43

## what percent of the distorted were each kind of marker
sig_marker_type <- marker_summary %>% 
  filter(sig == "yes") %>% 
  group_by(famID, Segregation) %>% 
  summarise(sum(n)) %>%
  mutate(famID_seg = paste(famID, Segregation, sep = "_")) %>%
  as.data.frame()
colnames(sig_marker_type)[3] <- "yes"


not_marker_type <- marker_summary %>% 
  filter(sig == "no") %>% 
  group_by(famID, Segregation) %>% 
  summarise(sum(n)) %>%
  mutate(famID_seg = paste(famID, Segregation, sep = "_")) %>%
  as.data.frame() %>%
  dplyr::select(-famID, -Segregation)
colnames(not_marker_type)[1] <- "no"


seg_join <- left_join(not_marker_type, sig_marker_type, by = "famID_seg") %>%
  dplyr::select(famID, Segregation, no, yes)

colnames(seg_join) <- c("famID", "Segregation", "not_distorted", "distorted") 
sig_marker_summary <- seg_join %>% mutate(total = not_distorted + distorted) %>% mutate(percent_distorted = distorted / total)
sig_marker_summary %>% group_by(Segregation) %>% summarise(mean(percent_distorted))


# plot for just families
ggplot(marker_summary, aes(fill=sig, y=n, x=Segregation)) + 
  geom_bar(position="stack", stat="identity") + 
  facet_wrap(~famID, nrow =2) + 
  theme_bw() +
  ylab("Markers (ct)") + 
  xlab("Marker Type") + 
  scale_fill_manual(values = c("gray70", "gray40"), "Significant Segregation Distortion:") + 
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size=20, face = "bold"), 
        strip.text.x = element_text(size = 15),
        strip.background.x = element_rect(fill = "white", colour = NA), 
        strip.text.y = element_text(size = 15),
        strip.background.y = element_rect(fill = "white", colour = NA), 
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.position="bottom")

# save the figure
ggsave("seg_dist_plot.tiff", 
       plot = last_plot(), 
       path = paste(dir, script, "/output", sep = ""),
       device = "tiff",
       scale = 1, 
       width = 15, 
       height = 10, 
       units = "in")



#### Step 5: Determine the Number of Base Pairs per cM ####
# get the maximum cM and pos for each LG
map_sum <- map %>% 
  mutate(chr = as.numeric(substr(marker, 4, 5)), 
         pos = substr(marker, 7, nchar(marker))) %>%
  group_by(chr) %>% 
  summarise(max_cm = max(cM), max_bp = max(as.numeric(as.character(pos)))) %>%
  mutate(max_bp / max_cm)

hist(map_sum$`max_bp/max_cm`) 
mean(map_sum$`max_bp/max_cm`) # 2969015 # bp threshold for 1 cm
median(map_sum$`max_bp/max_cm`)

(2969015 * 21.08) #21.08 cM, or the median LD decay rate according to the analysis
# threshold for 21 cM 
bp_thresh <- 62586836 

62586836 / mean(map_sum$max_bp)
# this is about 13% of a chromosome length

       
       