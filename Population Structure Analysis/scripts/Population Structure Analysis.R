# Project: IWG_NAM_Introduction
# Analysis - Population Structure Analysis
# Author: Kayla Altendorf 
# Date: 6/15/2020

# load necessary packages 
library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("broom")
library("rrBLUP")
library("tibble")
library("data.table") 
library("stringr")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("/Population Structure Analysis")

# set path for input files
in_path <- paste(dir, script, "/data", sep = "")
out_path <- paste(dir, script, "/output", sep = "")


#### Step 1: Read in Genotype Data ####
# read in GATK_NAM_snp_matrix.table.txt resulting from Additional Filtering and JoinMap File Prep Script
snp_matrix <- read.table(paste(dir, "/Additional Filtering and JoinMap File Prep/output/GATK_NAM_snp_matrix.table.txt", sep = ""), header = T)
snp_matrix[1:10, 1:10]

#### Step 2: Calculate Additive Relationship Matrix ####
# prepare matrix
snp_matrix_no_id <- snp_matrix[,-1:-5] # remove the SNP IDs
snp_matrix_no_id_t <- t(snp_matrix_no_id) # transpose the data

# using the A.mat function from rrBLUP
K <- A.mat(X = snp_matrix_no_id_t, min.MAF = 0, max.missing = 1) # default is 0.05 # no selection!

# run a PCA
pK <- prcomp(x = K)
summ <- summary(pK)
indiv_pc <- tidy(pK)

indiv_pc$famID <- substr(indiv_pc$row, 1, 5)
indiv_pc$type <- substr(indiv_pc$row, 6, 6)
indiv_pc$PC <- paste0("PC", indiv_pc$PC)
indiv_pc$donor <- indiv_pc$type == "P"

# subset the data
indiv_pc <- subset(indiv_pc, PC %in% c("PC1", "PC2"))
indiv_pc <- spread(data = indiv_pc, key = PC, value = value)

# rename common parent for figure
indiv_pc <- indiv_pc %>% mutate(famID = replace(famID, famID == "WGN59", "Common Parent"))
colnames(indiv_pc)[1] <- "data_id"
indiv_pc <- as.data.frame(indiv_pc)

# PCA plot 
pca_plot <- ggplot(data = indiv_pc, aes(x = PC1, y = PC2, color = famID, size = donor)) + 
  geom_point() + 
  theme_minimal() + 
  ylab("PC2 (15.5%)") + 
  xlab("PC1 (22.2%)") + 
  theme(axis.text.x = element_text(size = 25)) +
  theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.title.x = element_text(size = 30)) + 
  theme(axis.title.y = element_text(size = 30)) +
  theme(legend.text = element_text(size = 25)) + 
  theme(legend.title = element_text(size = 30)) + 
  guides(color=guide_legend(override.aes = list(size=4), title="Family")) + 
  guides(size = FALSE) + 
  scale_color_manual(values=c("#616365", "#c1d42f", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#770026",  "#68bbaf", "#D55E00", "#0b3d4c", "#CC79A7")) 

# save the figure
ggsave("pca_plot.tiff", 
       plot = last_plot(), 
       path = paste(dir, script, "/output", sep = ""),
       device = "tiff",
       scale = 1, 
       width = 15, 
       height = 10, 
       units = "in")
