# Project: IWG_NAM_Introduction
# Analysis - GWAS
# Author: Kayla R. Altendorf
# Date: 12/1/2020

# load packages
library("sommer")
library("vcfR")
library("dplyr")
install.packages("readr")
library("readr")
library("tibble")
library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("ape")
library("EMMREML")
library("compiler") #this library is already installed in R
library("scatterplot3d")
library("cowplot")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install(c("multtest"))

# download gapit
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

# location of the github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# folder we're on
folder <- "/GWAS"

#### Step 1: Load in Emmeans ####
# create a vector of traits
traits <- c("emergence_percent", "anthesis_score")
year <- c("2017", "2018", "2017", "2018")
loc <- c("STP", "STP", "TLI", "TLI")
env <- c("stp17", "stp18", "tli17", "tli18")


# read in backbone
backbone <- read.csv(paste(dir, "Phenotypic Data Analysis/Data/", "backbone.csv", sep = ""), header = T) %>% 
  dplyr::select(famID, parent, plantID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  dplyr::select(-plantID3, -famID) %>% 
  distinct()

# read in emmeans files for emergence_percent
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[1], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
emergence_percent <- list()
for (j in 1:length(files)) {
  emergence_percent[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

# read in emmeans files for anthesis_score
files <- list.files(paste(dir, "/Phenotypic Data Analysis/data/", traits[2], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
anthesis_score <- list()
for (j in 1:length(files)) {
  anthesis_score[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}



##### Step 2: Run GAPIT Using Imputed SNPs in  HapMap Format #####

# load in hapmap file
myG <- read_tsv(paste(dir, "Variant Filtering and Imputation/output/NAM_GATK_imputed.hmp.txt", sep = ""))

# update sample names
# read in key
key <- read.table(paste(dir, "Variant Calling/data/new_key.txt", sep = ""), header = T) %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  dplyr::select(GATK_Sample, Sample)

# remove extra cols
myG1 <- myG[,-1:-11]
myG1[1:10, 1:10]

# transpose
myG2 <- t(myG1) 
myG3 <- as.data.frame(myG2) %>% rownames_to_column(var = "GATK_Sample")
myG3[1:10, 1:10]

myG4 <- left_join(myG3, key, by = "GATK_Sample") %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-GATK_Sample) %>%
  t()

myG4[1:10, 1:10]



#### Step 3: Run GAPIT ####
# edit here only
trait <- emergence_percent
trait_name <- c("emergence_percent")

# set directory
dir.create(paste(dir, folder, "/output/GAPIT/", trait_name, sep = ""))
setwd(paste(dir, folder, "/output/GAPIT/", trait_name, sep = ""))

# run
for (i in 1:length(trait)) {
  phenotype <- trait[[i]]
  
  # choose only samples that are in common
  P <- phenotype$longID
  G <- colnames(myG4)
  common <- Reduce(intersect, list(G, P))
  
  # filter phenotype data
  P2 <- filter(phenotype, longID %in% common) %>% droplevels() %>% arrange(longID)
  
  # filter genotype data and put in order
  myG5 <- myG4[, colnames(myG4) %in% common]
  myG6 <- myG5[,order(colnames(myG5))] 
  
  # bind back unnecessary cols from hapmap
  myG7 <- cbind(myG[,1:11], myG6)
  
  # make sure column names align
  print(summary(colnames(myG7[,-1:-11]) == P2$longID))
  
  # get rid of hashes in column names
  colnames(myG7)[1] <- "rs"
  colnames(myG7)[6] <- "assembly"
  
  # rename phenotype header
  colnames(P2)[1] <- "Taxa"
  
  # rename phenotype to correct trait
  colnames(P2)[2] <- paste(trait_name, env_order[i], sep = "_")
  
  # make chromosome numeric
  myG7$chrom <- as.numeric(myG7$chrom)
  
  # write it out, and read it back in with header = F
  write.table(myG7, paste(dir, folder, "/output/myG7.hmp.txt", sep = ""), quote = F, row.names = F, sep = "\t") 
  myG8 <- read.table(paste(dir, folder, "/output/myG7.hmp.txt", sep = ""), head = F)

  # gapit command
  myGAPIT <- GAPIT(
    Y=P2,
    G=myG8,
    SNP.MAF=0.005,
    PCA.total = 0)
}


