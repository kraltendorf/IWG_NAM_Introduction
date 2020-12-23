# Project: IWG_NAM_Introduction
# Analysis - GWAS
# Author: Kayla R. Altendorf
# Date: 12/1/2020

# load packages
library("sommer")
library("vcfR")
library("dplyr")
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



##### Step 2: Prepare SNP Data in Numeric Format #####
# prepare imputed genotype file in numeric format
vcf <- read.vcfR(paste(dir, "Variant Filtering and Imputation/output/NAM_GATK_imputed.vcf", sep = ""), convertNA = TRUE)

# extract genotypes and id_frame
genotypes <- extract.gt(vcf, convertNA = FALSE)
id_frame <- as.data.frame(vcf@fix[,1:5]) # extract a dataframe with SNP ids -- this will come in handy later
id_frame <- id_frame %>% mutate(ID = paste(CHROM, "_", POS, sep = ""))

# update sample names from their ID in variant calling (e.g. flowcell, lane, barcode)
# to their sample names
# read in key
key <- read.table(paste(dir, "Variant Calling/data/new_key.txt", sep = ""), header = T) %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  dplyr::select(GATK_Sample, Sample)

# prepare genotypes
genotypes1 <- t(genotypes) # transpose
genotypes2 <- as.data.frame(genotypes1) %>% rownames_to_column(var = "GATK_Sample")

genotypes3 <- left_join(genotypes2, key, by = "GATK_Sample") %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-GATK_Sample) %>%
  t()

# change all | to / to remove phasing information, if any
genotypes3[genotypes3=="0|1"] <- "0/1"
genotypes3[genotypes3=="1|0"] <- "0/1"
genotypes3[genotypes3=="1|1"] <- "1/1"
genotypes3[genotypes3=="0|0"] <- "0/0"

genotypes_numeric <- genotypes3

for (i in 1:nrow(genotypes3)) {
  gt <- unlist(genotypes3[i,])
  gt1 <- gt
  gt1[gt == "0/0"] <- 0
  gt1[gt == "0/1"] <- 1
  gt1[gt == "1/1"] <- 2
  gt1[gt == "./."] <- NA
  genotypes_numeric[i,] <- gt1 
}

genotypes_numeric <- as.data.frame(t(genotypes_numeric)) %>% rownames_to_column(var = "taxa")

# extract snp information as a separate file
mdp_SNP_information <- id_frame %>% 
  dplyr::select(ID, CHROM, POS) %>%
  dplyr::rename(Name = ID, 
         Chromosome = CHROM, 
         Position = POS) %>%
  mutate(Chromosome = as.numeric(substr(Chromosome, 4, 5)))

# write out result
write.table(genotypes_numeric, paste(dir, folder, "/output/GATK_NAM_snp_matrix_imputed.table.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(mdp_SNP_information, paste(dir, folder, "/output/mdp_SNP_information.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

# read it back in 
myGD <- read.table(paste(dir, folder, "/output/GATK_NAM_snp_matrix_imputed.table.txt", sep = ""), head = T) 
myGM <- read.table(paste(dir, folder, "/output/mdp_SNP_information.txt", sep = ""), head = T) 


#### Step 3: Run GAPIT ####
# edit here only
trait <- anthesis_score
trait_name <- c("anthesis_score")

# set directory
dir.create(paste(dir, folder, "/output/GAPIT/", trait_name, sep = ""))
setwd(paste(dir, folder, "/output/GAPIT/", trait_name, sep = ""))

# run
i = 4
for (i in 1:length(trait)) {
  phenotype <- trait[[i]]
  
  # choose only samples that are in common
  P <- phenotype$longID
  G <- myGD$taxa
  common <- Reduce(intersect, list(G, P))
  
  # filter phenotype data
  P2 <- filter(phenotype, longID %in% common) %>% droplevels() %>% arrange(longID)
  
  myGD2 <- myGD %>% 
    filter(taxa %in% common) %>%
    arrange(taxa)

  # make sure sample names align
  print(summary(myGD2$taxa == P2$longID))
  
  # rename phenotype header
  colnames(P2)[1] <- "Taxa"
  
  # rename phenotype to correct trait
  colnames(P2)[2] <- paste(trait_name, env[i], sep = "_")
  
  # gapit command
  myGAPIT <- GAPIT(
    Y=P2,
    GD=myGD2,
    GM=myGM,
    SNP.MAF=0.005, # inserting a MAF here since imputation introduced a few low freq alleles post vcftools filtering
    PCA.total = 0)
}
 


