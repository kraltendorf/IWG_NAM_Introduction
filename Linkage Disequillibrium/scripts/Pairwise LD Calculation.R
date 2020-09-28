# Project: IWG_NAM_Introduction
# Analysis - Pairwise LD Calculation
# Author: Kayla R. Altendorf
# Date: 02/02/2020

# load required packages
library("vcfR")
library("genetics")
library("dplyr")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("/Linkage Disequillibrium")

# declare where you want the output to go 
out_path <- paste(dir, script, "/output", sep = "")
in_path <- paste(dir, script, "/data", sep = "")

#### Step 1: Prepare Data ####
# read in imported genotype data 
vcf <- read.vcfR(paste(dir, "Variant Filtering and Imputation/output/NAM_GATK_imputed.vcf", sep = ""), convertNA = TRUE) # 8003 variants

# extract genotypes and id_frame
genotypes <- extract.gt(vcf, convertNA = FALSE)
genotypes[1:10, 1:10]
id_frame <- as.data.frame(vcf@fix[,1:5]) 
id_frame <- id_frame %>% mutate(ID = paste(CHROM, "_", POS, sep = ""))

# update sample names
# read in key
key <- read.table(paste(dir, "Variant Calling/data/new_key.txt", sep = ""), header = T)
key <- key %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  dplyr::select(GATK_Sample, Sample)

# prepare genotypes
genotypes1 <- t(genotypes) 
genotypes2 <- as.data.frame(genotypes1) %>% rownames_to_column(var = "GATK_Sample")

genotypes3 <- left_join(genotypes2, key, by = "GATK_Sample") %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-GATK_Sample) %>%
  t()

# change all | to / to remove phasing information
genotypes3[genotypes3=="0|1"] <- "0/1"
genotypes3[genotypes3=="1|0"] <- "0/1"
genotypes3[genotypes3=="1|1"] <- "1/1"
genotypes3[genotypes3=="0|0"] <- "0/0"

# change data into base pairs, because that's what the genetics package uses
genotypes5 <- genotypes3

for (i in 1:nrow(genotypes3)) {
  ref <- id_frame$REF[i]
  alt <- id_frame$ALT[i]
  call <- unlist(genotypes3[i,])
  call1 <- call
  call1[call == "0/1"] <- paste(alt, ref, sep = "/") # replace genotypes with their ref and alt alleles
  call1[call == "0/0"] <- paste(ref, ref, sep = "/")
  call1[call == "1/1"] <- paste(alt, alt, sep = "/")
  call1[call == "./."] <- NA
  genotypes5[i,] <- call1 # add edited row back into new dataframe
}

# filter to markers that are on the genetic consensus map
# load in the map
map <- read.table(paste(dir, "JoinMap/output/Map Files for MapQTL/Formatted for R/no_group.map", sep = ""),  header = F)
colnames(map) <- c("marker", "cM")

genotypes6 <- as.data.frame(genotypes5) %>% rownames_to_column("ID") %>% filter(ID %in% map$marker) %>% mutate(chr = substr(ID, 4, 5))
genotypes6[1:10, 1:10]


#### Step 2: Calculate LD ####
# create a vector for iterating through chromosomes
chroms_one_digit <- 1:21
chroms_two_digit <- str_pad(chroms_one_digit, width = 2, side = c("left"), pad = "0")

# create output dataframe 
markers_by_chrom <- list()

# set a vector for famID
famID <- c("WGN07", "WGN15", "WGN26", "WGN36", "WGN38", "WGN39", "WGN45", "WGN46", "WGN55", "WGN63")

# set working directory
setwd(out_path) 

# the master loop to calculate pairwise ld between markers on a per family basis
# WARNING: this takes a while to run, at least 5 minutes per family 
for (k in 1:length(famID)) {
  
  # set time
  start <- Sys.time()
  
  # prepare family data
  id_frame <- genotypes6 %>% dplyr::select(ID, chr) # take out marker ids and chom
  genotypes_fam <- genotypes6[, grepl(famID[k], colnames(genotypes6))] # take out one family at a time
  id_genotypes_fam <- cbind(id_frame, genotypes_fam)

  # then nest each family's chromosome into a list
  for (j in 1:length(chroms_two_digit)) {
  data <- id_genotypes_fam %>% filter(chr == chroms_two_digit[j]) %>% dplyr::select(-chr) %>% column_to_rownames("ID")
  markers_by_chrom[[j]] <- data
  }
  
  # empty out r2 output list
  r2 <- list()
  
  # calculate LD per chromosome per family using the 'genetics' package commands
  for (l in 1:length(markers_by_chrom)) {
    data <- makeGenotypes(t(markers_by_chrom[[l]]))
    ld.data <- LD(data)
    r2[[l]] <- ld.data$`R^2`
    
    end <- Sys.time()
    time <- end - start  
    print(paste("ld calc complete for ", famID[k], " LG ", l, " - elapsed time: ", round(time, 2) , sep = ""))
  }
  
  # empty out marker distance list
  markers_distance <- list()
  
  # format from matrix to dataframe and calculate distance
  for (m in 1:length(r2)) {
    r3 <- as.data.frame(r2[[m]]) %>% rownames_to_column("marker1")
    r4 <- r3 %>% pivot_longer(-marker1, names_to = "marker2", values_to = "r2") %>% filter(! is.na(r2))
    marker1 <- r4 %>% dplyr::select(marker1, r2)
    marker2 <- r4 %>% dplyr::select(marker2, r2)
    colnames(map)[1] <- "marker1"
    marker_1 <- left_join(marker1, map, by = "marker1")
    colnames(map)[1] <- "marker2"
    marker_2 <- left_join(marker2, map, by = "marker2")
    markers <- cbind(marker_1, marker_2)
    colnames(markers) <- c("marker1", "r2", "cM_1", "marker2", "r2_delete", "cM_2")
    markers1 <- markers %>% dplyr::select(-r2_delete) %>% mutate(distance_cm = abs(cM_1 - cM_2))
    markers_distance[[m]] <- markers1
  }
  
  # write out each dataframe separately
  for (p in 1:length(markers_distance)) {
    write.table(markers_distance[[p]], paste(famID[k], "_chr_", p, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
  }
  
  end <- Sys.time()
  time <- end - start
  print(paste("family ", famID[k], " ALL COMPLETE ", "elapsed time: ", round(time,3), sep = ""))
}


#### Step 3: Format the Final Dataset ####
# read in the data
ld_fam <- list()
fam_list <- list()
fam_out <- list()

for (i in 1:length(famID)) {
  files <- list.files(out_path, full.names = T)
  files_fam <- files[grepl(famID[i], files)] 
  fam <- lapply(files_fam, read.table, header = T)
  
  # make new cols - chr, distance_bp, and fam
  for (j in 1:length(fam)) {
    fam_out[[j]] <- fam[[j]] %>% 
      mutate(marker1_pos = substr(marker1, 7, nchar(as.character(marker1))),
             marker2_pos = substr(marker2, 7, nchar(as.character(marker2))),
             distance_bp = abs(as.numeric(marker1_pos) - as.numeric(marker2_pos)),
             chr = as.numeric(substr(marker1, 4, 5)),
             family = famID[i])
  }
  # rbind into one dataframe
  ld_fam[[i]] <- do.call("rbind", fam_out)
}

# rbind all data together
ld_all <- do.call("rbind", ld_fam)

# write out the final result so it doesn't have to be re-run! 
write.table(ld_all, paste(out_path, "/ld_by_family_by_chrom.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
