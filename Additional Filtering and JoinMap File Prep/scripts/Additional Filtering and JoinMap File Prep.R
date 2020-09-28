# Project: IWG_NAM_Introduction
# Script 2 - Additional Filtering and JoinMap File Prep 
# Author: Kayla Altendorf 
# Date: 6/15/2020

# load required packages
library("vcfR")
library("dplyr")
library("tibble")
library("readr")
library("devtools")
library("stringr")
library("rlang")

# location of github directory 
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("/Additional Filtering and JoinMap File Prep")


#### Step 1: Prepare Genotype Data ####
# read in vcf file
# this was prefiltered to a minor allele frequency of 0.005, or 0.05 / 10 families
# selfs and outcrosses were removed file name is NAM_GATK_filtered_maf_selfs.vcf
vcf <- read.vcfR(paste(dir, "Variant Filtering and Imputation/output/NAM_GATK_filtered_maf_selfs.vcf", sep = ""), convertNA = TRUE)
  
# extract genotypes and id_frame
genotypes <- extract.gt(vcf, convertNA = FALSE)
id_frame <- as.data.frame(vcf@fix[,1:5]) # extract a dataframe with SNP ids -- this will come in handy later
id_frame <- id_frame %>% mutate(ID = paste(CHROM, "_", POS, sep = ""))

# update sample names from their ID in variant calling (e.g. flowcell, lane, barcode)
# to their sample names
# read in key
key <- read.table(paste(dir, "Variant Calling/data/new_key.txt", sep = ""), header = T)

key <- key %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  dplyr::select(GATK_Sample, Sample)

# prepare genotypes
genotypes1 <- t(genotypes) # transpose
genotypes2 <- as.data.frame(genotypes1) %>% rownames_to_column(var = "GATK_Sample")

genotypes3 <- left_join(genotypes2, key, by = "GATK_Sample") %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-GATK_Sample) %>%
  t()

genotypes[1:10, 1:10] # check to see that the sample names have been replaced

# change all | to / to remove phasing information
genotypes3[genotypes3=="0|1"] <- "0/1"
genotypes3[genotypes3=="1|0"] <- "0/1"
genotypes3[genotypes3=="1|1"] <- "1/1"
genotypes3[genotypes3=="0|0"] <- "0/0"

# extract parents and progeny
genotypes_progeny <- genotypes3[, !grepl("P", colnames(genotypes3))] # remove parents
genotypes_parents <- genotypes3[, grepl("P", colnames(genotypes3))] # keep parents


#### Step 2: Get Mode Parental Genotypes ####
# since parents were replicated in the GBS process
# here's a way of taking the most frequent call across samples
# set mode function  
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# create dataframe for results of parental genotype modes
parentID <- c("WGN07P02", "WGN15P01", "WGN26P02", "WGN36P02", "WGN38P01", "WGN39P01", "WGN45P02", "WGN46P01", "WGN55P02", "WGN63P02", "WGN59P20")
genotypes_parents_mode <- data.frame(matrix(NA, nrow = nrow(genotypes_parents), ncol = 11)) # create empty dataframe with room for all parents
colnames(genotypes_parents_mode) <- parentID

# run through each parent, take mode, output into genotypes_parents_mode
for (j in 1:length(parentID)) {
  parents <- genotypes_parents[, grep(parentID[j], colnames(genotypes_parents))]
  for (i in 1:nrow(parents)){
    call <- unlist(parents[i,])
    mode <- getmode(call)
    genotypes_parents_mode[i,j] <- mode
  }
}

# add rownames to keep snp ids
rownames(genotypes_parents_mode) <- rownames(genotypes_parents)

#### Step 3: Make SNP Matrix for Pre-filtering for JoinMap #### 
# using the pre imputation dataset, combine the parental modes back with the prongey
genotypes_all <- cbind(genotypes_progeny, genotypes_parents_mode)
genotypes_all <- as.data.frame(lapply(genotypes_all, as.character), stringsAsFactors = FALSE)

# convert to marker matrix
genotypes_all_matrix <- genotypes_all

for (i in 1:nrow(genotypes_all)) {
  gt <- unlist(genotypes_all[i,])
  gt1 <- gt
  gt1[gt == "0/0"] <- 1
  gt1[gt == "0/1"] <- 0
  gt1[gt == "1/1"] <- -1
  gt1[gt == "./."] <- NA
  genotypes_all_matrix[i,] <- gt1 
}

genotypes_all_matrix_id <- cbind(id_frame, genotypes_all_matrix) # add the id_frame back
# write out result
write.table(genotypes_all_matrix_id, paste(dir, script, "/output/GATK_NAM_snp_matrix.table.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

            
#### Step 4: Filter for Missing Data and Minor Allele Frequency within Families ####
# first nest the families into a list of dataframes
# set a vector of family names
famID <- c("WGN07", "WGN15", "WGN26", "WGN36", "WGN38", "WGN39", "WGN45", "WGN46", "WGN55", "WGN63")

# create an output dataframe
families <- list()

# subset families
for (j in 1:length(famID)) {
  progeny <- genotypes_progeny[, grep(famID[j], colnames(genotypes_progeny))] # extract one family at a time
  donor <- genotypes_parents_mode[, grep(famID[j], colnames(genotypes_parents_mode))]
  common <- genotypes_parents_mode[, grep("WGN59", colnames(genotypes_parents_mode))]
  family <- cbind(common, donor, progeny)
  families[[j]] <- family
}

# iterate through each family in the list, conduct filtering
for (j in 1:length(famID)) {
  family <- as.data.frame(families[[j]])
  family <- family %>% rownames_to_column("ID")
  family$missingness <- rowSums(family[,-1:-3] == "./.") / (ncol(family)-3) # calculate % missingness while ignoring parents
  family <- filter(family, missingness < 0.20) # filter out those that are greater than 20%
  present <- rowSums(! family[,-1:-3] == "./.")*2 # count number of present (non missing) and their possible allele sites (2)
  hets <- rowSums(family[,-1:-3] == "0/1")  # count the hets
  homo <- rowSums(family[,-1:-3] == "1/1")*2 # count the alternate homozygotes
  freq <- hets + homo # get SNP allele freq
  family$maf <- freq / present
  family_filtered <- filter(family, between(maf, 0.05, 0.95)) # filter 
  family_filtered <- dplyr::select(family_filtered, -missingness, -maf) # remove filtering cols
  families[[j]] <- family_filtered 
  print(nrow(family_filtered)) # print the count to the console to keep status on progress
}

# this is extra, but is interesting nonetheless
# calculate the number of segregations, or how many are polymorphic within each family OR are NA
for (j in 1:length(families)) {
  parents <- families[[j]][,2:3]
  parents$result <- NA # create new col for result
  for (i in 1:nrow(parents)) {
    if (parents$common[i] != "./." & parents$donor[i] != "./.") {
    if (parents$common[i] == parents$donor[i]) {parents$result[i] <- "FALSE"}
      else if (parents$common[i] != parents$donor[i]) {parents$result[i] <- "TRUE"}
    }
  }
  seg_detect <- filter(parents, result == "TRUE")
  print(nrow(seg_detect))
}

segs <- c(2743, 2741, 2559, 2690, 2801, 2816, 2613, 2543, 3065, 2661)
mean(segs)


#### Step 5: Edit Impossible Genotypes ####
# first step, determine segregation type
# create some output lists
families_seg <- list()
seg_frame <- list()

for (j in 1:length(families)) {
  family <- families[[j]] # extract one family at a time
  seg_type <- NA # empty out this vector
  for (i in 1:nrow(family)) { # for each row
    common <- family$common[i] # extract the common and donor parent genotypes
    donor <- family$donor[i]
    if (common == "0/0" & donor == "0/1") {family$seg_type[i] <- "<nnxnp>"} # if the common is homozygous, donor is het, than the seg type is nnxnp
    else if (common == "1/1" & donor == "0/1") {family$seg_type[i] <- "<nnxnp>"} # and so on. 
    else if (common == "0/1" & donor == "0/1") {family$seg_type[i] <- "<hkxhk>"}
    else if (common == "0/1" & donor == "1/1") {family$seg_type[i] <- "<lmxll>"}
    else if (common == "0/1" & donor == "0/0") {family$seg_type[i] <- "<lmxll>"}
    else {family$seg_type[i] <- "UN"} # if none of these options are TRUE, then the locus is UNINFORMATIVE, or "UN" (e.g. homozygous alt: 1/1 x 1/1)
  }
  family_seg <- family %>% 
    filter(seg_type != "UN") # remove the uninformative markers
  seg_frame[[j]] <- family_seg %>% dplyr::select(ID, seg_type) # add seg_type to its own list 
  #family_seg <- family_seg %>% select(-seg_type)
  families_seg[[j]] <- family_seg
}

# next, calcuate the frequency of "wrong" genotypes
# create a list for output
genotype_freq <- list()

# calculate genotype frequency
for (j in 1:length(families_seg)) { # iterating through families
  family <- families_seg[[j]] %>% select(-seg_type) # remove the seg type col
  present <- rowSums(family[,-1:-3] != "./.")  # declare how many are present, or nonmissing
  hom_ref <- rowSums(family[,-1:-3] == "0/0") / present # calc freq of homzoygous ref, and so on
  het <- rowSums(family[,-1:-3] == "0/1") / present
  hom_alt <- rowSums(family[,-1:-3] == "1/1") / present
  genotype_freq[[j]] <- cbind(family[,1:3], select(families_seg[[j]], seg_type), hom_ref, het, hom_alt)
  genotype_freq[[j]]$filter <- NA
}

# then make edits accordingly
# create another list for the output
families_seg_ed <- list()

for (j in 1:length(families_seg)) { # iterating through the families
  family <- as.data.frame(lapply(families_seg[[j]], as.character), stringsAsFactors = FALSE) %>% select(-seg_type) # change to as.character to avoid errors
  for (i in 1:nrow(family)) {
    common <- family$common[i] # extract common and donor parents
    donor <- family$donor[i]
    if (common == "0/0" & donor == "0/1") { # if common parent is homozygous, donor parent is heterozygous
      
      # and if the frequency of homozygous alternate alelle, which is not possible at this locus,
      # is really low, say between 0 and 0.05, we can say with pretty strong confidence that 
      # it's likely an error and that since the alternate alelle was picked up, it's likely a het
      if (genotype_freq[[j]]$hom_alt[i] != 0 & genotype_freq[[j]]$hom_alt[i] < 0.05) { 
        gt <- unlist(family[i,-1:-3])
        gt1 <- gt
        gt1[gt == "1/1"] <- c("0/1") # replace this site with het
        family[i, -1:-3] <- gt1
      }
      else if (genotype_freq[[j]]$hom_alt[i] != 0 & genotype_freq[[j]]$hom_alt[i] > 0.05) {
        
        # BUT if its erroneous at greater than 5% of the time, then drop it
        # because then it's more likely that the parent calls were wrong
        genotype_freq[[j]]$filter[i] <- "remove" 
      }
    }
    
    # moving on through the rest of the possible outcomes
    else if (common == "1/1" & donor == "0/1") {
      if (genotype_freq[[j]]$hom_ref[i] != 0 & genotype_freq[[j]]$hom_ref[i] < 0.05) {
        gt <- unlist(family[i,-1:-3])
        gt1 <- gt
        gt1[gt == "0/0"] <- c("0/1")
        family[i, -1:-3] <- gt1
      }
      else if (genotype_freq[[j]]$hom_ref[i] != 0 & genotype_freq[[j]]$hom_ref[i] > 0.05) {
        genotype_freq[[j]]$filter[i] <- "remove"
      }
    }
    else if (common == "0/1" & donor == "0/0") {
      if (genotype_freq[[j]]$hom_alt[i] != 0 & genotype_freq[[j]]$hom_alt[i] < 0.05) {
        gt <- unlist(family[i,-1:-3])
        gt1 <- gt
        gt1[gt == "1/1"] <- c("0/1")
        family[i, -1:-3] <- gt1
      }
      else if (genotype_freq[[j]]$hom_alt[i] != 0 & genotype_freq[[j]]$hom_alt[i] > 0.05) {
        genotype_freq[[j]]$filter[i] <- "remove"
      }
    }
    else if (common == "0/1" & donor == "1/1") {
      if (genotype_freq[[j]]$hom_ref[i] != 0 & genotype_freq[[j]]$hom_ref[i] < 0.05) {
        gt <- unlist(family[i,-1:-3])
        gt1 <- gt
        gt1[gt == "0/0"] <- c("0/1")
        family[i, -1:-3] <- gt1
      }
      else if (genotype_freq[[j]]$hom_ref[i] != 0 & genotype_freq[[j]]$hom_ref[i] > 0.05) {
        genotype_freq[[j]]$filter[i] <- "remove"
    }
    }
  }
  family_ed <- left_join(family, genotype_freq[[j]], by = "ID") %>% filter(is.na(filter)) %>% select(-common.y, -donor.y, seg_type, -hom_ref, -het, -hom_alt, -filter)
  families_seg_ed[[j]] <- family_ed
}


# now look at how things have changed
genotype_freq_ed <- list()

# first calculate genotype frequency of the edited file
for (j in 1:length(families_seg_ed)) {
  family <- families_seg_ed[[j]] %>% select(-seg_type)
  present <- rowSums(family[,-1:-3] != "./.") 
  hom_ref <- rowSums(family[,-1:-3] == "0/0") / present
  het <- rowSums(family[,-1:-3] == "0/1") / present
  hom_alt <- rowSums(family[,-1:-3] == "1/1") / present
  genotype_freq_ed[[j]] <- cbind(families_seg_ed[[j]][,1:3], select(families_seg_ed[[j]], seg_type), hom_ref, het, hom_alt)
}



#### Step 6: Change to JoinMap Format ####
# create a list for output
families_seg_ed_jm <- list()

# iterate through each edited family
# determine the seg type
# change formats to avoid errors

for (j in 1:length(families_seg_ed)) {
  seg_type <- as.data.frame(lapply(families_seg_ed[[j]], as.character), stringsAsFactors = FALSE) %>% select(ID, seg_type)
  family <- as.data.frame(lapply(families_seg_ed[[j]], as.character), stringsAsFactors = FALSE) %>% select(-seg_type)
  for (i in 1:nrow(family)) {
    if (seg_type$seg_type[i] == "<nnxnp>") { # if it's this seg type
      gt <- unlist(family[i, -1:-3]) # extract the genotypes
      gt1 <- sapply(X = gt, FUN = function(gen){ # apply this function to replace them with their joinmap codes
        if (gen == "0/0"){"nn"} 
        else if (gen == "1/1"){"nn"}
        else if (gen == "0/1"){"np"}
        else {"--"} # this is missing
      })
      family[i, -1:-3] <- gt1
    }
    else if (seg_type$seg_type[i] == "<hkxhk>") { # then the next seg type... 
      gt <- unlist(family[i, -1:-3])
      gt1 <- sapply(X = gt, FUN = function(gen){
        if (gen == "0/0"){"hh"} 
        else if (gen == "1/1"){"kk"}
        else if (gen == "0/1"){"hk"}
        else {"--"}
      })
      family[i, -1:-3] <- gt1
    }
    else if (seg_type$seg_type[i] == "<lmxll>") { # and the last seg type... 
      gt <- unlist(family[i, -1:-3])
      gt1 <- sapply(X = gt, FUN = function(gen){
        if (gen == "0/0"){"ll"} 
        else if (gen == "1/1"){"ll"}
        else if (gen == "0/1"){"lm"}
        else {"--"}
      })
      family[i, -1:-3] <- gt1
    }
  }
  family_progeny_only <- family %>% select(-common.x, -donor.x)
  family_progeny_only_ids <- family_progeny_only %>% column_to_rownames("ID")
  family_progeny_only_ids_sorted <- family_progeny_only_ids[ , order(names(family_progeny_only_ids))]
  family_progeny_only_ids_sorted_id <- family_progeny_only_ids_sorted %>% rownames_to_column("ID")
  family_progeny_only_seg <- left_join(seg_type, family_progeny_only_ids_sorted_id, by = "ID") %>% mutate(chrom = as.integer(substr(ID, 4, 5)))
  family_progeny_only_seg_rownames <- family_progeny_only_seg %>% column_to_rownames("ID")
  families_seg_ed_jm[[j]] <- family_progeny_only_seg_rownames
}




#### Step 7: Write out Text Files #####
setwd(paste(dir, script, "/output", sep = ""))

# write out family names
for (j in 1:length(families_seg_ed_jm)) {
  family <- as.data.frame(families_seg_ed_jm[[j]]) %>% rownames_to_column("ID")
  family_abridged <- family %>% select(-ID, -seg_type, -chrom)
  family_names <- colnames(family_abridged)
  write.table(family_names, paste(dir, script, "/output/", famID[j], "_names.txt", sep = ""), row.names = F, quote = F, col.names = F, sep = "\t")
}

# write out the joinmap files
for (j in 1:length(families_seg_ed_jm)) {
  family <- as.data.frame(families_seg_ed_jm[[j]]) %>% rownames_to_column("ID") # make the rownames a column with the SNP ID
  family_abridged <- family %>% select(-ID, -seg_type, -chrom) # remove extra cols
  family_names <- colnames(family_abridged) # get family names
  for (i in 1:21){ # and write out a file for each chromosome
    chr <- filter(family, chrom == i) %>% select(-chrom)
    filename <- paste(famID[j], "_chr", i, ".txt", sep = "") # determine file name
    write.table(paste("name = ", famID[j], "_chr", i, sep = "" ), # add each component of a joinmap file, appending along the way
              filename, row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(paste("popt = CP"), # cross pollinated population type
              filename, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    write.table(paste("nloc = ", nrow(chr), sep = ""), # the number of loci, which is the number of rows
              filename, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    write.table(paste("nind = ", (ncol(chr)-2), sep = ""), # number of individuals
              filename, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file= filename, append=TRUE) # a new line
    write.table(chr, filename, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file= filename, append=TRUE)
    write.table("individual names:", filename, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    write.table(colnames(chr)[-1:-2], filename, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  }
}

     
     