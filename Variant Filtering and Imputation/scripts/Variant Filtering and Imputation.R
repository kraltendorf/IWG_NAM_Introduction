# Project: IWG_NAM_Introduction
# Analysis - Variant Filtering and Imputation
# Author: Kayla Altendorf 
# Date: 6/15/2020

# load required packages
library("dplyr")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("Variant Filtering and Imputation")

#### Step 1: Filter Using VCFtools ####
# biallelic
# minumum allele depth of 5
# maximum missing data of 0.8
# minor allele frequency of 0.005

./vcftools --vcf /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/final_vcfs/NAM_GATK.vcf --min-alleles 2 --max-alleles 2 --minDP 5 --max-missing 0.8 --recode --recode-INFO-all --out /home/janderso/kaltendo/NAM_GATK_filtered.vcf
./vcftools --vcf /home/janderso/kaltendo/NAM_GATK/vcf_filtering/NAM_GATK_filtered.vcf --maf 0.005 --recode --recode-INFO-all --out /home/janderso/kaltendo/NAM_GATK_filtered_maf.vcf

# recalculate % missing per individual in vcf tools using maf filtered dataset
./vcftools --vcf /home/janderso/kaltendo/NAM_GATK/vcf_filtering/NAM_GATK_filtered_maf.vcf --missing-indv

# download file
missingness <- read.table(paste(dir, script, "/data/out.imiss", sep = ""), header = T) %>% arrange(-F_MISS)

# after viewing the file, these two have over 70 missing data
to_remove <- data.frame(GATK_Sample = c("CBEEGANXX_1_BC137_", "CBEEGANXX_3_BC4_"))


#### Step 2: Remove Selfs and Outcroses ####
# read in key file 
key <- read.table(paste(dir, "/Variant Calling/data/new_key.txt", sep = ""), header = T)
key <- key %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  select(GATK_Sample, Sample)

# read in selfs and unintended outcrosses - identified in a previous analysis
selfs <- read.table(paste(dir, script, "/data/selfs_outcrosses.txt", sep = ""), header = F)
colnames(selfs)[1] <- "Sample"

# join the two dataframes together
to_remove2 <- left_join(selfs, key, by = "Sample")
write.table(to_remove$GATK_Sample, paste(dir, script, "/output/samples_to_remove.txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")

# append 
write.table(to_remove2$GATK_Sample, paste(dir, script, "/output/samples_to_remove.txt", sep = ""), append = T, col.names = F, row.names = F, quote = F, sep = "\t")

# identify the parents, too, so they can be removed for imputation 
parents <- dplyr::filter(key, grepl("P", Sample)) # if the ID has a "P" in it, that refers to PARENT

# write out this too
write.table(parents$GATK_Sample, paste(dir, script, "/output/parents.txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")

# remove selfs and individuals with high missing data
# here is the vcf script, submitted as a job through terminal
./vcftools --vcf /home/janderso/kaltendo/NAM_GATK/vcf_filtering/NAM_GATK_filtered_maf.vcf \
--remove /home/janderso/kaltendo/NAM_GATK/vcf_filtering/samples_to_remove.txt \
--recode --recode-INFO-all --out /home/janderso/kaltendo/NAM_GATK/vcf_filtering/NAM_GATK_filtered_maf_selfs.vcf

# also remove the parents because we don't include them in imputation
./vcftools --vcf /home/janderso/kaltendo/NAM_GATK/vcf_filtering/NAM_GATK_filtered_maf_selfs.vcf \
--remove /home/janderso/kaltendo/NAM_GATK/vcf_filtering/parents.txt \
--recode --recode-INFO-all --out /home/janderso/kaltendo/NAM_GATK/vcf_filtering/NAM_GATK_filtered_maf_selfs_progeny.vcf

#### Step 3: Impute using LinkImpute ####
# see the script "linkimpute.sh"

#### Step 4: Finalize the Imputed File ####
# after imputation, it's required to change the header using bcftools since linkimpute messes it up
# so we'll extract the header from the pre-imputation file
bcftools view -h NAM_GATK_filtered_maf_selfs_progeny.vcf -o header.txt

# edit imputed file and put first line back so it knows it's a vcf file then reheader
# use vim to read the file, press "i" to make the edit, then esc :wq
bcftools reheader -h header.txt NAM_GATK_imputed.vcf -o imputed_new_header.vcf

# change format to hapmap file using Tassel
/home/janderso/kaltendo/tassel-5-standalone/run_pipeline.pl -Xms4G -Xmx4G -fork1 \
-vcf /home/janderso/kaltendo/NAM_GATK/vcf_filtering/NAM_GATK_imputed_header.vcf \
-export /home/janderso/kaltendo/NAM_GATK/vcf_filtering/NAM_GATK_imputed -exportType Hapmap

