#!/bin/bash -l
#PBS -l walltime=20:00:00,nodes=1:ppn=3,mem=62gb
#PBS -q small
#PBS -N combine_cohorts_1
#PBS -e combine_cohorts_1.error
#PBS -o combine_cohorts_1.output
#PBS -A janderso
#PBS -W group_list=janderso
#PBS -M kaltendo@umn.edu
#PBS -m abe

# combine the cohort gvcfs

# find location of reference
REFERENCE=/home/janderso/shared/IWG_v1_genome/annotated_v2_release/index/Thinopyrum_intermedium.mainGenome.fasta

# set gatk settings
GATK_SETTINGS='-DF NotDuplicateReadFilter -DF MappingQualityAvailableReadFilter'

# change to correct directory
cd /scratch.global/kaltendo/gatk_temp/NAM_GATK/HaplotypeCaller/NAM_GATK/GenotypeGVCF_attempt2/NAM_GATK/Raw/Combined_GVCFs

# load required tools
module load java
module load gatk

# run command
gatk --java-options "-Xmx50g" CombineGVCFs -R $REFERENCE \
--variant CBEEGANXX_1_cohort.g.vcf \
--variant CBEEGANXX_2_cohort.g.vcf \
--variant CBEEGANXX_3_cohort.g.vcf \
--variant CBEEGANXX_4_cohort.g.vcf ${GATK_SETTINGS} \
-O /home/janderso/kaltendo/NAM_GATK/GBarleyS/Pipeline/GATK_Pipeline/GenotypeGVCF/NAM_GATK/Raw_Variants/combine_cohorts_1.g.vcf
