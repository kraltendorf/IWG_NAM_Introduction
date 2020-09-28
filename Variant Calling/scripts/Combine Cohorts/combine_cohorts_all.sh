#!/bin/bash -l
#PBS -l walltime=50:00:00,nodes=1:ppn=8,mem=62gb
#PBS -q small
#PBS -N combine_cohorts_all
#PBS -e combine_cohorts_all.error
#PBS -o combine_cohorts_all.output
#PBS -A janderso
#PBS -W group_list=janderso
#PBS -M kaltendo@umn.edu
#PBS -m abe

# combine the cohort gvcfs

# find location of reference
REFERENCE=/home/janderso/shared/IWG_v1_genome/annotated_v2_release/index/Thinopyrum_intermedium.mainGenome.fasta

# set gatk settings
GATK_SETTINGS='-DF NotDuplicateReadFilter -DF MappingQualityAvailableReadFilter'

# from before when the files needed to be moved
# change to correct directory
#cd /home/janderso/kaltendo/NAM_GATK/GBarleyS/Pipeline/GATK_Pipeline/GenotypeGVCF/NAM_GATK/Raw_Variants

#mv combine_cohorts_1.g.vcf /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/all_cohorts
#mv combine_cohorts_2.g.vcf /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/all_cohorts
#mv combine_cohorts_3.g.vcf /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/all_cohorts
#mv combine_cohorts_4.g.vcf /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/all_cohorts

cd /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/all_cohorts

# load required tools
module load java
module load gatk

# run command
gatk --java-options "-Xmx50g" CombineGVCFs -R $REFERENCE \
--variant combine_cohorts_1.g.vcf \
--variant combine_cohorts_2.g.vcf \
--variant combine_cohorts_3.g.vcf \
--variant combine_cohorts_4.g.vcf ${GATK_SETTINGS} \
-O /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/all_cohorts/all_cohorts.g.vcf
