#!/bin/bash -l
#PBS -l walltime=20:00:00,nodes=1:ppn=8,mem=22gb
#PBS -q batch
#PBS -N gather_vcfs
#PBS -e gather_vcfs.error
#PBS -o gather_vcfs.output
#PBS -A janderso
#PBS -W group_list=janderso
#PBS -M kaltendo@umn.edu
#PBS -m abe

# gather vcfs across all chromosomes

# go to the correct directory
cd /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/final_vcfs

# load packages
module load picard-tools
module load java

# GatherVcfs command using picard tools
java -Xmx15g -jar /panfs/roc/msisoft/picard/2.18.16/picard.jar GatherVcfs \
INPUT=final_chr01.vcf \
INPUT=final_chr02.vcf \
INPUT=final_chr03.vcf \
INPUT=final_chr04.vcf \
INPUT=final_chr05.vcf \
INPUT=final_chr06.vcf \
INPUT=final_chr07.vcf \
INPUT=final_chr08.vcf \
INPUT=final_chr09.vcf \
INPUT=final_chr10.vcf \
INPUT=final_chr11.vcf \
INPUT=final_chr12.vcf \
INPUT=final_chr13.vcf \
INPUT=final_chr14.vcf \
INPUT=final_chr15.vcf \
INPUT=final_chr16.vcf \
INPUT=final_chr17.vcf \
INPUT=final_chr18.vcf \
INPUT=final_chr19.vcf \
INPUT=final_chr20.vcf \
INPUT=final_chr21.vcf \
OUTPUT=NAM_GATK.vcf
