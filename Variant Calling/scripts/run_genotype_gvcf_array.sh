#!/bin/bash -l
#PBS -l walltime=20:00:00,nodes=1:ppn=8,mem=22gb
#PBS -q batch
#PBS -N genotypegvcfs_array
#PBS -e genotypegvcfs_array.error
#PBS -o genotypegvcfs_array.output
#PBS -t 1-21
#PBS -A janderso
#PBS -W group_list=janderso
#PBS -M kaltendo@umn.edu
#PBS -m abe

# array job script for genotypegvcfs on a per chromosome basis

# load packages
module load java
module load gatk

# name reference
REFERENCE=/home/janderso/shared/IWG_v1_genome/annotated_v2_release/index/Thinopyrum_intermedium.mainGenome.fasta

# extract array number
array_num=$PBS_ARRAYID

# make it two digits
if [[ ${#array_num} -lt 2 ]] ; then
    array_num="00${array_num}"
    array_num="${array_num: -2}"
fi

# run command
gatk --java-options "-Xmx15g" GenotypeGVCFs \
-R $REFERENCE \
-V /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/all_cohorts/all_cohorts.g.vcf \
-DF NotDuplicateReadFilter \
-L Chr${array_num} \
-O /scratch.global/kaltendo/gatk_temp/NAM_GATK/GenotypeGVCF/final_vcfs/final_chr${array_num}.vcf 

