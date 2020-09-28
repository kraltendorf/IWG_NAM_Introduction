#!/bin/bash -l
#PBS -l walltime=05:00:00,nodes=1:ppn=8,mem=22gb
#PBS -q long
#PBS -N gatk_reference_prep
#PBS -e gatk_reference_prep.error
#PBS -o gatk_reference_prep.output
#PBS -A janderso
#PBS -W group_list=janderso
#PBS -M kaltendo@umn.edu
#PBS -m abe

# objective: create dictionary for the reference genome and index it with fadix

# load tools
module load gatk
module load samtools

# go to reference directory
cd /home/janderso/shared/IWG_v1_genome/annotated_v2_release/index

# gatk
gatk CreateSequenceDictionary -R Thinopyrum_intermedium.mainGenome.fasta

# samtools
samtools faidx Thinopyrum_intermedium.mainGenome.fasta

# end
