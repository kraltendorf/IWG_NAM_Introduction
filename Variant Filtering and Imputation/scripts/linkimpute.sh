#PBS -l walltime=20:00:00,nodes=1:ppn=2,mem=16gb
#PBS -q small
#PBS -N linkimpute_imputation
#PBS -e linkimpute_imputation.error
#PBS -o linkimpute_imputation.output
#PBS -A janderso
#PBS -W group_list=janderso
#PBS -M kaltendo@umn.edu
#PBS -m abe

# run linkimpute

# Set this variable to folder containing vcf files to impute
cd /home/janderso/kaltendo/NAM_GATK/vcf_filtering

# load java
module load java/openjdk-11.0.2

# call linkimpute to impute on all chromsomes at once
java -jar /home/janderso/kaltendo/software/LinkImpute.jar -v NAM_GATK_filtered_maf_selfs_progeny.vcf NAM_GATK_imputed.vcf