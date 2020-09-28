#!/bin/bash -l
#PBS -l walltime=04:00:00,nodes=1:ppn=2,mem=62gb
#PBS -q small
#PBS -N concat_fastq
#PBS -e concat_fastq.error
#PBS -o concat_fastq.output
#PBS -A janderso
#PBS -W group_list=janderso
#PBS -M kaltendo@umn.edu
#PBS -m abe

# go to the base dir
cd /scratch.global/kaltendo/gatk_temp/NAM_GATK/Quality_Control/Quality_Control_Fastq

# sort the barcode info file by sample name and place in working directory
sort -k4 /home/janderso/kaltendo/NAM_GATK/new_key.txt \
| tail -n +2 > /scratch.global/kaltendo/gatk_temp/NAM_GATK/Quality_Control/Quality_Control_Fastq/file_list.txt

# Use cut to put columns into separate arrays 
flowcell_array=($(cut -f1 file_list.txt))
lane_array=($(cut -f2 file_list.txt))
barcode_array=($(cut -f5 file_list.txt))

# loop over the index of both arrays
for i in ${!flowcell_array[@]}; do
	num=$(( $i % 2 ))
	if [ $num -eq 0 ]; then
	
		# get the ith elements of each array
		flowcell=${flowcell_array[$i]}
		lane=${lane_array[$i]}
		barcode=${barcode_array[$i]}
	
		# get the ith + 1 elements of each array
		flowcell2=${flowcell_array[$((i+1))]}
		lane2=${lane_array[$((i+1))]}
		barcode2=${barcode_array[$((i+1))]}
			
		# name the files 
		file1=/scratch.global/kaltendo/gatk_temp/NAM_GATK/Quality_Control/Quality_Control_Fastq/${flowcell}_${lane}_${barcode}_QC.fastq.gz 
		file2=/scratch.global/kaltendo/gatk_temp/NAM_GATK/Quality_Control/Quality_Control_Fastq/${flowcell2}_${lane2}_${barcode2}_QC.fastq.gz
		
		# check if they both exists
		if [ -e $file1 -a -e $file2 ]; then

				# concat the files, remove the old copies and change the name to the original then 
				cat "$file1" "$file2" > file_cat_QC.fastq.gz
				rm $file1
				rm $file2
				mv file_cat_QC.fastq.gz $file1
		fi
	fi
done
