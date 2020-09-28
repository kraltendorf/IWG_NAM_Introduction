#!/bin/bash -l
#PBS -l walltime=96:00:00,nodes=1:ppn=20,mem=5gb
#PBS -q small
#PBS -N structure2
#PBS -e structure2.error
#PBS -o structure2.output
#PBS -t 6,8
#PBS -A janderso
#PBS -W group_list=janderso
#PBS -M kaltendo@umn.edu
#PBS -m abe

### shell script to run Structure program

module load java
module load structure/2.3.4-console-CentOS7

# set array

arrayfile=$PBS_ARRAYID

cd /home/janderso/kaltendo/structure

# 5 is the number of replicates at each value of K.

for rep in {1..10} 
do 
structure -m mainparams -K ${arrayfile} -o /home/janderso/kaltendo/structure/outfiles_no_parents/outfile_k${arrayfile}_rep${rep}
done