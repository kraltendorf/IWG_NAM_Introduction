#!/bin/bash


###########################################
## Shell script to run GBSv2 in TASSEL 5 ##
## Adapted from Narinder Singh          ##
###########################################

#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END --mail-user=jcrain@ksu.edu
#SBATCH --mem-per-cpu=32G   # Memory per core, use --mem= for memory per node
#SBATCH --time=0-23:00:00   # Use the form DD-HH:MM:SS
#SBATCH --job-name=UMN2
#SBATCH --partition=ksu-plantpath-jpoland.q,batch.q,killable.q

# Set JAVA VM Version
module load Java

java -version #print java version

#Change the following settings

#export PATH="$PATH:{path to }/bowtie2-2.2.6"
#like 
export PATH="$PATH:/homes/jcrain/Software/bowtie2-2.2.6"


# Also use update the email id to get notification about the job.
name=UMN2 #set job name
memory=30 #set memory size

# change to directory to working directory 
#cd {where GBS pipeline will be ran from}
#example using directory structure as provided  
cd /homes/jcrain/IWG_Work_Analysis/gbs/NAM_Parent_ID/beocat/gbs

#keyFile={path to key file}/data/Intermediate_File/UMN_Key2_Filtered.txt
#if following directory structure 
keyFile=../../data/Intermediate_File/UMN_Key2_Filtered.txt #will work

#seqDir={path to directory with FASTQ files}
#can be located anywhere
#like 
seqDir=/bulk/jpoland/sequence

#dbPath={path to directory with indexed genome version}
#can be located anywhere
#like 
dbPath=/homes/jcrain/IWG/FASTA_V1/build/IWG_V1

#setting variable name to run TASSEL pipeline
#tasselPath={path to TASSEL pipeline}/run_pipeline.pl
#like 
tasselPath=/homes/jcrain/Software/tassel/tassel-5-standalone/run_pipeline.pl

#set path to VCF tools for filtering
#vcfPath={path to VCF tools}/vcftools
#like 
vcfPath=/homes/jcrain/bin/vcftools

## NO NEED TO CHANGE ANYTHING FROM HERE ON ##
 
directory=`pwd`


#Step 1
#have to request memory but aslo specify it calling the plugin

## GBSSeqToTagDBPlugin  - RUN Tags to DB  
$tasselPath -Xms${memory}G -Xmx${memory}G -fork1 -GBSSeqToTagDBPlugin -e PstI-MspI \
    -i ${seqDir} \
    -db ${name}.db \
    -k ${keyFile} \
    -kmerLength 64 -minKmerL 50 -mnQS 20 -mxKmerNum 300000000 \
    -endPlugin -runfork1 >> ${name}_pipeline.out

#Step 2
## TagExportToFastqPlugin  
$tasselPath -fork1 -TagExportToFastqPlugin \
    -db ${name}.db \
    -o ${name}_tagsForAlign.fa.gz -c 10 \
    -endPlugin -runfork1 >> ${name}_pipeline.out
    
#Step 3
## RUN BOWTIE #-S is write to SAM file -U is unparied reads to be aligned -x is aligned files -p is parellel cores 2 hours,  # -l -60
bowtie2 -p 20 --end-to-end -D 20 -R 3 -N 0 -L 20 -i S,1,0.30 --n-ceil C,0,0 --score-min L,-6,0 -k 2 \
    -x $dbPath \
    -U ${name}_tagsForAlign.fa.gz \
    -S ${name}.sam >> ${name}_pipeline.out
    
#3.5 only get unique snps
grep -v "XS:i" ${name}.sam | awk '$4!=0' > ${name}_unique.sam

#Step 4  
## SAMToGBSdbPlugin - SAM to DB 
$tasselPath -Xms${memory}G -Xmx${memory}G -fork1 -SAMToGBSdbPlugin \
    -i ${name}_unique.sam \
    -db ${name}.db \
    -aProp 0.0 -aLen 0 \
    -endPlugin -runfork1 >> ${name}_pipeline.out

#Step 5
## DiscoverySNPCaller
$tasselPath -Xms${memory}G -Xmx${memory}G -fork1 -DiscoverySNPCallerPluginV2 \
    -db ${name}.db \
    -mnLCov 0.1 -mnMAF 0.01 -deleteOldData true \
     -endPlugin -runfork1 >> ${name}_pipeline.out
  
#Step 6  
## SNPQualityProfilerPlugin - RUN QUALITY PROFILER 30 minutes 20GB
$tasselPath -Xms${memory}G -Xmx${memory}G -fork1 -SNPQualityProfilerPlugin \
    -db ${name}.db \
    -statFile ${name}_SNPqual_stats.txt \
    -endPlugin -runfork1 >> ${name}_pipeline.out
  
#Step 7    
## UpdateSNPPositionQualityPlugin - UPDATE DATABASE WITH QUALITY SCORE fast < 30 minutes 15GB
$tasselPath -Xms${memory}G -Xmx${memory}G -fork1 -UpdateSNPPositionQualityPlugin \
    -db ${name}.db \
    -qsFile ${name}_SNPqual_stats.txt \
    -endPlugin -runfork1 >> ${name}_pipeline.out

#Ends SNP discovery with database
#Use Production SNP caller to get SNPs and filter
#Step 8    
## ProductionSNPCallerPluginV2 - RUN PRODUCTION PIPELINE - output .vcf
$tasselPath -Xms${memory}G -Xmx${memory}G -fork1 -ProductionSNPCallerPluginV2 \
    -db ${name}.db \
    -i ${seqDir} \
    -k ${keyFile} \
    -o ${name}.vcf \
    -e PstI-MspI -kmerLength 64 \
    -endPlugin -runfork1 >>  ${name}_pipeline.out 
    

#Step 9 Extract tag sequences using code from GetTagSequenceFromDBPlugin    
$tasselPath -Xms${memory}G -Xmx${memory}G -fork1 -GetTagSequenceFromDBPlugin \
    -db ${name}.db \
    -o ${name}.tags.txt \
    -endPlugin -runfork1

##Step 10 and 11 make output that is probably not used, but included anyway
##Step 10 Use code from Liang Gao to get Tag position and merge with tag sequence
## extract tag sequences in the db (sqlite3)
## extract tagid, snpid, table
## Note the location of the script GBSv2_table_join_to_get_tagid.sql might be different for you if you do not use KSU computer cluster

sqlite3 -separator $'\t' ${name}.db < \
../../scripts/gbs/GBSv2_table_join_to_get_tagid.sql \
    > ${name}.joined.table.txt

##Step 11 Code from Liang Gao
## using the python script to link ID, sequences together. 
## Note the location of the script snp_to_tagid.py might be different for you if you donot use KSU computer cluster
python ../../scripts/gbs/snp_to_tagid.py  \
    -t ${name}.tags.txt \
    -s  ${name}.joined.table.txt \
    > ${name}.snpid.linked.seq.fa


#Get allele depth counts
$vcfPath --vcf ${name}.vcf --out ${name} --geno-depth

#Run allele call by depth
perl ../../scripts/gbs/Genotype_by_Depth4.pl  ${name}.vcf > ${name}_Depth_Call.vcf 

#Filter by missing and maf
#Filter on minimum of two reads and 70% max missing, MAF > 0.01, no indels or multiallelic calls.
$vcfPath --vcf ${name}_Depth_Call.vcf  --out ${name}_Depth_Call_Filtered  --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --minDP 2 --max-missing 0.3 --recode --recode-INFO-all

#rename
mv ${name}_Depth_Call_Filtered.recode.vcf ${name}_Depth_Call_Filtered.vcf

#recall final depth
$vcfPath --vcf ${name}_Depth_Call_Filtered.vcf --out ${name}_Depth_Call_Filtered --geno-depth

## Convert to Hapmap format
$tasselPath -Xms${memory}G -Xmx${memory}G -fork1 -vcf ${name}_Depth_Call_Filtered.vcf -export ${name} -exportType Hapmap >>  ${name}_pipeline.out 

#zip original files
gzip ${name}.sam

#zip unique sam file
gzip  ${name}_unique.sam

#zip full vcf file
gzip  ${name}.vcf

#zip depth call
gzip  ${name}_Depth_Call.vcf

#zipfull depth call
gzip ${name}.gdepth

#zip filtered final vcf
gzip ${name}_Depth_Call_Filtered.vcf