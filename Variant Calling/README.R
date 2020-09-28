#read me for github
#variant calling pipeline

##
Variant Calling Pipeline Using Genome Analysis Toolkit (GATK) v.4.1.2
This variant calling pipeline was developed by @neyhartj and described in detail on his Github page. Changes were made to adapt to the latest version and for intermediate wheatgrass.

The order of operations includes:
  
gatk_reference_prep.sh
run_fastx_barcode_splitter.sh
run_post_demultiplex_cleanup.sh
run_assess_quality.sh
run_quality_control.sh
run_assess_quality.sh
run_concat_command.sh # this was a custom script, not part of GBarleyS, that was used to combine reads that were associated with two different barcodes but those that originated from the same sample.
run_read_mapping.sh
run_samptools_processing.sh
run_haplotype_caller.sh
combine_cohorts (see folder)
run_genotype_gvcfs.sh
vcftools_filtering.sh
The folder also contains the original key file (gatk_key.txt) and the new_key.txt, which was created to concatenate reads from the same sample using the run_concat_command.sh script.