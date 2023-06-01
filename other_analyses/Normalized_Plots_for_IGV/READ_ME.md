Scripts for making normalized IGV plots for on-target and untargeted integration

These scripts require the tnseq-pipeline (or the typeV version) to have already been run, and uses the read-locations.csv, Output_log and the fasta files for target DNA(genomic/plasmid) for the analysis.

Notes before running:

The reads for target DNA coordinate (genome/plasmid) were normalized to the total transposon-end containing reads for that sample and scaled to the sample with highest transposon-end containing reads. All samples in a set had the same read allocation and were from the same sequencing run.