Scripts for making normalized bar plots of on-target and untargeted integration

These scripts require the tnseq-pipeline (or the typeV version) to have already been run, and uses the read-locations.csv, Output_log and the fasta files for target DNA(genomic/plasmid) for the analysis . 

Notes before running:

Check the position of the end of spacer and the strandedness of the spacer relative to input fasta file and whether the sample was 'genomic'(for E.  coli integration) or 'plasmid' (for biochemical integration). A window of 100 bp downstream of spacer was selected as the on-target window. Reads falling anywhere else in the E. coli genome was noted as untargeted. Reads falling elsewhere on pTarget and pDonor was noted as untargeted for biochemical samples. Reads due to contaminating pDonor and a potential PCR recombination artifact (at position 1198 and 2328 nt) were masked on the pDonor. 

The reads corresposoning to integration at on-target and untargeted sites were normalalized to the total transposon-end containing reads for that sample and scaled to the sample with highest transposon-end containing reads. All samples in a set had the same read allocation and were from the same sequencing run.