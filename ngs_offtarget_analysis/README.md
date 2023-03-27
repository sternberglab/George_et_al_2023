# Scripts for analyzing NGS read locations

### These scripts require the tnseq-pipeline (or the typeV version) to have already been run, and uses the read-locations.csv and input files from that pipeline to analyze off-target reads. 

# Notes before running:
Check the `parameters.py` file and comments there to specify your data source and other values that can be adjusted. 

In your input.csv file, if you do NOT provide a genbank_file in parameters.py, the target FASTA file paths should be absolute (eg. `/home/users/Chris/genome.fasta`) and not relative (`./genome.fasta`), since these scripts will be run from a different directory than the original pipeline. If a genbank file is provided it doesn't matter since the input fasta file is not used. A genbank file is required for essential gene analysis and for removing insertions into T7RNAP from other analyses. 


### Off-target Cascade binding
- `alternate_targets.py` finds the most likely Cascade binding site for a given read location in the genome. It evaluates 20 targets upstream, and 20 targets downstream of the read, centered around the expected protospacer location. The protospacer to intergration site canonical distance is customizable. 
- The likely site is determined by the highest match in the seed region (first 8 nucleotides), and, in case of ties, highest match for the entire spacer. There is an option to exclude every 5th basepair, if the system being examined has them flipped out. 


### AT enrichment
- `at_enrichment.py` contains functions for evaluating AT enrichemnt. 
-- `get_at_test_data` creates a baseline comparison dataset for a given genome and sequence length. It randomly picks a point in the genome and gets the highest AT percentage either upstream or downstream. Ex. For a sequence length of 30, if location 500 was randomly chosen, the AT percentage from 470-500 and 500-530 are calculated, and the higher value is saved. This process is repeated many times (set by `AT_PCT_TEST_SIZE`) and the results are saved. 


### Correlate off-targets
- `correlate_offtargets.py` contains functions to find common clusters of off-target reads across samples. The `score` is a custom metric combining the number of reads at the given bin with the number of samples containing off-target reads in that bin. 
