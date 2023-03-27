# Which analyses to run
# 'offtarget': each sample will have it's reads analyzed
# for Cascade off-target binding and AT richness
# 'correlation': only look for correlations in off-target reads
# between the listed samples
# 'both': both of the above
analysis_type = 'correlation'


# The directory where the _read_locations.csv files are located
# If unmoved, they will be in '{your tnseq-pipeline directory}/outputs/samples'
pipeline_outputs_directory = '../tnseq-pipeline/outputs/samples'

# The input information csv with information correlating 
# to the samples and the outputs produced
# This is the same file used to run Illumina-pipeline
info_file = ''

# If a genbank file is provided, then the target fasta 
# from the info file is ignored and this is used instead
# Genbank files are required for analyzing whether locations
# are intergenic, or entering which genes
# To leave blank, genbank_file = ''
# This will be used for ALL samples in the run, not the input 
# target fasta. This is also required to exclude T7 RNAP reads
# from AT weighted average enrichment calculations
genbank_file = ''

# The sample IDs (from the info_file) to analyze
# Leave empty (samples = []) to run on all samples
# in the info_file
samples = []

# Output files will be created here, in 'outputs', and under
# a subdirectory with this name. Change it to keep different
# results separate. 
output_subdirectory = 'ngs_results_1'

# The bp length around a site to look for AT enrichment: 
# More windows and larger windows will run slower
AT_enrichment_sizes = [40, 100]

# Size of target window in bp. Don't analyze potential 
# offtarget Cascade binding (mismatches) if the read is within
# this many bp of the anticipated target
# Eg. TARGET_WINDOW = 100 and anticipated target at 1152 would
# mean all reads from 1052-1252 would not be analyzed
TARGET_WINDOW = 100

# The expected integration site distance from the end of the protospacer
# Type 1 is usually 49 or 50bp
INTEGRATION_SITE_DISTANCE = 49

# The bin size to use for looking for areas of correlation
CORRELATIONS_BIN_SIZE = 100
