#!/Users/jerrin.t.george/opt/miniconda3/bin/python

import os
import time
import csv
from pathlib import Path


info_file = "./Output_NGS_pipeline.csv" # add path here
pipeline_input_directory = './Inputs'
#Use the output excel sheet from the Illumina Pipeline
samples = []
output_subdirectory = 'normalized output'



def read_info_file(filepath):
	try:
		with open(Path(info_file), 'r', encoding='utf-8-sig') as opened_info_file:
			reader = csv.DictReader(opened_info_file)
			all_rows = [row for row in reader]
	except:
		with open(Path(info_file), 'r', encoding='ISO-8859-1') as opened_info_file:
			reader = csv.DictReader(opened_info_file)
			all_rows = [row for row in reader]
	sample_info_map = {row["Sample"]: row for row in all_rows}
	return sample_info_map

def get_max_tn_reads(samples_to_run):
	all_samples = read_info_file(info_file)
	max_tn_end_reads = 0
	for sample in samples_to_run:
		sample_info = all_samples[sample]
		print(f"Transposon-end reads for {sample} is...")
		tn_end_reads = int(sample_info['Transposon end-containing reads'])
		print (tn_end_reads)
		max_tn_end_reads = max(max_tn_end_reads, tn_end_reads)
	return max_tn_end_reads

def read_ngs_output(sample):
	globname = f"*{sample}_target_read_locations.csv"
	files = Path(pipeline_input_directory).glob(globname)
	files = [f for f in files]
	if len(files) != 1:
		raise Exception(f"Found {len(files)} {globname} files in {str(Path(pipeline_input_directory).absolute())}")
	try:
		with open(files[0], 'r', encoding='utf-8-sig') as opened_reads_file:
			reader = csv.DictReader(opened_reads_file)
			reads = [r for r in reader]
	except:
		with open(files[0], 'r', encoding='ISO-8859-1') as opened_reads_file:
			reader = csv.DictReader(opened_reads_file)
			reads = [r for r in reader]
			#reads will have a list of dictionary for each given sample
	return reads 

def main():
	start_time = time.perf_counter()

	print("Starting the run...")
	outputs_directory = os.path.join(Path(__file__).parent.absolute(), 'outputs', output_subdirectory)
	os.makedirs(outputs_directory, exist_ok=True)

	all_samples = read_info_file(info_file)
	samples_to_run = samples if samples else list(all_samples.keys())
	results = {}

	max_tn_end_reads = get_max_tn_reads(samples_to_run)

	for sample in samples_to_run:
		sample_info = all_samples[sample]
		print(f"Processing {sample}...")
		tn_end_reads = int(sample_info['Transposon end-containing reads'])
		normalized_reads_list = []
		normalized_fwd_reads_list = []
		normalized_rev_reads_list = []

		reads = read_ngs_output(sample)


		for row in reads: #iterating through each row (i.e., key valeu pair) in the 'reads' list of dictionary

			actual_read = int(row['reads'])
			normalized_read = float((actual_read / tn_end_reads) * max_tn_end_reads)
			row['Normalized Reads'] = round(normalized_read, 2)
			normalized_reads_list.append(row)

			actual_fwd_read = int(float(row['fwd strand']))
			normalized_fwd_read = float((actual_fwd_read / tn_end_reads) * max_tn_end_reads)
			row['Normalized Fwd Reads'] = round(normalized_fwd_read, 2)
			normalized_fwd_reads_list.append(row)
			

			actual_rev_read = int(float(row['rev strand']))
			normalized_rev_read = float((actual_rev_read / tn_end_reads) * max_tn_end_reads)
			row['Normalized Rev Reads'] = round(normalized_rev_read, 2)
			normalized_rev_reads_list.append(row)
		

		sample_output_path = Path(os.path.join(outputs_directory, f"{sample}_normalized_data_all.csv"))
		os.makedirs(sample_output_path.parent, exist_ok=True)
		with open(sample_output_path, 'w', newline='') as sample_out:
			fieldnames = ['position', 'reads', 'tRL', 'tLR', 'fwd strand', 'rev strand',  'Normalized Reads', 'Normalized Fwd Reads', 'Normalized Rev Reads']
			writer = csv.DictWriter(sample_out, fieldnames=fieldnames)
			writer.writeheader()
			for row in normalized_reads_list:
				writer.writerow(row)


	print(f"Finished in {round(time.perf_counter()-start_time,2)} seconds")

main()
