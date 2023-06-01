#!/Users/jerrin.t.george/opt/miniconda3/bin/python

import os
import time
import csv
from pathlib import Path


info_file = "./Output_NGS_pipeline.csv" # add output file path here, edit output excel sheet to add Description
pipeline_input_directory = './Inputs'
target_type = 'genome'
#target_type = 'plasmid'
spacer_strand = 'reverse'
#spacer_strand = 'forward'
output_subdirectory = 'nontarget_insertion_results'
spacer_end = 334961 # check everytime

#spacer_end = 1013 # check everytimeOutput
samples = []


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

#Exclude on-target in Genome of ptarget taken
def exclude_form_target(files, start_position, end_position):
	with open(files[0], 'r', encoding='utf-8-sig') as opened_reads_file:
		reader = csv.DictReader(opened_reads_file)
		exclude_mask = [(int(row['position']) >= start_position) & (int(row['position']) <= end_position) for i, row in enumerate(reader)]
		opened_reads_file.seek(0)
		next(reader) # skip header row
		reads_excluding_range = sum(int(row['reads']) for i, row in enumerate(reader) if not exclude_mask[i])

	return reads_excluding_range

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

#Sum of reads in pTarget or Genome
def sum_total(files):
	
	with open(files[0], 'r', encoding='utf-8-sig') as opened_reads_file:
		reader = csv.DictReader(opened_reads_file)
		total_reads = sum([int(row['reads']) for row in reader])
	return total_reads

#Exclude contaminating reads in pDonor taken as the second plasmid here  

def exclude_form_second_target(files):
	with open(files[0], 'r', encoding='utf-8-sig') as opened_reads_file:
		reader = csv.DictReader(opened_reads_file)
		exclude_positions = [1198, 2333] # replace with the positions to exclude
		exclude_mask = []
		reads = []
		for row in reader:
			try:
				reads.append(int(row['reads']))
				if row['position'] != '':
					exclude_mask.append(int(row['position']) in exclude_positions)
			except ValueError:
				 # Skip rows with non-numeric 'reads'
				pass
		reads_excluding_positions = sum(reads[i] for i in range(len(reads)) if not exclude_mask[i])
	return reads_excluding_positions

def main():
	start_time = time.perf_counter()

	if target_type not in ['genome', 'plasmid']:
		raise Exception("The 'analysis_type' parameters must be 'genome', 'plasmid")

	print("Starting the run...")
	outputs_directory = os.path.join(Path(__file__).parent.absolute(), 'outputs', output_subdirectory)
	os.makedirs(outputs_directory, exist_ok=True)

	all_samples = read_info_file(info_file)
	samples_to_run = samples if samples else list(all_samples.keys())

	results = {}

	max_tn_end_reads = get_max_tn_reads(samples_to_run)

	for sample in samples_to_run:
		print(f"Analyzing {sample}...")

		sample_results = {'Sample': sample}
		sample_info = all_samples[sample]
		tn_end_reads = int(sample_info['Transposon end-containing reads'])


		if target_type in ['genome','plasmid']:
			if spacer_strand == 'forward':
				start_position = spacer_end
				end_position = spacer_end + 100
			else:
				start_position = spacer_end - 100
				end_position = spacer_end

		if target_type == 'genome':
			globname = f"*{sample}_target_read_locations.csv"
			files = Path(pipeline_input_directory).glob(globname)
			files = [f for f in files]
			if len(files) != 1:
				raise Exception(f"Found {len(files)} {globname} files in {str(Path(pipeline_input_directory).absolute())}. Remove the outputs you do not want to use from the folder. ")

			total_non_target = exclude_form_target(files, start_position, end_position)
			total_reads = sum_total(files)
			on_target = total_reads - total_non_target

			Normalized_on_target = 	float((on_target/tn_end_reads) * max_tn_end_reads)
			Normalized_total_non_target = float((total_non_target/tn_end_reads) * max_tn_end_reads)
			Normalized_total_reads = Normalized_on_target+Normalized_total_non_target

		elif target_type == 'plasmid':
			globname = f"*{sample}_target_read_locations.csv"
			files = Path(pipeline_input_directory).glob(globname)
			files = [f for f in files]
			if len(files) != 1:
				raise Exception(f"Found {len(files)} {globname} files in {str(Path(pipeline_input_directory).absolute())}. Remove the outputs you do not want to use from the folder. ")

			globname2 = f"*{sample}_second_target_read_locations.csv"
			files2 = Path(pipeline_input_directory).glob(globname2)
			files2 = [f for f in files2]
			if len(files2) != 1:
				raise Exception(f"Found {len(files2)} {globname2} files in {str(Path(pipeline_input_directory).absolute())}. Remove the outputs you do not want to use from the folder. ")

			non_target_plasmid = exclude_form_target(files, start_position, end_position)
			total_plasmid = sum_total(files)
			on_target = total_plasmid - non_target_plasmid
			total_second_plasmid = exclude_form_second_target(files2)
			total_non_target = non_target_plasmid + total_second_plasmid
			total_reads = total_plasmid + total_second_plasmid
			Normalized_on_target = 	float((on_target/tn_end_reads) * max_tn_end_reads)
			Normalized_total_non_target = float((total_non_target/tn_end_reads) * max_tn_end_reads)
			Normalized_total_reads = Normalized_on_target+Normalized_total_non_target


		else:
			print("Incorrect entry!")

		Description = sample_info['Description']
		sample_results['Description'] = Description
		sample_results['total_reads'] = total_reads
		if target_type == 'plasmid':
			sample_results['pTarget_non_target_reads'] = non_target_plasmid
			sample_results['pDonor_non_target_reads'] = total_second_plasmid


		sample_results['total_non_target_reads'] = total_non_target	
		sample_results['on_target_reads'] = on_target
		sample_results['Specificity'] = round((on_target/total_reads)*100, 1)
		sample_results['Transposon-end containing reads'] = tn_end_reads
		sample_results['Normalized total_non_target_reads'] = round(Normalized_total_non_target, 2)
		sample_results['Normalized on_target reads'] = round(Normalized_on_target, 2)
		sample_results['Normalized total reads'] = round(Normalized_total_reads, 2)
		sample_results['Fraction_total_non-target'] = round(((Normalized_total_non_target/Normalized_total_reads)*100), 2)
		sample_results['Fraction_on-target'] = round(((Normalized_on_target/Normalized_total_reads)*100), 2)


		
		print(f"Finished with sample {sample}: {time.perf_counter()-start_time} seconds")

		# Persist the sample-level results
		results[sample] = sample_results

	# Output the sample-level results
	with open(os.path.join(outputs_directory, 'analysis_summary.csv'), 'w') as outfile:
		fieldnames = results[samples_to_run[0]].keys()
		writer = csv.DictWriter(outfile, fieldnames=fieldnames)
		writer.writeheader()
		for (sample, summary_results) in results.items():
			writer.writerow(summary_results)

	print(f"Finished in {round(time.perf_counter()-start_time,2)} seconds")
main()