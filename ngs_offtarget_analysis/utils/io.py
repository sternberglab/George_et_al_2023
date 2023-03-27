from pathlib import Path
import csv
from parameters import info_file, genbank_file, pipeline_outputs_directory
from Bio import SeqIO

genomes = {}

def read_info_file(filepath):
	try:
		with open(Path(info_file), 'r', encoding='utf-8-sig') as opened_info_file:
			reader = csv.DictReader(opened_info_file)
			all_rows = [row for row in reader]
	except:
		with open(Path(info_file), 'r', encoding='utf-8-sig') as opened_info_file:
			reader = csv.DictReader(opened_info_file)
			all_rows = [row for row in reader]
	sample_info_map = {row["Sample"]: row for row in all_rows}
	return sample_info_map

def read_ngs_output(sample):
	globname = f"*{sample}_target_read_locations.csv"
	files = Path(pipeline_outputs_directory).glob(globname)
	files = [f for f in files]
	if len(files) > 1:
		raise Exception(f"Found {len(files)} {globname} files in {str(Path(pipeline_outputs_directory).absolute())}. Remove the outputs you do not want to use from the folder. ")
	if len(files) == 0:
		raise Exception(f"Found {len(files)} {globname} files in {str(Path(pipeline_outputs_directory).absolute())}.")
	
	with open(files[0], 'r', encoding='utf-8-sig') as opened_reads_file:
		reader = csv.DictReader(opened_reads_file)
		reads = [r for r in reader]
	
	return reads

def get_genome(sample_info):
	# Use genbank file if it's provided
	if genbank_file and len(genbank_file):
		if not genbank_file.endswith('.gb'):
			raise Exception(f"Provided genbank file {genbank_file} should end in .gb")
		if not Path(genbank_file).exists():
			raise Exception(f"Could not find genbank file at {str(Path(genbank_file).absolute())}")
		if genbank_file not in genomes:
			genomes[genbank_file] = SeqIO.read(genbank_file, 'gb')
		return genomes[genbank_file]

	# If no genbank file, use input target fasta file
	fasta_path = sample_info['Target fasta file']
	if not Path(fasta_path).exists():
		raise Exception(f"Could not find the input target fasta at {Path(fasta_path).absolute()}")
	if fasta_path not in genomes:
		genomes[fasta_path] = SeqIO.read(fasta_path, 'fasta')
	return genomes[fasta_path]

