import pandas as pd
from Bio import SeqIO
from .utils import output_path
from parameters import target_site_duplication_length as TSD

def extract_sequences_from_output(reads_output_csv, meta_info, seq_output_path):
	window = 50 # window to choose upstream and downstream
	
	data = pd.read_csv(reads_output_csv)
	seq = SeqIO.read(meta_info['Target fasta file'], 'fasta').seq.upper()
	# add edges to account for the window around the site
	seq = seq[-window:] + seq + seq[:window]

	spacer = meta_info['Spacer'].upper()
	genome_length=len(seq)

	output_seqs = data.drop(['reads'], axis=1)
	output_seqs = output_seqs.drop(['tLR'], axis=1)
	output_seqs = output_seqs.drop(['tRL'], axis=1)
	output_seqs['Seq_+'] = ""
	output_seqs['Seq_-'] = ""
	ontarget = seq.find(spacer)
	if ontarget < 1:
		ontarget = genome_length - seq.reverse_complement().find(spacer) - len(spacer) #define end of spacer position to omit on-target if any

	cutoff_reads=1

	ontarget_range = list(range(ontarget-1, ontarget+100))
	for row in range(0,data.shape[0]): #data.shape[0] will be the no. of rows in the csv file, iterating through each read value
		position = data.position[row] + window
		fw_ct = data.loc[row, 'fwd strand']
		rv_ct = data.loc[row, 'rev strand']
		# Correct for TSD. 
		# The remainder of this code assumes all insertions are oriented as LR
		if meta_info['transposon_end_side'].lower() == 'left':
			fw_position = position - TSD
			rv_position = position
		else:
			fw_position = position
			rv_position = position - TSD

		if fw_position not in ontarget_range and fw_ct >= cutoff_reads:
			fw_seq = seq[fw_position - window:fw_position + window]
			# make fwd strand, left primer reads RC bc its coming in as LR and reads are coming OUT of the tn
			if meta_info['transposon_end_side'].lower() == 'left':
				fw_seq = fw_seq.reverse_complement()
			output_seqs.loc[row, 'Seq_+'] = fw_seq
		if rv_position not in ontarget_range and rv_ct >= cutoff_reads:
			rv_seq = seq[rv_position - window:rv_position + window]
			if meta_info['transposon_end_side'].lower() != 'left':
				rv_seq = rv_seq.reverse_complement()
			output_seqs.loc[row, 'Seq_-'] = rv_seq

	output_seqs.to_csv(seq_output_path, index=False)
