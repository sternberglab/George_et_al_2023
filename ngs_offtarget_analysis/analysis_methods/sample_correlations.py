import itertools
import os
import csv
from Bio.SeqUtils import GC
import numpy as np

def convert_reads_to_bins(reads, binsize):
	# convert from {location: count, ...} to {bin: count, ...}
	positions = sorted(reads.keys())
	bins = {}
	for position in positions:
		bin_start = (position // binsize)*binsize
		bin_end = bin_start + binsize
		bin_key = f"{bin_start}_{bin_end}"
		bins[bin_key] = bins.get(bin_key, 0) + reads[position]
	return bins


def correlation_scores(dict_of_sample_locations, binsize=None, genome=None, outputs_directory=None):
	samples = sorted(list(dict_of_sample_locations.keys()))

	sort_index = 1
	# If binsize, then group the reads into those bins
	if binsize and binsize > 1:
		for (sample, reads) in dict_of_sample_locations.items():
			dict_of_sample_locations[sample] = convert_reads_to_bins(reads, binsize)
		if genome:
			sort_index = 2

	all_positions = [reads.keys() for (sample, reads) in dict_of_sample_locations.items()]
	all_positions = sorted(set(list(itertools.chain(*all_positions))))
	
	csv_rows = []
	for p in all_positions:
		reads = []
		for sample in samples:
			reads.append(dict_of_sample_locations[sample].get(p, 0))
		nonzero_items = sum([True if v else False for v in reads])
		score = pow(sum(reads), nonzero_items)
		if binsize and binsize > 1 and genome:
			start = int(p.split('_')[0])
			end = int(p.split('_')[1])
			at_pct = 100 - GC(genome.seq[start:end])
			csv_rows.append([p, at_pct, score] + reads)
		else:
			csv_rows.append([p, score] + reads)
	csv_rows = sorted(csv_rows, key=lambda x: x[sort_index], reverse=True)


	with open(os.path.join(outputs_directory, 'correlation_areas.csv'), 'w') as outf:
		writer = csv.writer(outf)
		if binsize and binsize > 1 and genome:
			writer.writerow(['pos', 'AT%', 'score'] + samples)
		else:
			writer.writerow(['pos', 'score'] + samples)
		for row in csv_rows:
			writer.writerow(row)
	return

def pearson_coefficient(dict_of_sample_locations, binsize=None, genome=None, outputs_directory=None):
	samples = sorted(list(dict_of_sample_locations.keys()))

	# If binsize, then group the reads into those bins
	if binsize and binsize > 1:
		for (sample, reads) in dict_of_sample_locations.items():
			dict_of_sample_locations[sample] = convert_reads_to_bins(reads, binsize)

	combinations = itertools.combinations(samples, 2)
	for combo in combinations:
		combo_positions = [reads.keys() for (sample, reads) in dict_of_sample_locations.items() if sample in combo]
		combo_positions = sorted(set(list(itertools.chain(*combo_positions))))

		array = []
		for sample in combo:
			reads = dict_of_sample_locations[sample]
			array.append([reads.get(p, 0) for p in all_positions])
		array = np.array(array)
		coeff = np.corrcoef(array)
	return coeff


