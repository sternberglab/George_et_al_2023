import random
from Bio.SeqUtils import GC
from utils.constants import AT_PCT_TEST_SIZE

at_test_datasets = {}

def get_at_test_data(genome, sequence_length):
	test_dataset_key = f"{genome.id}__{sequence_length}"
	if test_dataset_key in at_test_datasets:
		return at_test_datasets[test_dataset_key]

	genome_length = len(genome.seq)
	genome_expanded = genome[-sequence_length:].seq + genome.seq + genome[:sequence_length].seq
	at_opts = {}
	for i in range(AT_PCT_TEST_SIZE):
		start = random.randint(0, genome_length) + sequence_length
		downstream_at_pct = 100-GC(genome_expanded[start:start+sequence_length])
		upstream_at_pct = 100-GC(genome_expanded[start-sequence_length:start])
		best_at_pct = max([downstream_at_pct, upstream_at_pct])
		at_pct = round(best_at_pct, 2)
		if at_pct not in at_opts:
			at_opts[at_pct] = 0
		at_opts[at_pct] += 1
	at_test_datasets[test_dataset_key] = at_opts
	return at_test_datasets[test_dataset_key]

def get_at_likelihood(at_pct, random_sampling_dict):
	total_gte_at = 0
	for (k,v) in random_sampling_dict.items():
		if k >= at_pct:
			total_gte_at += v
	return total_gte_at / AT_PCT_TEST_SIZE

