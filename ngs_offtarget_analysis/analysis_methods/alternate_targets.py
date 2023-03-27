from utils.constants import TSD, GENOME_BUFFER

def get_likely_target(reads_site, genome, spacer, target_distance=50):
	# this function find the most likely target sequence for a given read
	# looking at windows 35-65bp on both sides of the integration site
	# using similarity to spacer as criteria, particularly in seed region
	read_position = reads_site+GENOME_BUFFER
	SPACER_LEN = len(spacer)

	debug=False

	best_match = None
	# evaluate targets
	for i in range(-10,10,1):
		# fwd target
		fwd_tgt_site = read_position - target_distance - SPACER_LEN - i
		tgt_seq = genome[fwd_tgt_site:fwd_tgt_site + SPACER_LEN]
		tgt_position = read_position-i-SPACER_LEN-GENOME_BUFFER
		best_match = evaluate_best_match(best_match, tgt_seq, tgt_position, True, spacer, debug=debug)
		# rv target
		rv_tgt_site = read_position - TSD + target_distance + i
		tgt_seq = genome[rv_tgt_site:rv_tgt_site + SPACER_LEN].reverse_complement()
		tgt_position = read_position+i-GENOME_BUFFER
		best_match = evaluate_best_match(best_match, tgt_seq, tgt_position, False, spacer, debug=debug)
	return best_match

def evaluate_best_match(existing_best, tgt_seq, position, is_fwd_strand, spacer, without_5ths=False, debug=False):
	effective_tgt = tgt_seq
	effective_spacer = spacer
	if without_5ths:
		effective_spacer = ''.join([s for (idx, s) in enumerate(spacer) if ((idx+1) % 6 != 0)])
		effective_tgt = ''.join([s for (idx, s) in enumerate(tgt_seq) if ((idx+1) % 6 != 0)])

	matches = list(tgt_bp == spacer_bp for tgt_bp, spacer_bp in zip(effective_tgt, effective_spacer))
	group_size = 7 if without_5ths else 8

	seed_match = sum(matches[:group_size])
	total_match = sum(matches)
	distal_match = sum(matches[-group_size:])
	best_match = existing_best


	if not best_match or \
		best_match['seed_match'] < seed_match or \
		(best_match['seed_match'] == seed_match and best_match['total_match'] < total_match):
		best_match = {
			'seed_match': seed_match,
			'total_match': total_match,
			'distal_match': distal_match,
			'seq': tgt_seq,
			'pos': position,
			'strand': 'fwd' if is_fwd_strand else 'rv'
		}
	return best_match