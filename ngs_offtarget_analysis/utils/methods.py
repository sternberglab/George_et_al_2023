from utils.constants import TSD

def get_target_site(sample_info, genome, target_distance=50):
	seq = genome.seq
	# Requires the spacer to exactly match the protospacer
	spacer = sample_info['Spacer'].upper()
	spacer_refseq = seq.upper().find(spacer)+1
	if spacer_refseq > 0:
		target_site = spacer_refseq + len(spacer) + target_distance
	else: 
		spacer_refseq = len(seq) - (seq.reverse_complement().upper().find(spacer) + len(spacer))
		target_site = spacer_refseq - target_distance + TSD
	return target_site

def add_cds_to_reads(reads, annotated_genome, essential_genes):
	# Check if the genome is annotated
	if not annotated_genome.features or len(annotated_genome.features) == 0:
		return reads

	# Annotate genbank features
	for row in reads:
		essential = False
		row_features = ""
		position = int(row['position'])
		for feature in annotated_genome.features:
			if feature.type != 'CDS':
				continue
			if position >= feature.location.start and position <= feature.location.end:
				row_features += f"{feature.qualifiers['label']}"
				try:
					full_label = feature.qualifiers['label'][0].split()
					if full_label[1] == 'CDS':
						gene_id = full_label[0]
						if gene_id and gene_id in essential_genes:
							essential = True
				except:
					continue
		if not len(row_features):
			row_features = 'intergenic / no CDS'
		row['CDS'] = row_features
		row['Essential'] = essential
	return reads