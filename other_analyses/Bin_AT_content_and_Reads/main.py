import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC

def main():

	bin_size = 46 # window to choose
	data = pd.read_csv('/SampleName_target_read_locations.csv')
	seq = SeqIO.read('./psl4043.fasta', 'fasta').seq.upper()
	# add edges to account for the window around the site
	Extend_seq = seq[-bin_size:] + seq + seq[:bin_size]

	genome_length=len(seq)+1 # pTarget has 2713 bases, 2714 bases were divided into 59 bins with each bin having 46 bp size

	max_bin=genome_length//bin_size + 1
	print(" Binning initiated... ")
	reads_per_bin = []


	dict_bins = {'Bin': [], 'Position': [], 'Reads': [], 'AT_content': []}
	Bin=1

	for i in range(bin_size, max_bin*bin_size, bin_size): #iterating through each bin
		temp = 0
		dict_bins['Bin'].append(Bin)
		AT_per=100-GC(Extend_seq[(i-bin_size):i])
		Position = i-(bin_size/2)
		dict_bins['Position'].append(Position)
		dict_bins['AT_content'].append(AT_per)

		for row in range(0,data.shape[0]): #data.shape[0] will be the no. of rows in the csv file, iterating through each read value
			if ((i-bin_size)<data.position[row]<=i):
				temp = temp + data.reads[row]

		dict_bins['Reads'].append(temp)
		print("Iterating through Bin ", Bin ,"/", max_bin)
		Bin +=1

	output_bin = pd.DataFrame(dict_bins, columns=['Bin', 'Position', 'Reads', 'AT_content'])
	output_bin.to_csv('Binned_Output_SampleName.csv', index = False)
	print("Done")
main()