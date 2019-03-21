import os
import re
import pandas as pd
import numpy as np
import optparse
import gzip

#Constants -- CHANGE THESE
sgRNA_RE = re.compile(r'GTGGAAAGGACGAAACACCG([ACGT]{20})')
sgRNA_RC_RE = re.compile(r'GCTATTTCTAGCTCTAAAAC([ACGT]{20})')
FASTQ_RE = re.compile(r'\@M037.*\s+.*\s+.*\s+\+\s+.*\s+')
FASTQ_READNAME_RE = re.compile(r'\@M037[\w:\-]+')

def make_rc(sequence):
	sequence = sequence[::-1]
	sequence = sequence.replace('G','c')
	sequence = sequence.replace('A','t')
	sequence = sequence.replace('C','g')
	sequence = sequence.replace('T','a')
	
	return sequence.upper()

#Constants
REMOVE_DUPLICATES = False
SAVE_ALL = False
COLLAPSE_TO_GENE = False

if __name__ == "__main__":
	print os.path.dirname(os.path.abspath(__file__))
	dsbra_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'dsbra')

	#Options for running analysis
	parser = optparse.OptionParser()

	parser.add_option('-i', '--input-filename',
		action="store", dest="input_filename",
		help="FASTQ file for analysis.", default=None)
		
	parser.add_option('-o', '--output-folder',
		action="store", dest="output_folder",
		help="Output folder for parsed files.", default=None)
		
	parser.add_option('-s', '--sgrna-list',
		action="store", dest="sgrna_list",
		help="Path to sgRNA list for parsing.", default=None)
		
	parser.add_option('-t', '--threads',
		action="store",dest="num_threads",
		help="Number of threads to use when writing to disk.", default=1)
		
	parser.add_option('-c', '--compress',
		action="store_true",dest="gene_only",
		help="Combine sgRNAs to gene level only.", default=False)
		
	options, args = parser.parse_args()

	print 'Input file:', options.input_filename
	print 'Output folder:', os.path.abspath(options.output_folder)
	print 'Number of threads:', options.num_threads

	sample_filename = options.input_filename
	#REF_SEQ = options.reference_filename
	output_dir = options.output_folder
	COMPRESS_TO_GENE = options.gene_only
	
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	CRISPR_db = pd.read_csv(options.sgrna_list,header=1)
	CRISPR_db.index = CRISPR_db['Sequence Name']
	CRISPR_db["Sequence (5'-3')"]
	sgRNA_Names = CRISPR_db.index
	sgRNA_seqs = np.array(CRISPR_db["Sequence (5'-3')"])
	sgRNA_seqs = [s[20:40] for s in sgRNA_seqs]
	gene_dict = {sgRNA_seqs[i]:sgRNA_Names[i] for i in range(len(sgRNA_seqs))}

	#Read fastq file to list
	#TODO: optimize this so you don't have to store in memory
	fastq_seqs = []

	if sample_filename:
		suffix = sample_filename[sample_filename.rfind('.')+1:]
		if suffix == 'gz':
			f = gzip.open(sample_filename,'rb')
		else:		
			f = open(sample_filename,'r')
			
		temp_lines = f.readlines(1000)
		while len(temp_lines)>0:
			fastq_seqs += FASTQ_RE.findall(''.join(temp_lines))
			temp_lines = f.readlines(1000)
		f.close()
		print "Done Reading"

		read_indices = {}
		for i in range(len(fastq_seqs)):
			read = fastq_seqs[i]
			sgSeq = None

			re_search = sgRNA_RE.search(read)
			if re_search is None:
				re_search = sgRNA_RC_RE.search(read)
				if re_search is not None:
					sgSeq = make_rc(re_search.groups(0)[0])
			else:
				sgSeq = re_search.groups(0)[0]        
			
			if sgSeq is not None:
				if sgSeq in read_indices:
					read_indices[sgSeq].append(i)
				else:
					read_indices[sgSeq] = [i]
					
		for sgRNA in read_indices.keys():
			if sgRNA in gene_dict:
				sgRNA_name = gene_dict[sgRNA]
								
				if COMPRESS_TO_GENE:
					sgRNA_name = sgRNA_name[:sgRNA_name.find('_')].rstrip()
				
				fname = os.path.join(output_dir,sgRNA_name+'.fq')
				with open(fname,'a') as f:
					for index in read_indices[sgRNA]:
						f.write(fastq_seqs[index])

		del fastq_seqs