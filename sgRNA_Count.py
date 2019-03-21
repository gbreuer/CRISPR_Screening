import os
import re
import pandas as pd
import numpy as np
import sys
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import platform
import pysam
import glob
import argparse
import gzip

sgRNA_RE = re.compile(r'GTGGAAAGGACGAAACACCG([ACGT]{20})')

if __name__ == '__main__':
	#Options for running analysis
	parser = argparse.ArgumentParser()

	parser.add_argument('-i', '--input-file',
		action="store", dest="input_filenames",
		help="Fastq files to be analyzed.", nargs='+', default=[])
		
	parser.add_argument('-o', '--output-filename',
		action="store", dest="output_filename",
		help="Output file for ranking.", default=None)
		
	parser.add_argument('-s', '--sgrna-list',
		action="store", dest="sgrna_list",
		help="Path to sgRNA list for parsing.", default=None)

	options = parser.parse_args()

	input_filenames = options.input_filenames
	output_filename = options.output_filename
	
	CRISPR_db = pd.read_csv(options.sgrna_list,header=1)
	CRISPR_db.index = CRISPR_db['Sequence Name']
	CRISPR_db["Sequence (5'-3')"]
	sgRNA_Names = CRISPR_db.index
	sgRNA_seqs = np.array(CRISPR_db["Sequence (5'-3')"])
	sgRNA_seqs = [s[20:40] for s in sgRNA_seqs]
	gene_dict = {sgRNA_seqs[i]:sgRNA_Names[i] for i in range(len(sgRNA_seqs))}
	
	sgRNA_Counts = {}
	
	for i in range(len(input_filenames)):	
		sample_filename = input_filenames[i]
		print sample_filename
		suffix = sample_filename[sample_filename.rfind('.')+1:]
		if suffix == 'gz':
			f = gzip.open(sample_filename,'rb')
		else:		
			f = open(sample_filename,'r')

		sgRNA_list = sgRNA_RE.findall(f.read())
		f.close()

		for sg in sgRNA_list:
			try:
				gene = gene_dict[sg]
			except:
				#SKIP ANYTHING THAT ISN'T EXPECTED
				continue
				gene = sg

			if not gene in sgRNA_Counts:
				sgRNA_Counts[gene] = [0]*len(input_filenames)
			sgRNA_Counts[gene][i] += 1

	sg_df = pd.DataFrame.from_dict(sgRNA_Counts)
	sg_df.index = input_filenames

	NORMALIZE_DATA = False

	sg_dft = sg_df.copy().transpose()
	if NORMALIZE_DATA:
		for c in sg_dft.columns:
			sg_dft[c] = sg_dft[c]/np.sum(sg_dft[c])
			
	sg_dft.to_csv(output_filename,sep='\t')