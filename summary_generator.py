import pandas as pd
import argparse
import os
import ast
import numpy as np
from glob import glob

SKIP_ZERO = True

if __name__ == '__main__':
	#Options for running analysis
	parser = argparse.ArgumentParser()

	parser.add_argument('-i', '--input-file',
		action="store", dest="input_filenames",
		help="FASTQ file to be trimmed.", nargs='+', default=[])
		
	parser.add_argument('-o', '--output-filename',
		action="store", dest="output_filename",
		help="Output file for ranking.", default=None)

	options = parser.parse_args()
	
	print options.input_filenames
	cleaned_filenames = []
	for fn in options.input_filenames:
		cleaned_filenames.append(os.path.abspath(fn))
	print cleaned_filenames
	output_filename = os.path.abspath(options.output_filename)
		
	#Deletions
	summary_dict = {}
	for fn in cleaned_filenames:
		data_dict = {}
		df = pd.read_csv(fn,header=1,sep='\t')
		deletions = df['deletion_lens']
		micro = df['deletion_is_micro']
		count = df['count']
		cleaned_dels = []
		micro_dels = []
		for i in range(len(deletions)):
			d = deletions[i]
			d_list = np.array(ast.literal_eval(d)).astype(int)
			if d_list.size == 0:
				if SKIP_ZERO:
					continue
				else:
					cleaned_dels += [0]*count[i]
					micro_dels += [False]*count[i]
			else:
				cleaned_dels += [np.max(d_list)]*count[i]
				
				if 'True' in micro[i]:
					micro_dels += [True]*count[i]
				else:
					micro_dels += [False]*count[i]
				
		data_dict['Deletion_Size'] = cleaned_dels
		data_dict['Mircohomology'] = micro_dels
		
		df = pd.DataFrame.from_dict(data_dict)
		df.to_csv(fn+'_dels.txt',sep='\t')
		summary_dict[fn] = cleaned_dels

	#Insertions
	for fn in cleaned_filenames:
		data_dict = {}
		df = pd.read_csv(fn,header=1,sep='\t')
		insertions = df['insertion_seqs']
		count = df['count']
		cleaned_ins = []
		for i in range(len(insertions)):
			d = insertions[i]
			d_list = np.array([len(ins) for ins in ast.literal_eval(d)]).astype(int)
			if d_list.size == 0:
				if SKIP_ZERO:
					continue
				else:
					cleaned_ins += [0]*count[i]
			else:
				cleaned_ins += [np.max(d_list)]*count[i]
				
		data_dict[fn] = cleaned_ins
		
		df = pd.DataFrame.from_dict(data_dict)
		df.to_csv(fn+'_ins.txt',sep='\t')
		
	out_df = pd.DataFrame.from_dict(summary_dict,orient='index').T
	out_df = out_df[sorted(out_df.columns)]
	out_df.to_csv(output_filename,sep='\t')