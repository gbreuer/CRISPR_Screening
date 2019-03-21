#draw indel picture
import pysam
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

HDR_DICT = {321:'A',333:'C',346:'A',356:'A'}

if __name__ == '__main__':
	#Options for running analysis
	parser = argparse.ArgumentParser()

	parser.add_argument('-i', '--input-file',
		action="store", dest="input_filenames",
		help="BAM files to be analyzed.", nargs='+', default=[])
		
	parser.add_argument('-o', '--output-file',
		action="store", dest="output_filename",
		help="TSV output for data.", default=None)

	options = parser.parse_args()
	
	file_dict = {}
	for fn in options.input_filenames:
		print fn
		hdr_counts = [0,0,0,0,0]
		samfile = pysam.AlignmentFile(fn,'rb')
		for read in samfile:
			hdr_count = 0
			query_locs, ref_locs = zip(*read.get_aligned_pairs())
			for location in HDR_DICT:
				try:
					query_loc = query_locs[ref_locs.index(location)]
				except:
					continue
					
				if query_loc == None:
					continue
				elif read.seq[query_loc] == HDR_DICT[location]:
					hdr_count += 1
					
			hdr_counts[hdr_count] += 1
			
		file_dict[fn] = hdr_counts
		
	final_output = pd.DataFrame.from_dict(file_dict).T
	final_output.to_csv(options.output_filename,sep='\t')