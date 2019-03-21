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
		
	parser.add_argument('-r', '--reference-filename',
		action="store", dest="reference_filename",
		help="Reference file that was used for alignment.", default=None)
		
	parser.add_argument('-s', '--start',
		action="store", dest="start_index",
		help="Start index for analysis.", default=None)
		
	parser.add_argument('-e', '--end',
		action="store", dest="end_index",
		help="End index for analysis.", default=None)
		
	parser.add_argument('-o', '--output-file',
		action="store", dest="output_filename",
		help="TSV output for data.", default=None)

	options = parser.parse_args()
	
	contig_name = None
	with open(options.reference_filename,'r') as f:
		contig_name = f.readline().rstrip()[1:]
		print contig_name
	
	start_index = int(options.start_index)
	end_index = int(options.end_index)
	
	all_dict = {}
	hdr_dict = {}
	for fn in options.input_filenames:
		print fn
		#CONTROL DF
		delsum_dict = {}
		mismatch_dict = {}

		samfile = pysam.AlignmentFile(fn,'rb')
		base_dict = {}
		for pileupcolumn in samfile.pileup(contig_name,start_index,end_index,max_depth=200000):
			if (pileupcolumn.pos < start_index) or (pileupcolumn.pos > end_index):
				continue

			basecount = {'A':0,'C':0,'G':0,'T':0,'Del':0}
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del:
					basecount[pileupread.alignment.query_sequence[pileupread.query_position]] += 1
				elif pileupread.is_del:
					basecount['Del'] +=1

			base_dict[pileupcolumn.pos] = basecount

		gene_df = pd.DataFrame.from_dict(base_dict).T
		gene_dels = gene_df['Del']/np.sum(gene_df,axis=1)
		gene_df = gene_df[['A','C','G','T']].divide(np.sum(gene_df[['A','C','G','T']],axis=1),axis=0)

		fn_hdr = {}
		for pos in HDR_DICT:
			fn_hdr[pos] = gene_df.loc[pos][HDR_DICT[pos]]
			
		hdr_dict[fn] = fn_hdr
		
		fig,ax = plt.subplots(figsize=(12,5))
		gene_dels.plot(kind='bar',alpha=0.25,width=0.9,ax=ax,color='red')
		all_dict[fn] = gene_dels
		#np.sum(np.abs(gene_df/2),axis=1).plot(kind='bar',alpha=0.25,width=0.9,ax=ax,color='blue')
		"""norm_df['A'].plot(kind='bar',alpha=0.25,width=0.8,ax=ax,color='blue')
		norm_df['C'].plot(kind='bar',alpha=0.25,width=0.8,ax=ax,color='green')
		norm_df['G'].plot(kind='bar',alpha=0.25,width=0.8,ax=ax,color='yellow')
		norm_df['T'].plot(kind='bar',alpha=0.25,width=0.8,ax=ax,color='orange')"""
		ax.set_ylim([0,0.4])
		xticks = np.arange(0,end_index-start_index,10)
		ax.set_xticks(xticks)
		ax.set_xticklabels(start_index+xticks)

		plt.title(fn)
		#plt.show()
		plt.savefig(fn+'.png')
		plt.close()
		del fig
		del ax
		del samfile
		
		
	df = pd.DataFrame.from_dict(all_dict).T
	df.to_csv(options.output_filename,sep='\t')
	
	hdr_df = pd.DataFrame.from_dict(hdr_dict).T
	hdr_df.to_csv(options.output_filename+'.hdr.txt',sep='\t')