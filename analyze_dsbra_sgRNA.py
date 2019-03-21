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
import optparse

sgRNA_RE = re.compile(r'ACACCG([ACGT]{20})')
	
def count_hdr(samfilename,hdr_dict):
    samfile = pysam.AlignmentFile(samfilename,'rb')
    start_pileup = np.min(hdr_dict.keys())
    end_pileup = np.max(hdr_dict.keys())
    coverage = 0
    
    read_dict = {}
    for pileupcolumn in samfile.pileup('LCv2_no_sgRNA',start_pileup,end_pileup,max_depth=500000):
        if not pileupcolumn.pos in hdr_dict.keys():
            continue
            
        for pileupread in pileupcolumn.pileups:
            if not pileupread.alignment.query_name in read_dict:
                read_dict[pileupread.alignment.query_name] = {k:False for k in hdr_dict.keys()}
                read_dict[pileupread.alignment.query_name]['Del'] = False
                coverage = pileupcolumn.n
                
            if pileupread.is_del:
                read_dict[pileupread.alignment.query_name]['Del'] = True
            elif pileupread.alignment.query_sequence[pileupread.query_position] == hdr_dict[pileupcolumn.pos]:
                read_dict[pileupread.alignment.query_name][pileupcolumn.pos] = True

    return (pd.DataFrame.from_dict(read_dict).T,coverage)
	
START_PILEUP = 300
END_PILEUP = 400
HDR_DICT = {321:'A',333:'C',346:'A',356:'A'}

if __name__ == '__main__':
	#Options for running analysis
	parser = optparse.OptionParser()

	parser.add_option('-i', '--input-directory',
		action="store", dest="input_directory",
		help="Directory containing files for analysis.", default=None)
		
	parser.add_option('-o', '--output-filename',
		action="store", dest="output_filename",
		help="Output file for parsed files.", default=None)

	options, args = parser.parse_args()

	input_directory = os.path.abspath(options.input_directory)
	output_filename = os.path.abspath(options.output_filename)
	print input_directory
	print output_filename
	
	samplesum_dict = {}
	for filename in glob.glob(os.path.join(input_directory,'*.sam.sorted.bam')):
		try:
			reads_dict,cov = count_hdr(filename,HDR_DICT)
			sum_hdr = np.sum(reads_dict[HDR_DICT.keys()],axis=1)
			gene_name = filename[filename.rfind('/')+1:filename.rfind('.sam.sorted.bam')]
			hdrsum_dict = np.sum(reads_dict[HDR_DICT.keys()],axis=0)
			hdrsum_dict['Coverage'] = cov
			hdrsum_dict['Deletions'] = np.sum(reads_dict['Del'])
			hdrsum_dict['HDR0'] = np.sum(sum_hdr==0)
			hdrsum_dict['HDR1'] = np.sum(sum_hdr==1)
			hdrsum_dict['HDR2'] = np.sum(sum_hdr==2)
			hdrsum_dict['HDR3'] = np.sum(sum_hdr==3)
			hdrsum_dict['HDR4'] = np.sum(sum_hdr==4)

			samplesum_dict[gene_name] = hdrsum_dict
		except:
			print "skipping "+filename
	
	pd.DataFrame.from_dict(samplesum_dict).T.to_csv(output_filename,sep='\t')