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

def score_ranklist(rownames,ranklist,weighted=False,showplot=False,ax=None,label=None):

    miss_score = np.sum(np.abs(ranklist.loc[rownames]))/(len(ranklist)-len(rownames))
    
    max_score = 0
    wheremax = 0
    score = 0.
    min_score = 0
    wheremin = 0
    
    scoretrack = []
    
    for i in range(len(ranklist.index)):
        if ranklist.index[i] in rownames:
            score += abs(ranklist[i])
        else:
            score -= miss_score
            
        scoretrack.append(score)
            
    #figure out actual score
    min_score = min(scoretrack)
    max_score = max(scoretrack)
    if abs(min_score) > abs(max_score):
        final_score = (min_score)
    else:
        final_score = (max_score)
        
    if showplot:
        if not ax:
            fig, ax = plt.subplots()
        p = ax.plot(np.array(scoretrack)/len(rownames),label=label)
        #ax.plot(ax.get_xlim(),[0,0],linestyle='--',color='black')
        return (final_score/len(rownames),p)
        
    return final_score/len(rownames)

if __name__ == '__main__':
	#Options for running analysis
	parser = optparse.OptionParser()

	parser.add_option('-i', '--input-file',
		action="store", dest="input_filename",
		help="CSV file containing sgRNAs, associated genes, and scoring.", default=None)
		
	parser.add_option('-o', '--output-filename',
		action="store", dest="output_filename",
		help="Output file for ranking.", default=None)

	options, args = parser.parse_args()

	input_filename = os.path.abspath(options.input_filename)
	output_filename = os.path.abspath(options.output_filename)
	
	new_df = pd.read_csv(input_filename,sep='\t',index_col=0)
	new_df.apply(pd.to_numeric, errors='ignore')
	values = new_df.iloc[:,1:]
	values[~np.isfinite(values)] = 0
	new_df.iloc[:,1:] = values
	
	genes = np.unique(new_df['Gene'])
	score_dict = {}
	for g in genes:
		score_dict[g] = {}
	for measure in new_df.columns[1:]:
		print measure
		data = new_df[measure].sort_values()
		data = (data-np.mean(data))/np.std(data)
		for g in genes: 
			score_dict[g][measure] = score_ranklist(new_df.index[new_df['Gene']==g], data)
			
	score_df = pd.DataFrame.from_dict(score_dict).T
	score_df.to_csv(output_filename,sep='\t')