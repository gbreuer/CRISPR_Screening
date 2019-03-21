#trim_fastq_file.py
import os
import re
import glob
import optparse

def make_rc(sequence):
	sequence = sequence[::-1]
	sequence = sequence.replace('G','c')
	sequence = sequence.replace('A','t')
	sequence = sequence.replace('C','g')
	sequence = sequence.replace('T','a')
	
	return sequence.upper()
	
if __name__ == '__main__':
	#Options for running analysis
	parser = optparse.OptionParser()

	parser.add_option('-i', '--input-file',
		action="store", dest="input_filename",
		help="FASTQ file to be trimmed.", default=None)
		
	parser.add_option('-s', '--end_seq',
		action="store", dest="end_seq",
		help="Regular expression for matching beginning of read.", default=None)
		
	parser.add_option('-c', '--revcomp',
		action="store_true", dest="rev_comp",
		help="Flag for reverse complement checking.", default=False)
		
	parser.add_option('-o', '--output-filename',
		action="store", dest="output_filename",
		help="Output file for ranking.", default=None)

	options, args = parser.parse_args()

	input_filename = os.path.abspath(options.input_filename)
	output_filename = os.path.abspath(options.output_filename)
	if options.rev_comp:
		read_re = re.compile('('+options.end_seq+'[ACTG]+[\r\n]+)|([ACTG]+'+make_rc(options.end_seq)+')')
	else:
		read_re = re.compile(options.end_seq)
	
	seqs_read = 0
	seqs_trimmed = 0
	with open(input_filename,'r') as input_file:
		with open(output_filename,'w') as output_file:
			lines = []
			for l in input_file:
				if len(lines) < 4:
					lines.append(l)
					
				if len(lines) == 4:
					name = lines[0]
					seq = lines[1]
					strand = lines[2]
					qual = lines[3]
					
					seqs_read += 1
					re_match  = read_re.search(seq)
					if re_match:
						start_pos = re_match.start()
						end_pos = re_match.end()
						seqs_trimmed += 1
						
						seq = seq[start_pos:end_pos]
						if seq[-1] != '\n':
							seq += '\n'
						qual = qual[start_pos:end_pos]
						if qual[-1] != '\n':
							qual += '\n'
							
					for line in [name,seq,strand,qual]:
						output_file.write(line)
						
					lines = []
				else:
					continue
					
	print "Changed %s of %s seqs in file."%(seqs_trimmed,seqs_read)