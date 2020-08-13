import sys
import os.path
import textwrap
import gzip

import argparse
from argparse import RawTextHelpFormatter

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def get_args():
	usage = 'phanotate.py [-opt1, [-opt2, ...]] infile'
	parser = argparse.ArgumentParser(description='PHANOTATE: A phage genome annotator', formatter_class=RawTextHelpFormatter, usage=usage)

	parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')

	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	parser.add_argument('-f', '--outfmt', action="store", default="tabular", dest='outfmt', \
				help='format of the output [tabular]', choices=['tabular','genbank','fasta'])
	parser.add_argument('-w', '--windowsize', action="store", default=40, type=int, help='The size of the scrolling window [40]')

	args = parser.parse_args()

	return args


def read_fasta(filepath):
	my_contigs = dict()
	name = ''
	seq = ''
	taxa = ''
	with gzip.open(filepath , 'rt') as f:
		for line in f:
			if(line.startswith(">")):
				my_contigs[name] = seq
				name, taxa = line.split(' ', 1)
				try:
					name = taxa.split('[')[1].split('|')[0].strip()
				except:
					pass
				seq = ''
			else:
				seq += line.replace("\n", "").upper()
		my_contigs[name] = seq

	if '' in my_contigs:
		del my_contigs['']

	return my_contigs

def read_gff(filepath):
	my_frames = dict()	
	with gzip.open(filepath , 'rt') as f:
		for line in f:
			if not line.startswith(('#','\n')):
				column = line.split()
				# forward direction
				if column[6] == '+':
					beg = int(column[3])
					end = int(column[4]) - 2
					frame = ((beg - 1) % 3 ) + 1
					other = ((end - 1) % 3 ) + 1
					if frame == other:
						for i in range(beg, end):
							my_frames[i] = frame
				# reverse
				elif column[6] == '-':
					beg = int(column[3])
					end = int(column[4]) - 2
					frame = ((beg - 1) % 3 ) + 1
					other = ((end - 1) % 3 ) + 1
					if frame == other:
						for i in range(beg, end):
							my_frames[i] = -frame
				# unknown
				else:
					raise ValueError("A gene was found that is neither forward or reverse")
	return my_frames
			




