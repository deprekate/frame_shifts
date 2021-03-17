#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import argparse
from argparse import RawTextHelpFormatter
from collections import Counter

class Translate:
	def __init__(self):
		nucs = ['T', 'C', 'A', 'G']
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		amino_acids = 'FFLLSSSSYY#+CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		self.translate = dict(zip(codons, amino_acids))
		self.amino_acids = sorted(set(amino_acids))

	def codon(self, codon):
		codon = codon.upper()
		if codon in self.translate:
			return self.translate[codon]
		else:
			return ''

	def seq(self, seq):
		aa = ''
		for i in range(0, len(seq), 3):
			aa += self.codon(seq[i:i+3])
		return aa

	def seq_rev(self, seq):
		aa = ''
		for i in range(0, len(seq), 3):
			aa += self.codon(self.rev_comp(seq[i:i+3]))
		return aa[::-1]

	def rev_comp(self, seq):
		seq_dict = {'A':'T','T':'A','G':'C','C':'G',
					'N':'N',
					'R':'Y','Y':'R','S':'S','W':'W','K':'M','M':'K',
					'B':'V','V':'B','D':'H','H':'D'}
		return "".join([seq_dict[base] for base in reversed(seq)])

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def same_frame(a,b):
	return (a)%3 == (b-2)%3

def gc_content(seq):
	g = seq.count('G')
	c = seq.count('C')
	a = seq.count('A')
	t = seq.count('T')
	return round( (g+c) / (g+c+a+t) , 3)
	
if __name__ == '__main__':
	usage = 'make_train.py [-opt1, [-opt2, ...]] infile'
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file in genbank format')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	parser.add_argument('-w', '--window', action="store", type=int, default=120,  help='The size of the window')
	parser.add_argument('-v', '--verbose', action="store_true", help=argparse.SUPPRESS)
	parser.add_argument('--ids', action="store", help=argparse.SUPPRESS)
	args = parser.parse_args()

	translate = Translate()

	dna = ''
	flag = False
	pairs = dict()

	with open(args.infile) as fp:
		for line in fp:
			if line.startswith('     CDS '):
				m = re.findall(r"\d+\.\.\d+", line)
				for pair in m:
					pairs[tuple(map(int, pair.split('..')))] = -1 if 'complement' in line else 1
			elif line.startswith('ORIGIN'):
				dna = '\n'
			elif dna:
				line = line[10:].replace(' ','')
				dna += line.upper()

	dna = dna.replace('\n', '')

	gc = gc_content(dna) 

	# get which frame is the coding frame
	coding_frame = dict()
	for (left,right),forward in pairs.items():
		if same_frame(left,right):
			for i in range(left, right, 3):
				coding_frame[ +(i - 1) * forward ] = True
				coding_frame[ +(i + 0) * forward ] = False
				coding_frame[ +(i + 1) * forward ] = False
				coding_frame[ -(i - 1) * forward ] = False
				coding_frame[ -(i + 0) * forward ] = False
				coding_frame[ -(i + 1) * forward ] = False
		else:
			raise ValueError('not same')
	
	# get the aminoacid frequency window
	half = int(args.window / 2)
	for i in range(0, len(dna)-2, 3):
		for f in [0,1,2]:
			n = (i+f)
			window = dna[max(0+f, n-half) : n+half]
			#print(n, window)
			print(n+1, coding_frame.get(n, None), gc, sep='\t', end='')
			freq = Counter(translate.seq(window))
			total = sum(freq.values())
			for aa in translate.amino_acids:
				print('\t', end='')
				print(round(freq.get(aa,0)/total, 3), end='')
			print()	
			#print(freq)
			print(-(n+1), coding_frame.get(-n, None), gc, sep='\t', end='')
			freq = Counter(translate.seq_rev(window))
			for aa in translate.amino_acids:
				print('\t', end='')
				print(round(freq.get(aa,0)/total, 3), end='')
			print()	




