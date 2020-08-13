from __future__ import division
from collections import deque
from decimal import Decimal
from math import log
import sys
import itertools
import copy

from codon import Codon
from codon_probability import CodonProbability

class CodingEntropy(dict):
	def __init__(self, nucleotides,  window = 120):
		self.window = window//3
		self.frames = itertools.cycle([1, 2, 3])
		self.frame = None

		self.codon = deque(['-','-','-'])

		self.codons = [None] * 4
		self.codons[1] = deque([tuple(['-','-','-'])] * self.window)
		self.codons[2] = deque([tuple(['-','-','-'])] * self.window)
		self.codons[3] = deque([tuple(['-','-','-'])] * self.window)		

		self.frequency = [None] * 4
		self.frequency[1] = self._init_dict()
		self.frequency[2] = self._init_dict()
		self.frequency[3] = self._init_dict()
		
		self[0] = deque([])
		self[1] = deque([])
		self[2] = deque([])
		self[3] = deque([])
		self[-1] = deque([])
		self[-2] = deque([])
		self[-3] = deque([])

		nucs = ['T', 'C', 'A', 'G']
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		self.amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		self.translate = dict(zip(codons, self.amino_acids))

		self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-'}
		
		self.codon_probs = CodonProbability(nucleotides)

		for i, nucl in enumerate(nucleotides,1):
			# find nucleotide frequency
			self.frame = next(self.frames)
			self.add_base(nucl)
			position = (i-1) // 3

			self[0].append(self.frequency[self.frame].copy())

			# find shannon entropy
			#se = self.trinucleotide_entropy(self.frequency[self.frame])
			#se = self.dinucleotide_entropy(self.frequency[self.frame])
			#print(i, nucl, self.codon)	
			se = self.peptide_entropy(self.frequency[self.frame])
			self[self.frame].append(se)
			#try:
			#	self[position].append(se)
			#except:
			#	self.append([])
			#	self[position].append(se)
				
			
			#self[self.frame].append(se)

			#se = self.dinucleotide_entropy(self.reverse_frequencies(self.frequency[self.frame]))
			se = self.peptide_entropy(self.reverse_frequencies(self.frequency[self.frame]))
			self[-self.frame].append(se)
			#try:
			#	self[position].append(se)
			#except:
			#	self.append(list)
			#	self[position].append(se)
			#self[-self.frame].append(se)
		self.end()
			

	def end(self):
		#remove half of first
		self[0].popleft()
		for _ in range(self.window//2):
			for frame in [1,2,3]:
				self[frame].popleft()
				self[-frame].popleft()
				self[0].popleft()
		#add half of last
		for _ in range(self.window//2):
			for frame in [1,2,3]:
 				# pop oldest codon
				popped_codon = self.codons[frame].popleft()
				self.frequency[frame][popped_codon] -= 1
				se = self.peptide_entropy(self.frequency[frame])
				self[frame].append(se)
				se = self.peptide_entropy(self.reverse_frequencies(self.frequency[frame]))
				self[-frame].append(se)
			self[0].append(self.frequency[frame].copy())

	def reverse_frequencies(self, dictionary):
		new_dict = dict()
		for key in dictionary:
			new_key = tuple([self.complement[t] for t in key[::-1]])
			new_dict[new_key] = dictionary[key]
		return new_dict

	def translate_dict(self, dictionary):
		aminoacid_dictionary = dict()
		for codon in dictionary:
			if '-' not in codon and dictionary[codon]:
				aa = self.translate[''.join(codon)]
				count = aminoacid_dictionary.get(aa, 0)
				aminoacid_dictionary[aa] = dictionary[codon] + count 
		return aminoacid_dictionary

	def add_base(self, base):

		# set newest codon
		self.codon.popleft()
		self.codon.append(base)

 		# add newest codon
		self.codons[self.frame].append(tuple(self.codon))
		self.frequency[self.frame][tuple(self.codon)] += 1
		# misc
		#aa = self.translate.get("".join(self.codon), "-")
		#self.aa[self.frame][aa] += 1

 		# pop oldest codon
		popped_codon = self.codons[self.frame].popleft()
		self.frequency[self.frame][popped_codon] -= 1
		# misc
		#aa = self.translate.get("".join(popped_codon), "-")
		#self.aa[self.frame][aa] -= 1

	def peptide_entropy(self, dictionary):
		se = 0;
		new_dict = dict()
		total = 0
		for key in dictionary:
			if '-' not in key and dictionary[key]:
				aa = self.translate[''.join(key)]
				count = new_dict.get(aa, 0)
				new_dict[aa] = dictionary[key] + count 
		for key in new_dict:
			new_dict[key] = new_dict[key] * self.codon_probs.probability(aa) #/ self.amino_acids.count(key)
			#new_dict[key] = new_dict[key] / self.amino_acids.count(key)
			total += new_dict[key]
		for key in new_dict:
			p = new_dict[key] / total
			#p = new_dict[key] / 22
			se += -p * log(p)
		return se

	def dinucleotide_entropy(self, dictionary):
		se = 0;
		new_dict = dict()
		for key in dictionary:
			if '-' not in key and dictionary[key]:
				count = new_dict.get(key[0:2] , 0)
				new_dict[key[0:2]] = count + 1
		for key in new_dict:
			p = new_dict[key] / 16
			se += -p * log(p)
		return se

	def trinucleotide_entropy(self, dictionary, direction = 'forward'):
		se = 0;
		for key in dictionary:
			if '-' not in key and dictionary[key]:
				p = dictionary[key] / 64
				se += -p * log(p)
		return se

	def _init_dict(self):
		freq_dict = dict()
		for first in 'ATGC-':
			for second in 'ATGC-':
				for third in 'ATGC-':
					freq_dict[(tuple([first, second, third]))] = 0
		return freq_dict
				
	def _init_aadict(self):
		freq_dict = dict()
		amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-'
		for aa in amino_acids:
			freq_dict[aa] = 0
		return freq_dict

