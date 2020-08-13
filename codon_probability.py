
class CodonProbability:

	def __init__(self, dna):
		base_counts = {'A':0, 'T':0, 'G':0, 'C':0}
		for base in dna:
			base_counts[base] += 1
		A = base_counts['A'] / sum(base_counts.values())
		T = base_counts['T'] / sum(base_counts.values())
		G = base_counts['G'] / sum(base_counts.values())
		C = base_counts['C'] / sum(base_counts.values())
		self.codons = {
			'F' : T*T*T + T*T*C,
			'L' : C*T*T + C*T*C + C*T*A + C*T*G + T*T*A + T*T*G,
			'I' : A*T*T + A*T*C,
			'M' : A*T*G,
			'V' : G*T*T + G*T*C + G*T*A + G*T*G,
			'S' : T*C*T + T*C*C + T*C*A + T*C*G + A*G*T + A*G*C,
			'P' : C*C*T + C*C*C + C*C*A + C*C*G,
			'T' : A*C*T + A*C*C + A*C*A + A*C*G,
			'A' : G*C*T + G*C*C + G*C*A + G*C*G,
			'Y' : T*A*T + T*A*C,
			'*' : T*A*A + T*A*G + T*G*A,
			'H' : C*A*T + C*A*C,
			'Q' : C*A*A + C*A*G,
			'N' : A*A*T + A*A*C,
			'K' : A*A*A + A*A*G,
			'D' : G*A*T + G*A*C,
			'E' : G*A*A + G*A*G,
			'C' : T*G*T + T*G*C,
			'W' : T*G*G,
			'R' : C*G*T + C*G*C + C*G*A + C*G*G + A*G*A + A*G*G,
			'G' : G*G*T + G*G*C + G*G*A + G*G*G
		}

	def probability(self, codon):
		return self.codons[codon]

#my_probs = CodonProbability('ATGC')
#print(my_probs.probability('M'))
#print(my_probs.probability('S'))
		

