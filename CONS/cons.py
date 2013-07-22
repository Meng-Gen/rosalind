import sys

def get_dna_strings(lines):
	dna_strings = {}
	dna_length = 0
	is_first_dna = True
	id = ''
	dna = ''
	for line in lines:
		line = line.strip();
		if line[0] == '>':
			if is_first_dna:
				is_first_dna = False
			else:
				dna_strings.update({id : dna})
			id = line[1:]
			dna = ''
		else:
			dna += line
		dna_strings.update({id : dna})
		dna_length = len(dna)
	return dna_strings, dna_length

def get_freqs(list):
	freqs = {}
	for c in list:
		freqs[c] = freqs.get(c, 0) + 1
	max_freq = 0
	max_char = ''
	for c in freqs:
		if freqs[c] > max_freq:
			max_freq = freqs[c]
			max_char = c
	return freqs, max_char
	
def main():
	dna_strings, dna_length = get_dna_strings(sys.stdin.readlines())
	consensus = []
	profile_matrix = []
	for i in range(dna_length):
		freqs, max_char = get_freqs([dna_strings[id][i] for id in dna_strings])
		consensus.append(max_char)
		profile_matrix += [freqs]
	print(''.join(consensus))
	for c in 'ACGT':
		print('%c: ' % c, end='')
		for i in range(dna_length):
			if c in profile_matrix[i]:
				print(profile_matrix[i][c], end='')
			else:
				print(0, end='')
			if i == dna_length - 1:
				print()
			else:
				print(' ', end='')

if __name__ == '__main__':
    sys.exit(main())
