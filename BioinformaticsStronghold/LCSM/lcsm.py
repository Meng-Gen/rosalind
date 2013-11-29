import sys

def get_dna_strings(lines):
	dna_strings = []
	is_first_dna = True
	dna = ''
	for line in lines:
		line = line.strip();
		if line[0] == '>':
			if is_first_dna:
				is_first_dna = False
			else:
				dna_strings.append(dna)
			dna = ''
		else:
			dna += line
	dna_strings.append(dna)
	return dna_strings

def has_common_substring(substring, dna_strings):
	for i in range(1, len(dna_strings)):
		if dna_strings[i].find(substring) == -1:
			return False
	return True

def main():
	dna_strings = get_dna_strings(sys.stdin.readlines())
	base = dna_strings[0]
	substrings = set()
	for i in range(len(base)):
		for j in range(i, len(base)):
			substrings.add(base[i:j+1])
	max_len = 0
	max_substring = ''
	for substring in substrings:
		if len(substring) < max_len:
			continue
		if has_common_substring(substring, dna_strings):
			max_len = len(substring)
			max_substring = substring
	print(max_substring)

if __name__ == '__main__':
    sys.exit(main())
