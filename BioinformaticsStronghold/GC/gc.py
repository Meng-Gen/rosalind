import sys

def get_dna_strings(lines):
	dna_strings = {}
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
	return dna_strings

def get_gc_content_percentage(dna):
	return sum([dna.count(i) for i in 'GC']) * 100 / len(dna)
	
def main():
	dna_strings = get_dna_strings(sys.stdin.readlines())
	max_id_so_far = ''
	max_so_far = 0
	is_first_id = True
	for id in dna_strings:
		gc_content = get_gc_content_percentage(dna_strings[id])
		if is_first_id:
			max_id_so_far = id
			max_so_far = gc_content
			is_first_id = False
		elif gc_content > max_so_far:
			max_id_so_far = id
			max_so_far = gc_content
	print(max_id_so_far)
	print('%.6f' % max_so_far)

if __name__ == '__main__':
    sys.exit(main())
