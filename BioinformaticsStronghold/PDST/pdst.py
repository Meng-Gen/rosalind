import sys

def get_fasta_records(lines):
    dna_strings = []
    is_first_dna = True
    id, dna = '', ''
    for line in lines:
        line = line.strip();
        if line[0] == '>':
            if is_first_dna:
                is_first_dna = False
            else:
                dna_strings.append([id, dna])
            id, dna = line[1:], ''
        else:
            dna += line
    dna_strings.append([id, dna])
    return dna_strings

def p_distance(s, t):
    num_diff = sum(1 if _[0] != _[1] else 0 for _ in zip(s, t))
    return num_diff / len(s)
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    for s in fasta_records:
        row_p_distance = []
        for t in fasta_records:
            row_p_distance.append(p_distance(s[1], t[1]))
        print(' '.join(map(str, row_p_distance)))
            
if __name__ == '__main__':
    sys.exit(main())
