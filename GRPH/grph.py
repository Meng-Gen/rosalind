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

def is_adjacent(s, t, overlap_length):
    return s != t and s[-overlap_length:] == t[:overlap_length]
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    for s in fasta_records:
        for t in fasta_records:
            if is_adjacent(s[1], t[1], 3):
                print(s[0], t[0])

if __name__ == '__main__':
    sys.exit(main())
