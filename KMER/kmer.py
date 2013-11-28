import itertools
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

def kmer_composition(dna, k):
    freq = {}
    for kmer in itertools.product('ACGT', repeat=k):
        freq[''.join(kmer)] = 0
    for i in range(0, len(dna) - k + 1):
        freq[dna[i:i+4]] += 1
    return [freq[''.join(_)] for _ in itertools.product('ACGT', repeat=k)]
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    dna = fasta_records[0][1]
    print(' '.join(map(str, kmer_composition(dna, 4))))

if __name__ == '__main__':
    sys.exit(main())
