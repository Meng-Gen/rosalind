import math
import sys

def get_fasta_records(lines):
    dna_strings = []
    is_first_rna = True
    id, rna = '', ''
    for line in lines:
        line = line.strip();
        if line[0] == '>':
            if is_first_rna:
                is_first_rna = False
            else:
                dna_strings.append([id, rna])
            id, rna = line[1:], ''
        else:
            rna += line
    dna_strings.append([id, rna])
    return dna_strings

def num_perfect_matchings(rna):
    symbol_count = [rna.count(i) for i in 'ACGU']
    return math.factorial(symbol_count[0]) * math.factorial(symbol_count[1])
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    rna = fasta_records[0][1]
    print(num_perfect_matchings(rna))

if __name__ == '__main__':
    sys.exit(main())
