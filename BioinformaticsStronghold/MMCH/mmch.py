import math
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
 
def get_num_maximum_matchings(rna):
    au_pair = [rna.count(_) for _ in 'AU']
    cg_pair = [rna.count(_) for _ in 'CG']
    return get_num_pair_matchings(au_pair) * get_num_pair_matchings(cg_pair)

def get_num_pair_matchings(pair):
    pair.sort()
    return math.factorial(pair[1]) // math.factorial(pair[1] - pair[0])    
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    rna = fasta_records[0][1]
    print(get_num_maximum_matchings(rna))
    
if __name__ == '__main__':
    sys.exit(main())
