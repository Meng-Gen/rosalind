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
    
def dna_to_rna(dna):
    return dna.replace("T", "U")

def rna_to_protein(rna):
    rna_codon_table = {
        'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
        'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
        'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
        'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
        'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
        'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
        'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
        'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
        'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
        'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
        'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
        'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
        'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
        'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
        'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
        'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G', 
    }
    protein = []
    for i in range(0, len(rna), 3):
        if rna_codon_table[rna[i:i+3]] == 'Stop':
            break
        else:
            protein.append(rna_codon_table[rna[i:i+3]])
    return protein
            
def dna_to_protein(dna):
    return rna_to_protein(dna_to_rna(dna))

def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    dna = fasta_records[0][1]
    for i in range(1, len(fasta_records)):
        dna = dna.replace(fasta_records[i][1], '')
    print(''.join(dna_to_protein(dna)))

if __name__ == '__main__':
    sys.exit(main())
