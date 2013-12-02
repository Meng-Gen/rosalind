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
    return dna.replace('T', 'U')

def reverse_complement_dna(dna):
    return dna.replace("A", "a").replace("T", "A").replace("a", "T").replace("C", "c").replace("G", "C").replace("c", "G")[::-1]
    
def rna_codon(rna):
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
    valid_rna = rna[:len(rna) - (len(rna) % 3)]
    protein = []
    for i in range(0, len(valid_rna), 3):
        protein.append(rna_codon_table[valid_rna[i:i+3]])
    return protein

def candidate_proteins(dna):
    rna = dna_to_rna(dna)
    reverse_complement_rna = dna_to_rna(reverse_complement_dna(dna))
    return [rna_codon(rna[_:]) for _ in range(3)] + [rna_codon(reverse_complement_rna[_:]) for _ in range(3)]
    
def open_reading_frames(protein):
    start_codon_pos = [i for i in range(len(protein)) if protein[i] == 'M']
    stop_codon_pos = [i for i in range(len(protein)) if protein[i] == 'Stop']
    orfs = set()
    visited_i = set()
    for j in stop_codon_pos:
        for i in start_codon_pos:
            if i > j:
                break
            if i in visited_i:
                continue
            visited_i.add(i)
            orfs.add(''.join(protein[i:j]))
    return orfs    
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    dna = fasta_records[0][1]
    
    all_orfs = set()
    for protein in candidate_proteins(dna):
        all_orfs |= open_reading_frames(protein)
    print('\n'.join(list(all_orfs)))

if __name__ == '__main__':
    sys.exit(main())
