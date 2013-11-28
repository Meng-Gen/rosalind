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

def is_reverse_palindrome(dna):
    return dna == dna_reverse_complement(dna)
    
def dna_reverse_complement(dna):
    c = dna.replace("A", "a").replace("T", "A").replace("a", "T")
    c = c.replace("C", "c").replace("G", "C").replace("c", "G")
    return c[::-1]

def locate_restriction_sites(dna):
    n = len(dna)
    for i in range(n):
        for j in range(4, 12 + 1):
            if i + j > n:
                break
            if is_reverse_palindrome(dna[i:i+j]):
                yield i + 1, j
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    dna = fasta_records[0][1]
    for position, length in locate_restriction_sites(dna):
        print(position, length)

if __name__ == '__main__':
    sys.exit(main())
