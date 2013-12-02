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

def kmp_failure_array(s):
    m, i, j = len(s), 1, 0
    p = [0] * m
    while i < m:
        if s[j] == s[i]:
            p[i] = j + 1
            i, j = i+1, j+1
        elif j != 0:
            j = p[j-1]
        else: 
            i += 1
    return p
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    dna = fasta_records[0][1]
    print(' '.join(list(map(str, kmp_failure_array(dna)))))
            
if __name__ == '__main__':
    sys.exit(main())
