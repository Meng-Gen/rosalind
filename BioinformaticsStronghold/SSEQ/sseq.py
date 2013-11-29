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

def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    s = fasta_records[0][1]
    t = fasta_records[1][1]

    sequence_indexes = [] # 1-based
    j = 0
    for i in range(len(s)):
        if j == len(t):
            break
        if s[i] == t[j]:
            sequence_indexes.append(i + 1)
            j += 1
    print(' '.join(map(str, sequence_indexes)))
            
if __name__ == '__main__':
    sys.exit(main())
