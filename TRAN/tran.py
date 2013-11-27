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

def num_transitions(s, t):
    assert(len(s) == len(t))
    transitions = ['AG', 'GA', 'CT', 'TC']
    return sum([1 if ''.join(_) in transitions else 0 for _ in zip(s, t)])

def num_transversions(s, t):
    assert(len(s) == len(t))
    transversions = ['AC', 'AT', 'GC', 'GT', 'CA', 'CG', 'TA', 'TG']
    return sum([1 if ''.join(_) in transversions else 0 for _ in zip(s, t)])
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    s = fasta_records[0][1]
    t = fasta_records[1][1]
    print(num_transitions(s, t)/num_transversions(s, t))

if __name__ == '__main__':
    sys.exit(main())
