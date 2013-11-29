import sys

def get_permutation_pairs(lines):
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
    pairs = get_permutation_pairs(sys.stdin.readlines())
    print(pairs)

if __name__ == '__main__':
    sys.exit(main())
