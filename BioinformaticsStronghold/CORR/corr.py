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

def reverse_complement_dna(dna):
    return dna.replace("A", "a").replace("T", "A").replace("a", "T").replace("C", "c").replace("G", "C").replace("c", "G")[::-1]

def hamming_distance(s, t):
    return len([_ for _ in zip(s, t) if _[0] != _[1]])
    
def get_corrections(dna_records):
    reverse_complement_dna_records = [reverse_complement_dna(_) for _ in dna_records]
    total_records = dna_records + reverse_complement_dna_records
    correct_dna_set = set()
    for correct_dna in total_records:
        is_correct = len([_ for _ in total_records if _ == correct_dna]) >= 2
        if is_correct:
            correct_dna_set.add(correct_dna)
    incorrect_dna_set = set()
    for wrong_dna in dna_records:
        if wrong_dna in correct_dna_set:
            continue
        wrt_correct_dna_list = [] # w.r.t = with respect to 
        for correct_dna in correct_dna_set:
            if hamming_distance(wrong_dna, correct_dna) == 1:
                wrt_correct_dna_list.append(correct_dna)
        if len(wrt_correct_dna_list) == 1:
            yield wrong_dna, wrt_correct_dna_list[0]
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    dna_records = [_[1] for _ in fasta_records]
    for correction in get_corrections(dna_records):
        print('%s->%s' % (correction[0], correction[1]))
            
if __name__ == '__main__':
    sys.exit(main())
