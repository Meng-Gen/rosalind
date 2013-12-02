import sys

def is_adjacent(s, t, overlap_length):
    return s != t and s[-overlap_length:] == t[:overlap_length]

def dna_rc(dna):
    rc = dna.replace("A", "a").replace("T", "A").replace("a", "T")
    rc = rc.replace("C", "c").replace("G", "C").replace("c", "G")
    return rc[::-1]

def main():
    dna_records = [dna.strip() for dna in sys.stdin.readlines()]
    dna_rc_records = [dna_rc(dna) for dna in dna_records]
    for dna in set(dna_records + dna_rc_records):
        print('({begin}, {end})'.format(begin=dna[:-1], end=dna[1:]))    

if __name__ == '__main__':
    sys.exit(main())
