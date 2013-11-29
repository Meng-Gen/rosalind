import sys

dna = sys.stdin.read().strip()
print(' '.join([str(dna.count(i)) for i in 'ACGT']))
