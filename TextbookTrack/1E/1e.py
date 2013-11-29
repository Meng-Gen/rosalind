import itertools
import sys

def get_skew_prefix(genome):
    contributions = [0]
    for nucleotide in genome:
        if nucleotide == 'C':
            contributions.append(-1)
        elif nucleotide == 'G':
            contributions.append(1)
        else:
            contributions.append(0)
    return list(itertools.accumulate(contributions))

def get_min_skew_pos(skew_prefix):
    min_skew = min(skew_prefix)
    return [i for i in range(len(skew_prefix)) if skew_prefix[i] == min_skew]
    
def main():
    genome = sys.stdin.readline().strip()
    skew_prefix = get_skew_prefix(genome)
    min_skew_pos = get_min_skew_pos(skew_prefix)
    print(' '.join(map(str, min_skew_pos)))
    
if __name__ == '__main__':
    sys.exit(main())
