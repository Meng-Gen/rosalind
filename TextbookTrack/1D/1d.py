import sys

def is_clump(kmer, L, t, genome):
    k = len(kmer)
    all_positions = []
    for i in range(len(genome)-k+1):
        if kmer == genome[i:i+k]:
            all_positions.append(i)
    pos_count = len(all_positions)
    if pos_count < t:
        return False
    for i in range(pos_count-t+1):
        if all_positions[i+t-1] - all_positions[i] + k - 1 <= L:
            return True
    return False

def main():
    genome = sys.stdin.readline().strip()
    k, L, t = map(int, sys.stdin.readline().strip().split())
    clumps = set()
    visited_kmer = set()
    for i in range(len(genome)-k+1):
        kmer = genome[i:i+k]
        # do not check duplicated k-mers
        if kmer in visited_kmer:
            continue
        visited_kmer.add(kmer)
        if is_clump(kmer, L, t, genome):
            clumps.add(kmer)
    print(' '.join(list(clumps)))
    
if __name__ == '__main__':
    sys.exit(main())
