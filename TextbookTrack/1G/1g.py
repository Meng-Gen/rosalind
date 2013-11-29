import itertools
import sys

def all_possible_matching_kmers(genome, k, d):
    all_possible_kmers = set()
    for i in range(len(genome)-k+1):
        motif = genome[i:i+k]
        kmers = approximating_matching_kmers(motif, k, d)
        all_possible_kmers |= kmers
        print(kmers)
    return all_possible_kmers
    
def approximating_matching_kmers(motif, k, d):
    all_kmers = set()
    all_mismatch_patterns = set()
    for mismatch_pattern in itertools.permutations(['D']*d + ['M']*(k-d), k):
        all_mismatch_patterns.add(''.join(mismatch_pattern))
    for pattern in all_mismatch_patterns:
        for nucleotides in itertools.product('ACGT', repeat=d):
            kmer = ''
            mismatch_count = 0
            for i in range(k):
                if pattern[i] == 'M':
                    kmer += motif[i]
                else:
                    kmer += nucleotides[mismatch_count]
                    mismatch_count += 1
            all_kmers.add(kmer)
    return all_kmers

def approximating_matching_count(genome, kmer, d):
    count = 0
    m, n = len(kmer), len(genome)
    for i in range(n):
        motif = genome[i:i+m]
        if len(motif) == m and is_approximating_matching(motif, kmer, d):
            count += 1
    return count
    
def is_approximating_matching(p, q, d):
    return sum([1 if p[i] != q[i] else 0 for i in range(len(p))]) <= d
    
def main():
    genome = sys.stdin.readline().strip()
    k, d = map(int, sys.stdin.readline().strip().split())
    
    max_count, max_count_kmers = 0, []
    for kmer in all_possible_matching_kmers(genome, k, d):
        count = approximating_matching_count(genome, kmer, d)
        print(kmer, count)
        if count == max_count:
            max_count_kmers.append(kmer)
        elif count > max_count:
            max_count, max_count_kmers = count, [kmer]
    print(' '.join(max_count_kmers))

if __name__ == '__main__':
    sys.exit(main())
