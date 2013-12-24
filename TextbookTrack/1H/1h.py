import itertools
import sys

class ApproximateMatchingProblem():
    def __init__(self, genome, k, d):
        self.genome = genome
        self.k = k
        self.d = d
        self.approximating_motifs = None
    
    def get_most_frequent_motifs(self):
        self.approximating_motifs = {}
        genome_rc = self.__get_reverse_complement(self.genome)  
        m = len(self.genome) - self.k + 1
        for i in range(m):
            piece = self.genome[i:i+self.k]
            self.__update_approximating_motifs(piece)
            piece = genome_rc[i:i+self.k]
            self.__update_approximating_motifs(piece)
        return self.__get_most_frequent_motifs_impl()
    
    def __get_reverse_complement(self, genome):
        complement = genome[::-1]
        complement = complement.replace("A", "a").replace("T", "A").replace("a", "T")
        complement = complement.replace("C", "c").replace("G", "C").replace("c", "G")
        return complement
    
    def __update_approximating_motifs(self, piece):
        for pattern in itertools.product([False, True], repeat=self.k):
            count = pattern.count(False)
            if count > self.d:
                continue
            approximating_patterns = []
            for i in range(self.k):
                if not pattern[i]:
                    possible_nucleotides = list('ACGT')
                    possible_nucleotides.remove(piece[i])
                    approximating_patterns.append(possible_nucleotides)
                else:
                    approximating_patterns.append([piece[i]])
            for x in itertools.product(*approximating_patterns):
                motif = ''.join(x)
                if motif not in self.approximating_motifs:
                    self.approximating_motifs[motif] = 0
                self.approximating_motifs[motif] += 1
    
    def __get_most_frequent_motifs_impl(self):
        max_freq_so_far = 0
        max_freq_motifs_so_far = set()
        for motif in self.approximating_motifs:
            freq = self.approximating_motifs[motif]
            if freq > max_freq_so_far:
                max_freq_so_far = freq
                max_freq_motifs_so_far = set()
                max_freq_motifs_so_far.add(motif)
            elif freq == max_freq_so_far:
                max_freq_motifs_so_far.add(motif)
        return max_freq_motifs_so_far
        
def main():
    genome = sys.stdin.readline().strip()
    k, d = map(int, sys.stdin.readline().strip().split())
    problem = ApproximateMatchingProblem(genome, k, d)
    print(' '.join(list(problem.get_most_frequent_motifs())))
    
if __name__ == '__main__':
    sys.exit(main())
