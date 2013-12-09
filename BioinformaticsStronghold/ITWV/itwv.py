import itertools
import sys

class DisjointMotifsFinder():
    def __init__(self, gene, motifs):
        self.gene = gene
        self.motifs = motifs    
        self.solution_matrix = [[None for _ in motifs] for _ in motifs]
        
    def find(self):
        for i, j in itertools.product(range(len(self.motifs)),repeat=2):
            self.solution_matrix[i][j] = self.are_interwoven(self.motifs[i], self.motifs[j])
        return self.solution_matrix
        
    def are_interwoven(self, t, u):
        n = len(t) + len(u)
        for i in range(len(self.gene) - n + 1):
            if self.are_interwoven_exact(self.gene[i:i+n], t, u):
                return True
        return False
        
    def are_interwoven_exact(self, s, t, u):
        subproblem_table = [[None for j in range(len(u) + 1)] for i in range(len(t) + 1)]
        subproblem_table[0][0] = True
        for i in range(len(t)):
            subproblem_table[i+1][0] = subproblem_table[i][0] and (t[i] == s[i])
        for j in range(len(u)):
            subproblem_table[0][j+1] = subproblem_table[0][j] and (u[j] == s[j])
        
        for i, j in itertools.product(range(len(t)), range(len(u))):
            subproblem_table[i+1][j+1] = subproblem_table[i][j+1] and (t[i] == s[i+j+1])
            subproblem_table[i+1][j+1] |= subproblem_table[i+1][j] and (u[j] == s[i+j+1])
        return subproblem_table[len(t)][len(u)]
        
def main():
    dna_strings = [dna.strip() for dna in sys.stdin.readlines()]
    gene = dna_strings[0]
    motifs = dna_strings[1:]
    finder = DisjointMotifsFinder(gene, motifs)
    for matrix_row in finder.find():
        print(' '.join(['1' if _ is True else '0' for _ in matrix_row]))
    
if __name__ == '__main__':
    sys.exit(main())
