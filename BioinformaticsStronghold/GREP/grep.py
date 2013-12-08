import sys

class GenomeAssembly():
    def __init__(self, dna_strings):
        self.dna_strings = dna_strings
        self.all_genome_strings = set()
        self.dna_count = len(dna_strings)
        self.k = len(dna_strings[0])
        self.solution_vector = []
        self.visited_edges = [False for _ in dna_strings]
        
    def assemble(self):
        self.solution_vector.append(self.dna_strings[0])
        self.visited_edges[0] = True
        self.backtrack(1)
        return self.all_genome_strings
        
    def backtrack(self, dimension):
        if not self.is_valid_assembly():
            return
        if len(self.solution_vector) == self.dna_count:
            self.build_genome_string()
            return 
        for i in range(len(self.dna_strings)):
            if self.visited_edges[i] is True:
                continue
            self.solution_vector.append(self.dna_strings[i])
            self.visited_edges[i] = True
            self.backtrack(dimension + 1)
            self.visited_edges[i] = False
            self.solution_vector.pop(-1)
        pass
        
    def is_valid_assembly(self):
        n = len(self.solution_vector)
        if n == 1:
            return True
        suffix = self.solution_vector[-2][-self.k+1:]
        prefix = self.solution_vector[-1][:self.k-1]
        if suffix != prefix:
            return False
        if n == self.dna_count:
            suffix = self.solution_vector[-1][-self.k+1:]
            prefix = self.solution_vector[0][:self.k-1]
            if suffix != prefix:
                return False
        return True

    def build_genome_string(self):
        genome_string = self.solution_vector[0]
        for i in range(1, self.dna_count - self.k + 1):
            genome_string += self.solution_vector[i][-1]
        self.all_genome_strings.add(genome_string)
      
def main():
    dna_strings = [dna.strip() for dna in sys.stdin.readlines()]
    assembly = GenomeAssembly(dna_strings)
    for genome_string in assembly.assemble():
        print(genome_string)
    
if __name__ == '__main__':
    sys.exit(main())
