import sys

class SpectrumGraph():
    def __init__(self, spectrum_records):
        self.MONOISOTOPIC_MASS_TABLE = {
            'A': 71.03711,
            'C': 103.00919,
            'D': 115.02694,
            'E': 129.04259,
            'F': 147.06841,
            'G': 57.02146,
            'H': 137.05891,
            'I': 113.08406,
            'K': 128.09496,
            'L': 113.08406,
            'M': 131.04049,
            'N': 114.04293,
            'P': 97.05276,
            'Q': 128.05858,
            'R': 156.10111,
            'S': 87.03203,
            'T': 101.04768,
            'V': 99.06841,
            'W': 186.07931,
            'Y': 163.06333,
        }
        self.adjacency_matrix = None 
        self.num_nodes = None      
        self.protein_strings = set()
        
        self.__build_graph(spectrum_records)
    
    def get_longest_protein_string(self):
        for begin_pos in range(self.num_nodes):
            self.__get_protein_strings(begin_pos, '')
        longest_len_so_far, longest_protein_so_far = 0, None
        for protein in self.protein_strings:
            if len(protein) > longest_len_so_far:
                longest_len_so_far = len(protein)
                longest_protein_so_far = protein
        return longest_protein_so_far
    
    def __get_protein_strings(self, begin_pos, protein_string):
        self.protein_strings.add(protein_string)
        for end_pos in range(self.num_nodes):
            protein = self.adjacency_matrix[begin_pos][end_pos]
            if not protein:
                continue
            self.__get_protein_strings(end_pos, protein_string + protein)
                
    def __build_graph(self, spectrum_records):
        spectrum_records.sort()
        n = len(spectrum_records)
        self.adjacency_matrix = [[None for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                mass = spectrum_records[j] - spectrum_records[i]
                protein = self.__get_possible_protein(mass)
                if protein:
                    self.adjacency_matrix[i][j] = protein
        self.num_nodes = n
        
    def __get_possible_protein(self, mass):
        for protein in self.MONOISOTOPIC_MASS_TABLE:
            if abs(self.MONOISOTOPIC_MASS_TABLE[protein] - mass) < 1e-4:
                return protein
        return None
    
def main():
    spectrum_records = [float(_.strip()) for _ in sys.stdin.readlines()]
    graph = SpectrumGraph(spectrum_records)
    print(graph.get_longest_protein_string())
    
if __name__ == '__main__':
    sys.exit(main())
