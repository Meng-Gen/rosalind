import sys

class RigidMotzkinNumberFactory():
    def __init__(self):
        self.rigid_threshold = 4
        self.memorized_table = {
            '' : 1, 
        }
        self.vaild_edges = { 
            'AU', # base pair
            'UA', # base pair
            'CG', # base pair
            'GC', # base pair
            'GU', # wobble base pair
            'UG', # wobble base pair
        }
        
    def create(self, rna):
        return self.__get_number(rna)
        
    def __get_number(self, rna):
        if rna in self.memorized_table:
            return self.memorized_table[rna]
        n = len(rna)
        sum = self.__get_number(rna[1:])
        for i in range(self.rigid_threshold, n):
            connected_edge = rna[0] + rna[i]
            if connected_edge not in self.vaild_edges:
                continue
            rest_rna_pieces = [rna[1:i], rna[i+1:]]
            pieces_numbers = [self.__get_number(_) for _ in rest_rna_pieces]
            sum += pieces_numbers[0] * pieces_numbers[1]
        self.memorized_table[rna] = sum
        self.memorized_table[rna[::-1]] = sum
        return sum
    
def main():
    rna = sys.stdin.read().strip()
    factory = RigidMotzkinNumberFactory()
    num = factory.create(rna)
    print(num)
    
if __name__ == '__main__':
    sys.exit(main())
