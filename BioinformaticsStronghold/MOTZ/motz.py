import sys

class MotzkinNumberFactory():
    def __init__(self):
        self.memorized_table = { 
            '' : 1,
            'A' : 1,
            'C' : 1,
            'G' : 1,
            'U' : 1,
            'AU' : 2, 
            'UA' : 2,
            'CG' : 2,
            'GC' : 2,
        }
        self.vaild_edges = { 
            'AU', 
            'UA',
            'CG',
            'GC',
        }
        
    def create(self, rna):
        return self.__get_number(rna)
        
    def __get_number(self, rna):
        if rna in self.memorized_table:
            return self.memorized_table[rna]
            
        n = len(rna)
        sum = self.__get_number(rna[1:])
        connected_edge = rna[:2]
        if connected_edge in self.vaild_edges:
            rest_rna = rna[2:]
            sum += self.__get_number(rest_rna)
        connected_edge = rna[n-1] + rna[0]
        if connected_edge in self.vaild_edges:
            rest_rna = rna[1:n-1]
            sum += self.__get_number(rest_rna)
        for i in range(2, n-1):
            connected_edge = rna[0] + rna[i]
            if connected_edge in self.vaild_edges:
                rest_rna_pieces = [rna[1:i], rna[i+1:]]
                pieces_numbers = [self.__get_number(_) for _ in rest_rna_pieces]
                sum += pieces_numbers[0] * pieces_numbers[1]        
        self.memorized_table[rna] = sum
        self.memorized_table[rna[::-1]] = sum         
        return sum
        
def get_fasta_records(lines):
    dna_strings = []
    is_first_dna = True
    id, dna = '', ''
    for line in lines:
        line = line.strip();
        if line[0] == '>':
            if is_first_dna:
                is_first_dna = False
            else:
                dna_strings.append([id, dna])
            id, dna = line[1:], ''
        else:
            dna += line
    dna_strings.append([id, dna])
    return dna_strings  
    
def main():
    fasta_records = get_fasta_records(sys.stdin.readlines())
    rna = fasta_records[0][1]
    factory = MotzkinNumberFactory()
    num = factory.create(rna)
    print(num % 1000000)
    
if __name__ == '__main__':
    sys.exit(main())
