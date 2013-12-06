import sys

class CatalanNumberFactory():
    def __init__(self):
        self.memorized_table = { 
            'AU' : 1,
            'UA' : 1,
            'CG' : 1,
            'GC' : 1,
        }
        
    def create(self, rna):
        return self.__get_catalan_number(rna)
        
    def __get_catalan_number(self, rna):
        if rna in self.memorized_table:
            return self.memorized_table[rna]
            
        counts = [rna.count(i) for i in 'ACGU']
        if counts[0] != counts[3] or counts[1] != counts[2]:
            return 0
            
        n = len(rna)
        sum = 0
        connected_edge = rna[:2]
        if connected_edge in self.memorized_table:
            rest_rna = rna[2:]
            sum += self.__get_catalan_number(rest_rna)
        connected_edge = rna[n-1] + rna[0]
        if connected_edge in self.memorized_table:
            rest_rna = rna[1:n-1]
            sum += self.__get_catalan_number(rest_rna)
        for i in range(2, n-1):
            connected_edge = rna[0] + rna[i]
            if connected_edge in self.memorized_table:
                rest_rna_pieces = [rna[1:i], rna[i+1:]]
                pieces_numbers = [self.__get_catalan_number(_) for _ in rest_rna_pieces]
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
    factory = CatalanNumberFactory()
    num = factory.create(rna)
    print(num % 1000000)
    
if __name__ == '__main__':
    sys.exit(main())
