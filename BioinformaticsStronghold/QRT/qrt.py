import itertools
import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    taxa = lines[0].split()
    partial_character_table = lines[1:]
    return taxa, partial_character_table
    
class QuartetBuilder():
    def __init__(self, partial_character_table, taxa):
        self.partial_character_table = partial_character_table
        self.taxa = taxa
        
        self.id_quartet_set = {}
        self.taxa_quartet_strings = []
        
    def build(self):
        for row in self.partial_character_table:
            self.__build_id_quartets_from_row(row)
        self.__build_taxa_quartets()
        return self.taxa_quartet_strings
    
    def __build_id_quartets_from_row(self, row):
        n = len(self.partial_character_table[0])
        on_set, off_set = set(), set()
        for i in range(n):
            if row[i] == '0':
                off_set.add(i)
            elif row[i] == '1':
                on_set.add(i)
        for a in itertools.combinations(on_set, 2):
            for b in itertools.combinations(off_set, 2):
                self.__insert_quartet(a, b)
    
    def __insert_quartet(self, a, b):
        a_array, b_array = list(a), list(b)
        a_array.sort()
        b_array.sort()
        if a_array[0] > b_array[0]:
            a_array, b_array = b_array, a_array
        hash = '''%d-%d-%d-%d''' % (a_array[0], a_array[1], b_array[0], b_array[1])
        if hash in self.id_quartet_set:
            return
        self.id_quartet_set[hash] = [a, b]
        
    def __build_taxa_quartets(self):
        for id_quartet in self.id_quartet_set:
            a, b = self.id_quartet_set[id_quartet]
            quartet_string = '''{%s, %s} {%s, %s}''' % \
                    (self.taxa[a[0]], self.taxa[a[1]], self.taxa[b[0]], self.taxa[b[1]])
            self.taxa_quartet_strings.append(quartet_string)
        
def main():
    taxa, partial_character_table = read_dataset()
    builder = QuartetBuilder(partial_character_table, taxa)
    print('\n'.join(builder.build()))
    
if __name__ == '__main__':
    sys.exit(main())
