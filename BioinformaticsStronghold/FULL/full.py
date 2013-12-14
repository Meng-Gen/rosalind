import sys

def get_records():
    dataset = [float(_.strip()) for _ in sys.stdin.readlines()]
    parent_mass = dataset[0]
    spectrum_records = dataset[1:]
    return parent_mass, spectrum_records

class SpectrumAnalyzer():
    def __init__(self):
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
        self.parent_mass = None
        self.spectrum_records = None
        
        self.n = None
        self.protein_len = None
        self.solution_spectrum_vector = None
        self.solution_protein_vector = None
        self.used_records = None

    def analyze(self, parent_mass, spectrum_records):
        self.__init_analyzing(parent_mass, spectrum_records)

        self.solution_spectrum_vector = [self.spectrum_records[0]]
        self.solution_protein_vector = []
        self.used_records[0] = True 
        self.used_records[-1] = True      
        self.__backtrack()   
    
    def __init_analyzing(self, parent_mass, spectrum_records):
        self.parent_mass = parent_mass
        self.spectrum_records = spectrum_records
        self.spectrum_records.sort()
        self.n = len(spectrum_records)
        self.protein_len = (self.n - 1)//2
        self.solution_spectrum_vector = []
        self.solution_protein_vector = []
        self.used_records = [False for _ in self.spectrum_records]
    
    def __backtrack(self):        
        if len(self.solution_protein_vector) == self.protein_len:            
            print(''.join(self.solution_protein_vector))
            sys.exit()
        
        curr_spectrum = self.solution_spectrum_vector[-1]
        for protein in self.MONOISOTOPIC_MASS_TABLE:
            mass = self.MONOISOTOPIC_MASS_TABLE[protein]
            expected_spectrum = mass + curr_spectrum
            
            pos = self.__get_actual_spectrum_pos(expected_spectrum)
            if not pos:
                continue
            if self.used_records[pos] is True:
                continue
            self.solution_spectrum_vector.append(expected_spectrum)
            self.solution_protein_vector.append(protein)
            self.used_records[pos] = True
            self.used_records[self.n - 1 - pos] = True
            self.__backtrack()
            self.used_records[self.n - 1 - pos] = False
            self.used_records[pos] = False
            self.solution_protein_vector.pop(-1)
            self.solution_spectrum_vector.pop(-1)
            
    def __get_actual_spectrum_pos(self, expected_spectrum):
        for i in range(self.n):
            if abs(expected_spectrum - self.spectrum_records[i]) < 1e-4:
                return i
        return None
        
def main():
    parent_mass, spectrum_records = get_records()
    analyzer = SpectrumAnalyzer()
    analyzer.analyze(parent_mass, spectrum_records)
    
if __name__ == '__main__':
    sys.exit(main())
