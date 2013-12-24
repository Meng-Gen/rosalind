import itertools
import sys

def read_dataset():
    lines = sys.stdin.readlines()
    num_protein_strings = int(lines[0].strip())
    proteins = [lines[i].strip() for i in range(1, num_protein_strings + 1)]
    unknown_spectrums = [float(lines[i].strip()) for i in range(num_protein_strings + 1, len(lines))]
    return proteins, unknown_spectrums
   
class SpectrumMatcher():
    def __init__(self, proteins, unknown_spectrums):
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
        self.unknown_spectrums = unknown_spectrums
        self.proteins = proteins
    
    def match(self):
        max_multiplicity_so_far = 0
        max_protein_so_far = None
        for protein in self.proteins:
            spectrums = self.__get_complete_spectrums(protein)
            difference = self.__minkowski_difference(self.unknown_spectrums, spectrums)
            frequency = self.__convert_to_frequency(difference)
            max_multiplicity = self.__get_max_multiplicity(frequency)
            if max_multiplicity > max_multiplicity_so_far:
                max_multiplicity_so_far = max_multiplicity
                max_protein_so_far = protein
        return max_multiplicity_so_far, max_protein_so_far
                
    def __get_complete_spectrums(self, protein):
        mass_array = [self.MONOISOTOPIC_MASS_TABLE[c] for c in protein]
        # prefix
        spectrums = []
        acc_mass = 0.0
        for mass in mass_array[:-1]:
            acc_mass += mass
            spectrums.append(acc_mass)
        # suffix
        acc_mass = 0.0
        for mass in mass_array[::-1]:
            acc_mass += mass
            spectrums.append(acc_mass)
        return spectrums

    def __minkowski_difference(self, s, t):
        return [x-y for x, y in itertools.product(s, t)]

    def __convert_to_frequency(self, multiset):
        frequency = {}
        for x in multiset:
            x_literal = "%0.5f" % x # hack for performance
            if x_literal not in frequency:
                frequency[x_literal] = 0
            frequency[x_literal] += 1
        return frequency        
        
    def __get_max_multiplicity(self, frequency):
        max_so_far = 0
        for protein in frequency:
            mass = frequency[protein]
            if mass > max_so_far:
                max_so_far = mass
        return max_so_far        
    
def main():
    proteins, unknown_spectrums = read_dataset()
    matcher = SpectrumMatcher(proteins, unknown_spectrums)
    print('\n'.join(map(str, matcher.match())))
    
if __name__ == '__main__':
    sys.exit(main())
