import sys

MONOISOTOPIC_MASS_TABLE = {
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

def get_prefix_spectrum_records():
    records = [float(_.strip()) for _ in sys.stdin.readlines()]
    records.sort()
    return records

def get_possible_protein(weight):
    for protein in MONOISOTOPIC_MASS_TABLE:
        if abs(MONOISOTOPIC_MASS_TABLE[protein] - weight) < 1e-4:
            return protein
        
def main():
    records = get_prefix_spectrum_records()
    protein_string = ''
    for i in range(1, len(records)):
        weight = records[i] - records[i-1]
        protein = get_possible_protein(weight)
        protein_string += protein
    print(protein_string)
    
if __name__ == '__main__':
    sys.exit(main())
