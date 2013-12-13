import sys

def get_translation_table(coding_dna, protein):
    from Bio.Seq import translate
    available_table = [_ for _ in range(1, 7)] + [_ for _ in range(9, 16)]
    for i in available_table:
        if translate(coding_dna, table=i, to_stop=True) == protein:
            return i
            
def main():
    coding_dna, protein = [line.strip() for line in sys.stdin.readlines()]
    print(get_translation_table(coding_dna, protein))
    
if __name__ == '__main__':
    sys.exit(main())
