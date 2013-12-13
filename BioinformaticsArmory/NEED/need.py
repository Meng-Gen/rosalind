import sys

from Bio import Entrez
from Bio import SeqIO

def get_records(ids):
    Entrez.email = 'plover@gmail.com'
    query_string = ', '.join(ids)
    handle = Entrez.efetch(db="nucleotide", id=[', '.join(ids)], rettype="fasta")
    return list(SeqIO.parse(handle, "fasta"))
    
def main():
    ids = sys.stdin.read().strip().split()
    for record in get_records(ids):
        print(record.seq)

    # Then go to http://www.ebi.ac.uk/Tools/psa/emboss_needle/nucleotide.html
    #
    # STEP 1 - Enter your nucleotide sequences
    # Easy
    #
    # STEP 2 - Set your pairwise alignment options
    # MATRIX: DNAfull
    # GAP OPEN: 10
    # GAP EXTEND: 1.0
    # OUTPUT FORMAT: pair
    # END GAP PENALTY: true
    # END GAP OPEN: 10
    # END GAP EXTEND: 1.0
        
if __name__ == '__main__':
    sys.exit(main())
