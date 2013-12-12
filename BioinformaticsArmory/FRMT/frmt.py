import sys

from Bio import Entrez
from Bio import SeqIO

def get_records(ids):
    Entrez.email = 'plover@gmail.com'
    query_string = ', '.join(ids)
    handle = Entrez.efetch(db="nucleotide", id=[', '.join(ids)], rettype="fasta")
    return list(SeqIO.parse(handle, "fasta"))

def get_shortest_record(records):
    shortest_len_so_far = len(records[0].seq)
    shortest_fasta_pos_so_far = 0
    for i in range(1, len(records)):
        curr_len = len(records[i].seq)
        if curr_len < shortest_len_so_far:
            shortest_len_so_far = curr_len
            shortest_fasta_pos_so_far = i
    return records[shortest_fasta_pos_so_far]
    
def main():
    ids = sys.stdin.read().strip().split()
    records = get_records(ids)
    record = get_shortest_record(records)
    print(record.format("fasta"))
    
if __name__ == '__main__':
    sys.exit(main())
