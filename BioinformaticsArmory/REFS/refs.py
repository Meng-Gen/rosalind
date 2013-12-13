import sys

def read_dataset():
    return [line.strip() for line in sys.stdin.readlines()]
    
def get_count(organism, seq_len_from, seq_len_to, publication_date_to):
    from Bio import Entrez
    Entrez.email = 'plover@gmail.com'
    query_string = '''"%s"[Organism] AND %s:%s[Sequence Length] AND srcdb_refseq[Properties] AND 1986/1/1:%s[Publication Date]''' \
            % (organism, seq_len_from, seq_len_to, publication_date_to)
    handle = Entrez.esearch(db="nucleotide", term=query_string)
    record = Entrez.read(handle)
    return record["Count"]
    
def main():
    organism, seq_len_from, seq_len_to, publication_date_to = read_dataset()
    print(get_count(organism, seq_len_from, seq_len_to, publication_date_to))
    
if __name__ == '__main__':
    sys.exit(main())
