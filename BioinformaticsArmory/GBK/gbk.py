import sys

def read_dataset():
    return [line.strip() for line in sys.stdin.readlines()]
    
def get_count(organism, publication_begin_date, publication_end_date):
    from Bio import Entrez
    Entrez.email = 'plover@gmail.com'
    query_string = '''"%s"[Organism] AND %s:%s[Publication Date]''' \
            % (organism, publication_begin_date, publication_end_date)
    handle = Entrez.esearch(db="nucleotide", term=query_string)
    record = Entrez.read(handle)
    return record["Count"]
    
def main():
    organism, publication_begin_date, publication_end_date = read_dataset()
    print(get_count(organism, publication_begin_date, publication_end_date))
    
if __name__ == '__main__':
    sys.exit(main())
