import sys

def is_adjacent(s, t, overlap_length):
    return s != t and s[-overlap_length:] == t[:overlap_length]

def dna_rc(dna):
    rc = dna.replace("A", "a").replace("T", "A").replace("a", "T")
    rc = rc.replace("C", "c").replace("G", "C").replace("c", "G")
    return rc[::-1]

def get_edges(dna_records):
    edges = {}
    dna_rc_records = [dna_rc(dna) for dna in dna_records]
    for dna in set(dna_records + dna_rc_records):
        begin, end = dna[:-1], dna[1:]
        edges[begin] = end
    return edges
    
def make_circular_string(edges, k):
    init_node = get_init_node(edges)
    curr_node = edges[init_node]
    circular_string = init_node
    while curr_node != init_node:
        circular_string += curr_node[-1]
        curr_node = edges[curr_node]       
    return circular_string[:-k+1]
    
def get_init_node(edges):
    for node in edges:
        return node   
    
def main():
    dna_records = [dna.strip() for dna in sys.stdin.readlines()]
    edges = get_edges(dna_records)
    k = len(dna_records[0])
    print(make_circular_string(edges, k-1))

if __name__ == '__main__':
    sys.exit(main())
