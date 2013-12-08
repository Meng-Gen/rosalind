import sys

def dna_rc(dna):
    rc = dna.replace("A", "a").replace("T", "A").replace("a", "T")
    rc = rc.replace("C", "c").replace("G", "C").replace("c", "G")
    return rc[::-1]

def get_de_bruijn_set(dna_records, k):
    dna_rc_records = [dna_rc(dna) for dna in dna_records]
    n = len(dna_records[0])
    de_bruijn_set = set()
    for dna in set(dna_records + dna_rc_records):
        for i in range(n - k + 1):
            de_bruijn_set.add(dna[i:i+k])
    return de_bruijn_set
    
def get_de_bruijn_graph(de_bruijn_set):
    graph = {}
    for dna in de_bruijn_set:
        begin, end = dna[:-1], dna[1:]
        graph[begin] = end
    return graph

def has_cycles(graph):    
    node_degree = {}
    for node in graph:
        if node not in node_degree:
            node_degree[node] = [0, 0] # indegree, outdegree
        if graph[node] not in node_degree:
            node_degree[graph[node]] = [0, 0]
        node_degree[node][0] += 1
        node_degree[graph[node]][1] += 1
    for node in node_degree:
        if node_degree[node][0] != 1 or node_degree[node][1] != 1:
            return False
    return True
    
def equals_to_cycles(graph, expected_count):
    if not has_cycles(graph):
        return False
    actual_count = 0
    visited_nodes = set()
    for node in graph:
        if node in visited_nodes:
            continue
        curr_node = node
        while curr_node not in visited_nodes:
            visited_nodes.add(curr_node)
            curr_node = graph[curr_node]
        actual_count += 1
    return actual_count == expected_count

def make_circular_string(graph, k):
    init_node = get_init_node(graph)
    curr_node = graph[init_node]
    circular_string = init_node
    while curr_node != init_node:
        circular_string += curr_node[-1]
        curr_node = graph[curr_node]       
    return circular_string[:-k+1]
    
def get_init_node(graph):
    for node in graph:
        return node   

def main():
    dna_records = [dna.strip() for dna in sys.stdin.readlines()]
    n = len(dna_records[0])
    max_de_bruijn_graph, max_k = None, None
    for k in range(n, 1, -1):
        de_bruijn_set = get_de_bruijn_set(dna_records, k)
        de_bruijn_graph = get_de_bruijn_graph(de_bruijn_set)
        if equals_to_cycles(de_bruijn_graph, 2):
            max_de_bruijn_graph, max_k = de_bruijn_graph, k
            break
    print(make_circular_string(max_de_bruijn_graph, max_k - 1))

if __name__ == '__main__':
    sys.exit(main())
