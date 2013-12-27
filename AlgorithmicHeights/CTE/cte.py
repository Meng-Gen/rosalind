import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    n = int(lines[0])
    edge_list_array = []
    is_begin = True
    num_edge = None
    edge_list = None
    for line in lines[1:]:
        if not line:
            continue
        edge = list(map(int, line.split()))
        if is_begin:
            is_begin = False
            edge_list = [edge]
            num_edge = edge[1]
            continue    
        edge_list.append(edge)
        if len(edge_list) - 1 == num_edge:
            is_begin = True
            edge_list_array.append(edge_list)
    assert(n == len(edge_list_array))
    return edge_list_array
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.infinity = 10**20
        
    def get_shortest_cycle_len(self):
        n = self.edge_list[0][0]
        d = [self.infinity for i in range(n + 1)]
        p, q, pq_weight = self.edge_list[1]
        d[q] = 0
        for i in range(1, n):
            for u, v, w in self.edge_list[2:]:
                if d[v] > d[u] + w:
                    d[v] = d[u] + w
        if d[p] * 2 > self.infinity:
            return -1
        else:
            return d[p] + pq_weight
   
def main():
    edge_list_array = read_dataset()
    rv = []
    for edge_list in edge_list_array:
        graph = Graph(edge_list)
        rv.append(graph.get_shortest_cycle_len())
    print(' '.join(map(str, rv)))
    
if __name__ == '__main__':
    sys.exit(main())
