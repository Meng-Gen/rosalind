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
        
    def has_negative_weight_cycle(self):
        n = self.edge_list[0][0]
        d = [self.infinity for i in range(n + 1)]
        d[1] = 0
        for i in range(1, n - 1):
            for u, v, w in self.edge_list[1:]:
                if d[v] > d[u] + w:
                    d[v] = d[u] + w
        for u, v, w in self.edge_list[1:]:
            if d[v] > d[u] + w:
                return True
        return False
   
def main():
    edge_list_array = read_dataset()
    rv = []
    for edge_list in edge_list_array:
        graph = Graph(edge_list)
        if graph.has_negative_weight_cycle():
            rv.append(1)
        else:
            rv.append(-1)
    print(' '.join(map(str, rv)))
    
if __name__ == '__main__':
    sys.exit(main())
