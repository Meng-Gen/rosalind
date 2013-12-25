import sys

def read_dataset():
    return [list(map(int, line.strip().split())) for line in sys.stdin.readlines()]
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.infinity = 10**20
        
    def get_shortest_distances(self, source):
        n = self.edge_list[0][0]
        d = [self.infinity for i in range(n + 1)]
        d[source] = 0
        for i in range(1, n - 1):
            for u, v, w in self.edge_list[1:]:
                if d[v] > d[u] + w:
                    d[v] = d[u] + w
        rv = []
        for i in range(1, n + 1):
            if d[i] * 2 > self.infinity:
                rv.append(-1)
            else:
                rv.append(d[i])
        return rv
   
def main():
    edge_list = read_dataset()
    graph = Graph(edge_list)
    print(' '.join(map(str, graph.get_shortest_distances(1))))
    
if __name__ == '__main__':
    sys.exit(main())
