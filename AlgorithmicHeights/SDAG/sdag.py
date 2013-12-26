import sys

def read_dataset():
    return [list(map(int, line.strip().split())) for line in sys.stdin.readlines()]
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.graph_impl = None
        self.visited_nodes = None
        self.sorted_nodes = None
        self.infinity = 10**20
        self.init_graph()
        
    def get_shortest_distances(self):
        self.topological_sort()
        n = self.edge_list[0][0]
        d = [None, 0] + [self.infinity for i in range(n - 1)]
        for u in self.sorted_nodes:
            for v, weight in self.graph_impl[u]:
                if d[v] > d[u] + weight:
                    d[v] = d[u] + weight
        rv = []
        for u in range(1, n + 1):
            if d[u] * 2 > self.infinity:
                rv.append('x')
            else:
                rv.append(d[u])
        return rv
        
    def topological_sort(self):
        self.visited_nodes = set()
        self.sorted_nodes = []
        n = self.edge_list[0][0]
        for u in range(1, n+1):
            if u not in self.visited_nodes:
                self.depth_first_search(u)
        self.sorted_nodes = self.sorted_nodes[::-1]
        
    def depth_first_search(self, u):
        if u in self.visited_nodes:
            return
        self.visited_nodes.add(u)
        for v, weight in self.graph_impl[u]:
            self.depth_first_search(v)
        self.sorted_nodes.append(u)
        
    def init_graph(self):
        self.graph_impl = {}
        n = self.edge_list[0][0]
        for node in range(1, n + 1):
            self.graph_impl[node] = []
        for u, v, weight in self.edge_list[1:]:
            self.graph_impl[u].append([v, weight])        
   
def main():
    edge_list = read_dataset()
    graph = Graph(edge_list)
    print(' '.join(map(str, graph.get_shortest_distances())))
    
if __name__ == '__main__':
    sys.exit(main())
