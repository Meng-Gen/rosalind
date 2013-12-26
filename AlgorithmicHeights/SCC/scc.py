import sys

def read_dataset():
    return [list(map(int, line.strip().split())) for line in sys.stdin.readlines()]
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.graph_impl = None
        self.graph_transpose_impl = None
        self.visited_nodes = None
        self.discovery_time = None
        self.finishing_time = None
        self.time = None
        self.has_cycle = None
        self.__build_graph()
    
    def num_strongly_connected_components(self):
        self.depth_first_search()
        return self.depth_first_search_transpose()
    
    def depth_first_search(self):
        self.visited_nodes = set()
        self.discovery_time = {}
        self.finishing_time = {}
        self.time = 0
        n = self.edge_list[0][0]
        for u in range(1, n+1):
            self.discovery_time[u] = None
            self.finishing_time[u] = None
        for u in range(1, n+1):
            if u not in self.visited_nodes:
                self.depth_first_search_visit(u)
   
    def depth_first_search_visit(self, u):
        self.time += 1
        self.discovery_time[u] = self.time
        self.visited_nodes.add(u)
        for v in self.graph_impl[u]:
            if v not in self.visited_nodes:
                self.depth_first_search_visit(v)
        self.time += 1
        self.finishing_time[u] = self.time
        
    def depth_first_search_transpose(self):
        self.visited_nodes = set()
        count = 0
        for u, f in sorted(self.finishing_time.items(), key=lambda x: x[1])[::-1]:
            if u not in self.visited_nodes:
                count += 1
                self.depth_first_search_visit_transpose(u)
        return count
                
    def depth_first_search_visit_transpose(self, u):
        self.visited_nodes.add(u)
        for v in self.graph_transpose_impl[u]:
            if v not in self.visited_nodes:
                self.depth_first_search_visit_transpose(v)
        
    def __build_graph(self):
        self.graph_impl = {}
        self.graph_transpose_impl = {}
        n = self.edge_list[0][0]
        for node in range(n):
            self.graph_impl[node + 1] = []
            self.graph_transpose_impl[node + 1] = []
        for x, y in self.edge_list[1:]:
            self.graph_impl[x].append(y)
            self.graph_transpose_impl[y].append(x)
        
def main():
    edge_list = read_dataset()
    graph = Graph(edge_list)
    print(graph.num_strongly_connected_components())
    
if __name__ == '__main__':
    sys.exit(main())
