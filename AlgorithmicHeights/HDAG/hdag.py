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
        self.graph_impl = None
        self.sorted_nodes = None 
        self.visited_nodes = None
        self.init_graph()
    
    def get_hamiltonian_path(self):
        self.topological_sort()
        n = self.edge_list[0][0]
        for i in range(n - 1):
            if self.sorted_nodes[i+1] not in self.graph_impl[self.sorted_nodes[i]]:
                return None
        return self.sorted_nodes
    
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
        for v in self.graph_impl[u]:
            self.depth_first_search(v)
        self.sorted_nodes.append(u)
    
    def init_graph(self):
        self.graph_impl = {}
        n = self.edge_list[0][0]
        for u in range(1, n+1):
            self.graph_impl[u] = []
        for u, v in self.edge_list[1:]:
            self.graph_impl[u].append(v)
        
def main():
    edge_list_array = read_dataset()
    for edge_list in edge_list_array:
        graph = Graph(edge_list)
        path = graph.get_hamiltonian_path()
        if path:
            print(1, ' '.join(map(str, path)))
        else:
            print(-1)
    
if __name__ == '__main__':
    sys.exit(main())
