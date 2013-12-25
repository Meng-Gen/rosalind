import sys

def read_dataset():
    return [list(map(int, line.strip().split())) for line in sys.stdin.readlines()]
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.graph_impl = None
        self.visited_nodes = None
        self.sorted_nodes = None
        self.__build_graph()
    
    def topological_sort(self):
        self.visited_nodes = set()
        self.sorted_nodes = []
        n = self.edge_list[0][0]
        for node in range(1, n+1):
            if node not in self.visited_nodes:
                self.depth_first_search(node)
        return self.sorted_nodes[::-1]
        
    def depth_first_search(self, top_node):
        if top_node in self.visited_nodes:
            return
        self.visited_nodes.add(top_node)
        for next_node in self.graph_impl[top_node]:
            self.depth_first_search(next_node)
        self.sorted_nodes.append(top_node)
        
    def __build_graph(self):
        self.graph_impl = {}
        n = self.edge_list[0][0]
        for node in range(n):
            self.graph_impl[node + 1] = []
        for x, y in self.edge_list[1:]:
            self.graph_impl[x].append(y)
        
def main():
    edge_list = read_dataset()
    graph = Graph(edge_list)
    print(' '.join(map(str, graph.topological_sort())))
    
if __name__ == '__main__':
    sys.exit(main())
