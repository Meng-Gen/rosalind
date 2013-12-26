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
        self.visited_nodes = None
        self.has_cycle = None
        self.__build_graph()
    
    def is_acyclic(self):
        self.visited_nodes = {}
        n = self.edge_list[0][0]
        for node in range(1, n+1):
            if node not in self.visited_nodes:
                self.depth_first_search(node)
        return not self.has_cycle
        
    def depth_first_search(self, top_node):
        if top_node in self.visited_nodes:
            if self.visited_nodes[top_node] == 'BackEdge':
                self.has_cycle = True
            return
                
        self.visited_nodes[top_node] = 'BackEdge'
        for next_node in self.graph_impl[top_node]:
            self.depth_first_search(next_node)
        self.visited_nodes[top_node] = 'NotBackEdge'
        
    def __build_graph(self):
        self.graph_impl = {}
        n = self.edge_list[0][0]
        for node in range(n):
            self.graph_impl[node + 1] = []
        for x, y in self.edge_list[1:]:
            self.graph_impl[x].append(y)
        
def main():
    edge_list_array = read_dataset()
    rv = []
    for edge_list in edge_list_array:
        graph = Graph(edge_list)
        if graph.is_acyclic():
            rv.append(1)
        else:
            rv.append(-1)
    print(' '.join(map(str, rv)))
    
if __name__ == '__main__':
    sys.exit(main())
