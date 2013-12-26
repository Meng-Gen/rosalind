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
        # Graph stuff
        self.edge_list = edge_list
        self.graph_impl = None
        self.graph_transpose_impl = None
        # SCC graph
        self.components_graph = None 
        # SCCs
        self.components = None
        # Topological sort of SCC graph
        self.sorted_nodes = None 
        # Depth first search stuff
        self.visited_nodes = None
        self.curr_component = None
        self.discovery_time = None
        self.finishing_time = None
        self.time = None
        
        self.__build_graph()
    
    def is_semi_connected(self):
        self.build_component_graph()
        self.topological_sort()
        return self.is_semi_connected_components()
    
    def build_component_graph(self):
        self.depth_first_search()
        self.depth_first_search_transpose()
        self.components_graph = {}
        n = len(self.components)
        for i in range(n):
            self.components_graph[i] = []
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                v_set = set(self.components[j])
                for u in self.components[i]:
                    if set(self.graph_impl[u]) & v_set:
                        self.components_graph[i].append(j)
                        break
    
    def topological_sort(self):
        self.visited_nodes = set()
        self.sorted_nodes = []
        n = len(self.components)
        for u in range(n):
            if u not in self.visited_nodes:
                self.depth_first_search_topological_sort(u)
        self.sorted_nodes = self.sorted_nodes[::-1]
    
    def is_semi_connected_components(self):
        n = len(self.components)
        for i in range(n - 1):
            if i+1 not in self.components_graph[i]:
                return False
        return True
    
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
        self.components = []
        for u, f in sorted(self.finishing_time.items(), key=lambda x: x[1])[::-1]:
            if u not in self.visited_nodes:
                self.curr_component = []
                self.depth_first_search_visit_transpose(u)
                self.components.append(self.curr_component)
                
    def depth_first_search_visit_transpose(self, u):
        self.visited_nodes.add(u)
        self.curr_component.append(u)
        for v in self.graph_transpose_impl[u]:
            if v not in self.visited_nodes:
                self.depth_first_search_visit_transpose(v)
    
    def depth_first_search_topological_sort(self, u):
        if u in self.visited_nodes:
            return
        self.visited_nodes.add(u)
        for v in self.components_graph[u]:
            self.depth_first_search_topological_sort(v)
        self.sorted_nodes.append(u)
    
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
    edge_list_array = read_dataset()
    rv = []
    for edge_list in edge_list_array:
        graph = Graph(edge_list)
        if graph.is_semi_connected():
            rv.append(1)
        else:
            rv.append(-1)
    print(' '.join(map(str, rv)))
    
if __name__ == '__main__':
    sys.exit(main())
