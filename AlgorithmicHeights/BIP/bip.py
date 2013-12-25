import sys

def read_dataset():
    lines = [line.strip() for line in sys.stdin.readlines()]
    n = int(lines[0])
    edge_list = []
    edge_list_array = []
    for line in lines[2:]:
        if not line:
            edge_list_array.append(edge_list)
            edge_list = []
        else:
            edge_list.append(list(map(int, line.split())))
    edge_list_array.append(edge_list)
    return edge_list_array
    
class Graph():
    def __init__(self, edge_list):
        self.edge_list = edge_list
        self.num_nodes = edge_list[0][0]
        self.graph_impl = {}
        self.__build_graph()

    def is_bipartite(self):
        color_array = [None for _ in range(self.num_nodes)]
        # Consider disconnected graphs
        for node in range(1, self.num_nodes + 1):
            if color_array[node - 1] is not None:
                continue
            init_state = [node, True]
            bfs_queue = [init_state]
            while bfs_queue:
                top_node, top_color = bfs_queue[0]
                bfs_queue.pop(0)
                expected_color = color_array[top_node - 1]
                if expected_color is None:
                    color_array[top_node - 1] = top_color
                elif top_color != expected_color:
                    return False
                else:
                    continue
                for next_node in self.graph_impl[top_node]:
                    next_color = not top_color
                    bfs_queue.append([next_node, next_color])
        return True
        
    def __build_graph(self):
        for node in range(self.num_nodes):
            self.graph_impl[node + 1] = []
        for x, y in self.edge_list[1:]:
            self.graph_impl[x].append(y)
            self.graph_impl[y].append(x)
        
def main():
    edge_list_array = read_dataset()
    rv = []
    for edge_list in edge_list_array:
        graph = Graph(edge_list)
        if graph.is_bipartite():
            rv.append(1)
        else:
            rv.append(-1)
    print(' '.join(map(str, rv)))
    
if __name__ == '__main__':
    sys.exit(main())
